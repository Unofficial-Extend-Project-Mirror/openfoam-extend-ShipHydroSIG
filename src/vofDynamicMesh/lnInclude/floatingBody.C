/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "floatingBody.H"
#include "transform.H"
#include "transformField.H"
#include "fixedValueTetPolyPatchFields.H"

#include "incompressible/incompressibleTwoPhaseMixture/twoPhaseMixture.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::floatingBody::check() const
{
    forAll (hullPatches_, patchI)
    {
        if (!hullPatches_[patchI].active())
        {
            FatalErrorIn("void floatingBody::check() const")
                << "Hull patch named " << hullPatches_[patchI].name()
                << " not found.  Available patches: "
                << mesh_.boundaryMesh().names()
                << abort(FatalError);
        }
    }
}


Foam::tmp<Foam::volSymmTensorField> Foam::floatingBody::devRhoReff() const
{
    if (mesh_.foundObject<incompressible::RASModel>("RASProperties"))
    {
        const incompressible::RASModel& ras =
            mesh_.lookupObject<incompressible::RASModel>("RASProperties");

        return ras.devReff();
    }
    else if (mesh_.foundObject<incompressible::LESModel>("LESProperties"))
    {
        const incompressible::LESModel& les =
            mesh_.lookupObject<incompressible::LESModel>("LESProperties");

        return les.devBeff();
    }
    else if (mesh_.foundObject<twoPhaseMixture>("transportProperties"))
    {
        const twoPhaseMixture& twoPhaseProperties =
            mesh_.lookupObject<twoPhaseMixture>("transportProperties");

        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        return twoPhaseProperties.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorIn("floatingBody::devRhoReff()")
            << "No valid model for viscous stress calculation."
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::vector Foam::floatingBody::dragForce() const
{
    // Note: Return zero if there is no flow solution
    // HR, 6/Jan/2010
    if
    (
        !mesh_.foundObject<volVectorField>("U")
        || !mesh_.foundObject<volScalarField>("p")
    )
    {
        WarningIn
        (
            "floatingBody::dragForce() const"
        )   << "No flow solution found. Returning zero."
            << endl;

        return vector::zero;
    }

    const volScalarField& alpha = mesh_.lookupObject<volScalarField>("alpha1");
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    const surfaceVectorField::GeometricBoundaryField& Sfb =
        mesh_.Sf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
    const volSymmTensorField::GeometricBoundaryField& devRhoReffb
        = tdevRhoReff().boundaryField();

    vector pressureForce = vector::zero;
    vector viscousForce = vector::zero;

    forAll (hullPatches_, patchI)
    {
        const label curPatch = hullPatches_[patchI].index();

        // Note: calculating forces on the wetted surface only
        // HJ, 9/Mar/2009

        pressureForce +=
            gSum
            (
                alpha.boundaryField()[curPatch]*
                p.boundaryField()[curPatch]*Sfb[curPatch]
            );

        viscousForce +=
            gSum
            (
                alpha.boundaryField()[curPatch]*
                (Sfb[curPatch] & devRhoReffb[curPatch])
            );
    }

    vector totalForce = pressureForce + viscousForce;

    Info<< "Pressure force = " << pressureForce << nl
        << "Viscous force = " << viscousForce << nl
        << "Total force = " << totalForce << endl;

    return totalForce;
}


Foam::vector Foam::floatingBody::dragMoment() const
{
    // Note: Return zero if there is no flow solution
    // HR, 6/Jan/2010
    if
    (
        !mesh_.foundObject<volVectorField>("U")
        || !mesh_.foundObject<volScalarField>("p")
    )
    {
        WarningIn
        (
            "floatingBody::dragForce() const"
        )   << "No flow solution found. Returning zero."
            << endl;

        return vector::zero;
    }

    const volScalarField& alpha = mesh_.lookupObject<volScalarField>("alpha1");
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    const surfaceVectorField::GeometricBoundaryField& Sfb =
        mesh_.Sf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
    const volSymmTensorField::GeometricBoundaryField& devRhoReffb
        = tdevRhoReff().boundaryField();

    vector pressureMoment = vector::zero;
    vector viscousMoment = vector::zero;

    forAll (hullPatches_, patchI)
    {
        // Get index of current patch
        const label curPatch = hullPatches_[patchI].index();

        // Note: calculating forces on the wetted surface only
        // HJ, 9/Mar/2009

        pressureMoment +=
            gSum
            (
                (mesh_.Cf().boundaryField()[curPatch] - sixDOF_.X().value())
              ^ (
                    alpha.boundaryField()[curPatch]*
                    p.boundaryField()[curPatch]*Sfb[curPatch]
                )
            );

        viscousMoment +=
            gSum
            (
                (mesh_.Cf().boundaryField()[curPatch] - sixDOF_.X().value())
              ^ (
                    alpha.boundaryField()[curPatch]*
                    (Sfb[curPatch] & devRhoReffb[curPatch])
                )
            );
    }

    vector totalMoment = pressureMoment + viscousMoment;

    Info<< "Pressure moment = " << pressureMoment << nl
        << "Viscous moment = " << viscousMoment << nl
        << "Total moment = " << totalMoment << endl;

    return totalMoment;
}


void Foam::floatingBody::fixTranslation(vector& v) const
{
    if (fixedSurge_)
    {
        v.x() = 0;
    }

    if (fixedSway_)
    {
        v.y() = 0;
    }

    if (fixedHeave_)
    {
        v.z() = 0;
    }
}


void Foam::floatingBody::fixRotation(vector& rot) const
{
    if (fixedRoll_)
    {
        rot.x() = 0;
    }

    if (fixedPitch_)
    {
        rot.y() = 0;
    }

    if (fixedYaw_)
    {
        rot.z() = 0;
    }

}


Foam::vector Foam::floatingBody::rotCentreVelocity() const
{
    return
    (
        // New position of rotation centroid in absolute coordinate system
        transform(sixDOF_.toAbsolute(), centreOfSlider_)
        // Motion of centre of gravity
      + sixDOF_.X().value()
        // Old position of rotation centroid
      - rotCentre_
    )/mesh_.time().deltaT().value();
}


Foam::tmp<Foam::vectorField> Foam::floatingBody::boundaryVelocity
(
    const vectorField& boundaryPoints
) const
{
    return
        // Translation of centroid
        rotCentreVelocity()
        // Rotation around centroid
      + (
            // New points,scaled around centroid
            transform
            (
                sixDOF_.rotIncrementTensor().T(),
                boundaryPoints - rotCentre_
            )
            // Old points, scaled around centroid
          - (boundaryPoints - rotCentre_)
        )/mesh_.time().deltaT().value();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floatingBody::floatingBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),

    hullPatches_(),
    hullSlider_(word(dict.lookup("hullSlider")), mesh_.boundaryMesh()),
    hullSliderZone_(word(dict.lookup("hullSliderZone")), mesh_.faceZones()),
    fixedSlider_(word(dict.lookup("fixedSlider")), mesh_.boundaryMesh()),
    fixedSliderZone_(word(dict.lookup("fixedSliderZone")), mesh_.faceZones()),
    centreOfSlider_(dict.lookup("centreOfSlider")),

    addPropulsionForce_(dict.lookup("addPropulsionForce")),
    propulsionDirection_(dict.lookup("propulsionDirection")),
    propulsionForceCentroid_(dict.lookup("propulsionForceCentroid")),

    sixDOF_
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    ),
    solver_
    (
        ODESolver::New
        (
            dict.lookup("odeSolver"),
            sixDOF_
        )
    ),
    epsilon_(readScalar(dict.lookup("epsilon"))),

    fixedSurge_(dict.lookup("fixedSurge")),
    fixedSway_(dict.lookup("fixedSway")),
    fixedHeave_(dict.lookup("fixedHeave")),
    fixedRoll_(dict.lookup("fixedRoll")),
    fixedPitch_(dict.lookup("fixedPitch")),
    fixedYaw_(dict.lookup("fixedYaw")),

    steadyState_(dict.lookup("steadyState")),
    curTimeIndex_(-1),
    rotCentre_(sixDOF_.X().value()),
    oldForce_(vector::zero),
    oldMoment_(vector::zero),
    // Dealing with writing a motion file in parallel.  HJ, 14/Aug/2008
    of_(mesh_.time().caseConstant()/".."/name_ + ".dat")
{
    // Read hull patch names
    wordList hullPatchNames(dict.lookup("hullPatches"));

    hullPatches_.setSize(hullPatchNames.size());

    forAll (hullPatchNames, patchI)
    {
        hullPatches_.set
        (
            patchI,
            new polyPatchID
            (
                hullPatchNames[patchI],
                mesh_.boundaryMesh()
            )
        );
    }

    // Rescale propulsion direction vector
    if (mag(propulsionDirection_) < SMALL)
    {
        FatalErrorIn
        (
            "floatingBody::floatingBody(const word& name, const fvMesh& mesh,"
            "const dictionary& dict)"
        )   << " Incorrect propulsion direction: " << propulsionDirection_
            << abort(FatalError);
    }
    else
    {
        propulsionDirection_ /= mag(propulsionDirection_);
    }

    // Block (initial) motion components as requested

    // Translation
    fixTranslation(sixDOF_.U().value());

    // Rotation
    fixRotation(sixDOF_.omega().value());

    check();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::floatingBody::setMotion(tetPointVectorField& motionU)
{
    // Calculate forces and moments
    vector curForce = dragForce();
    vector curMoment = dragMoment();

    if (addPropulsionForce_)
    {
        // Here, curForce = drag.  HJ, 11/Mar/2009

        const vector propulsionForce =
            - propulsionDirection_*(propulsionDirection_ & curForce);

        const vector propulsionMoment =
            (
                (propulsionForceCentroid_ - sixDOF_.X().value())
              ^ propulsionForce
            );

        Info<< "Propulsion force = " << propulsionForce
            << " Propulsion moment = " << propulsionMoment << endl;

        curForce += propulsionForce;
        curMoment += propulsionMoment;
    }

    const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");

    // Grab initial force for the first time-step
    if (curTimeIndex_ == -1)
    {
        oldForce_ = curForce;
        oldMoment_ = curMoment;
        curTimeIndex_ = mesh_.time().timeIndex();
    }

    // Set force and moment
    sixDOF_.force().value() =
        0.5*(curForce + oldForce_)
      + sixDOF_.mass().value()*g.value();

    sixDOF_.moment().value() =
        0.5*(curMoment + oldMoment_);

    // Manage multiple calls to mesh motion within a single time-step
    if (curTimeIndex_ < mesh_.time().timeIndex())
    {
        oldForce_ = curForce;
        oldMoment_ = curMoment;
        curTimeIndex_ = mesh_.time().timeIndex();
    }

    // Manage steady-state: reset velocities to zero
    if (steadyState_)
    {
        sixDOF_.U().value() = vector::zero;
        sixDOF_.omega().value() = vector::zero;
    }
    else
    {
        // Block motion components as requested

        // Translation
        fixTranslation(sixDOF_.force().value());

        // Rotation
        fixRotation(sixDOF_.moment().value());
    }

    // Grab centre of rotation before the ODE solution
    rotCentre_ = sixDOF_.X().value()
        + transform(sixDOF_.toAbsolute(), centreOfSlider_);

    // Solve for motion
    solver_->solve
    (
        mesh_.time().value(),
        mesh_.time().value() +  mesh_.time().deltaT().value(),
        epsilon_,
        mesh_.time().deltaT().value()
    );

    Info<< "body: " << name() << nl
        << "f = " << sixDOF_.force().value() << tab
        << "mg = " << sixDOF_.mass().value()*g.value() << tab
        << "M = " << sixDOF_.moment().value() << nl
        << "x = " << sixDOF_.X().value() << tab
        << "u = " << sixDOF_.Uaverage().value() << nl
        << "rot vector = " << sixDOF_.rotVector() << tab
        << "rot angle = " << sixDOF_.rotAngle().value() << tab
        << "omega = " << sixDOF_.omegaAverageAbsolute().value() << endl;

    // Set boundary motion on the hull patches
    forAll (hullPatches_, patchI)
    {
        fixedValueTetPolyPatchVectorField& motionUHull =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[hullPatches_[patchI].index()]
            );

        // Set motion based on current hull geometry
        motionUHull ==
            boundaryVelocity
            (
                motionUHull.patch().localPoints()
            );
    }

    // Set motion on slider surfaces
    if (hullSlider_.active())
    {
        fixedValueTetPolyPatchVectorField& motionUHullSlider =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[hullSlider_.index()]
            );

        motionUHullSlider ==
            boundaryVelocity
            (
                motionUHullSlider.patch().localPoints()
            );
    }

    if (fixedSlider_.active())
    {
        fixedValueTetPolyPatchVectorField& motionUFixedSlider =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[fixedSlider_.index()]
            );

        motionUFixedSlider == rotCentreVelocity();
    }


    // Report position and rotation.  Write only if master
    if (Pstream::master())
    {
        of_ << mesh_.time().value() << tab;

        // Translation
        if (!fixedSurge_)
        {
            of_ << sixDOF_.X().value().x() << tab;
        }

        if (!fixedSway_)
        {
            of_ << sixDOF_.X().value().y() << tab;
        }

        if (!fixedHeave_)
        {
            of_ << sixDOF_.X().value().z() << tab;
        }

        // Rotation
        if (!fixedRoll_ || !fixedPitch_ || !fixedYaw_)
        {
            of_ << sixDOF_.rotAngle().value() << tab;
        }

        of_ << endl;
    }
}


void Foam::floatingBody::correctMotion(vectorField& newAllPoints) const
{
    // Adjust motion of slider zone points
    if (hullSliderZone_.active())
    {
        // Hull slider zone takes translation and rotation
        const pointField newP =
            boundaryVelocity
            (
                mesh_.faceZones()[hullSliderZone_.index()]().localPoints()
            );

        const labelList& addr =
            mesh_.faceZones()[hullSliderZone_.index()]().meshPoints();

        forAll (addr, i)
        {
            if (addr[i] >= mesh_.nPoints())
            {
                newAllPoints[addr[i]] += newP[i]*mesh_.time().deltaT().value();
            }
        }
    }

    if (fixedSliderZone_.active())
    {
        // Fixed slider zone takes translation only
        const vector rcv = rotCentreVelocity();

        const labelList& addr =
            mesh_.faceZones()[fixedSliderZone_.index()]().meshPoints();

        forAll (addr, i)
        {
            if (addr[i] >= mesh_.nPoints())
            {
                newAllPoints[addr[i]] += rcv*mesh_.time().deltaT().value();
            }
        }
    }
}


// ************************************************************************* //
