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

\*---------------------------------------------------------------------------*/

#include "threeDofMotionSubset.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "cellSet.H"
#include "uniformDimensionedFields.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "fixedValueTetPolyPatchFields.H"

#include "incompressible/incompressibleTwoPhaseMixture/twoPhaseMixture.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(threeDofMotionSubset, 0);

    addToRunTimeSelectionTable(dynamicFvMesh, threeDofMotionSubset, IOobject);
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::threeDofMotionSubset::devRhoReff() const
{
    if (this->foundObject<incompressible::RASModel>("RASProperties"))
    {
        const incompressible::RASModel& ras =
            this->lookupObject<incompressible::RASModel>("RASProperties");

        return ras.devReff();
    }
    else if (this->foundObject<incompressible::LESModel>("LESProperties"))
    {
        const incompressible::LESModel& les =
            this->lookupObject<incompressible::LESModel>("LESProperties");

        return les.devBeff();
    }
    else if (this->foundObject<twoPhaseMixture>("transportProperties"))
    {
        const twoPhaseMixture& twoPhaseProperties =
            this->lookupObject<twoPhaseMixture>("transportProperties");

        const volVectorField& U = this->lookupObject<volVectorField>("U");

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


Foam::vector Foam::threeDofMotionSubset::hullForce() const
{
    const fvMesh& mesh = *this;

    const volScalarField& alpha = mesh.lookupObject<volScalarField>("alpha1");
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
    const volSymmTensorField::GeometricBoundaryField& devRhoReffb
        = tdevRhoReff().boundaryField();

    // Get index of current patch
    const label curPatch = hullPatch_.index();

    vector pressureForce =
        gSum
        (
            alpha.boundaryField()[curPatch]*
            p.boundaryField()[curPatch]*
            mesh.Sf().boundaryField()[curPatch]
        );

    vector viscousForce =
        gSum
        (
            alpha.boundaryField()[curPatch]*
            (mesh.Sf().boundaryField()[curPatch] & devRhoReffb[curPatch])
        );

    vector totalForce = pressureForce + viscousForce;

    Info<< "Pressure force = " << pressureForce << nl
        << "Viscous force = " << viscousForce << nl
        << "Total force = " << totalForce << endl;

    return totalForce;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threeDofMotionSubset::threeDofMotionSubset(const IOobject& io)
:
    dynamicFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    equation_
    (
        IOobject
        (
            dict_.lookup("name"),
            time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    ),
    solver_
    (
        ODESolver::New
        (
            dict_.lookup("solver"),
            equation_
        )
    ),
    epsilon_(readScalar(dict_.lookup("epsilon"))),
    subsetMesh_
    (
        IOobject
        (
            "motion",
            io.time().constant(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    motionPtr_(NULL),
    hullPatch_(word(dict_.lookup("hullPatchName")), boundaryMesh()),
    motionPatch_(dict_.lookup("motionPatchName")),
    centroid_(dict_.lookup("centreOfGravity")),
    curTimeIndex_(-1),
    oldForce_(vector::zero)
{
    // Check that the hull patch has been found
    if (!hullPatch_.active())
    {
        FatalErrorIn
        (
            "threeDofMotionSubset::threeDofMotionSubset(const IOobject& io)"
        )   << "Cannot find hull patch "<< hullPatch_.name() << nl
            << "Valid patches are: " << boundaryMesh().names()
            << abort(FatalError);
    }

    // Create subset
    word setName = dict_.lookup("set");

    cellSet currentSet
    (
        IOobject
        (
            setName,
            io.time().constant(),
            polyMesh::meshSubDir/"sets",
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    subsetMesh_.setLargeCellSubset(currentSet, -1);
    subsetMesh_.subMesh().write();

    // Create motion solver on the subset
    motionPtr_ = motionSolver::New(subsetMesh_.subMesh());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::threeDofMotionSubset::~threeDofMotionSubset()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::threeDofMotionSubset::update()
{
    // Grab force
    vector curForce = hullForce();

    const uniformDimensionedVectorField& g =
        this->lookupObject<uniformDimensionedVectorField>("g");

    equation_.force().value() = 0.5*(curForce + oldForce_)
        + equation_.mass().value()*g.value();

    if (curTimeIndex_ == -1 || curTimeIndex_ < time().timeIndex())
    {
        Info << "Storing current force as old force" << endl;
        oldForce_ = curForce;
        curTimeIndex_ = time().timeIndex();
    }

    solver_->solve
    (
        time().value(),
        time().value() + time().deltaT().value(),
        epsilon_,
        time().deltaT().value()
    );

    // Set the motion onto the motion patch
    tetPointVectorField& motionU =
        const_cast<tetPointVectorField&>
        (
            subsetMesh_.subMesh().objectRegistry::
            lookupObject<tetPointVectorField>
            (
                "motionU"
            )
        );

    // Find moving patch index on the subset mesh
    polyPatchID movingPatchID
    (
        motionPatch_,
        subsetMesh_.subMesh().boundaryMesh()
    );

    if (!movingPatchID.active())
    {
        FatalErrorIn("bool threeDofMotionSubset::update()")
            << "Problem finding patch " << hullPatch_.name()
            << " on subset" << abort(FatalError);
    }

    fixedValueTetPolyPatchVectorField& motionUHullPatch =
        refCast<fixedValueTetPolyPatchVectorField>
        (
            motionU.boundaryField()[movingPatchID.index()]
        );

    Info << "Boundary velocity = " << equation_.U().value() << endl;

    motionUHullPatch == equation_.U().value();

    // Get the points from the moving part
    pointField subsetPoints = motionPtr_->newPoints();

    //- Get copy of mesh points
    pointField p = points() + equation_.U().value()*time().deltaT().value();

    //- Map the moving part
    const labelList& subsetPointAddr = subsetMesh_.pointMap();

    forAll (subsetPoints, subsetI)
    {
        p[subsetPointAddr[subsetI]] = subsetPoints[subsetI];
    }

    // Move subset points
    subsetMesh_.subMesh().movePoints(subsetPoints);

    // Move mesh points
    fvMesh::movePoints(p);

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //
