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

#include "threeDofMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "surfaceFields.H"
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
    defineTypeNameAndDebug(threeDofMotion, 0);

    addToRunTimeSelectionTable(dynamicFvMesh, threeDofMotion, IOobject);
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField> Foam::threeDofMotion::devRhoReff() const
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


Foam::vector Foam::threeDofMotion::hullForce() const
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

Foam::threeDofMotion::threeDofMotion(const IOobject& io)
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
    motionPtr_(motionSolver::New(*this)),
    hullPatch_(word(dict_.lookup("hullPatchName")), boundaryMesh()),
    curTimeIndex_(-1),
    oldForce_(vector::zero)
{
    // Check that the hull patch has been found
    if (!hullPatch_.active())
    {
        FatalErrorIn("threeDofMotion::threeDofMotion(const IOobject& io)")
            << "Cannot find hull patch "<< hullPatch_.name() << nl
            << "Valid patches are: " << boundaryMesh().names()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::threeDofMotion::~threeDofMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::threeDofMotion::update()
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
            this->objectRegistry::lookupObject<tetPointVectorField>("motionU")
        );

    fixedValueTetPolyPatchVectorField& motionUHullPatch =
        refCast<fixedValueTetPolyPatchVectorField>
        (
            motionU.boundaryField()[hullPatch_.index()]
        );

    Info << "Boundary velocity = " << equation_.U().value() << endl;

    motionUHullPatch == equation_.U().value();

    fvMesh::movePoints(motionPtr_->newPoints());

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //
