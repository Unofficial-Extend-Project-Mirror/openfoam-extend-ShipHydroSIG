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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "surfaceFields.H"
#include "rasVofForceAndTorqueFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

#include "incompressible/incompressibleTwoPhaseMixture/twoPhaseMixture.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rasVofForceAndTorqueFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        rasVofForceAndTorqueFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::rasVofForceAndTorqueFunctionObject::devRhoReff() const
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    if (mesh.foundObject<incompressible::RASModel>("RASProperties"))
    {
        const incompressible::RASModel& ras =
            mesh.lookupObject<incompressible::RASModel>("RASProperties");

        return ras.devReff();
    }
    else if (mesh.foundObject<incompressible::LESModel>("LESProperties"))
    {
        const incompressible::LESModel& les =
            mesh.lookupObject<incompressible::LESModel>("LESProperties");

        return les.devBeff();
    }
    else if (mesh.foundObject<twoPhaseMixture>("transportProperties"))
    {
        const twoPhaseMixture& twoPhaseProperties =
            mesh.lookupObject<twoPhaseMixture>("transportProperties");

        const volVectorField& U = mesh.lookupObject<volVectorField>("U");

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rasVofForceAndTorqueFunctionObject::
rasVofForceAndTorqueFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    patchNames_(dict.lookup("patches")),
    origin_(dict.lookup("origin")),
    of_(time_.path()/word(dict.lookup("file")))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating rasVofForceAndTorqueFunctionObject for patches "
        << patchNames_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rasVofForceAndTorqueFunctionObject::start()
{
    return true;
}


bool Foam::rasVofForceAndTorqueFunctionObject::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volScalarField& alpha = mesh.lookupObject<volScalarField>("alpha1");
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    const surfaceVectorField::GeometricBoundaryField& Sfb =
        mesh.Sf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
    const volSymmTensorField::GeometricBoundaryField& devRhoReffb
        = tdevRhoReff().boundaryField();

    vectorField pressureForces(patchNames_.size(), vector::zero);
    vectorField pressureMoments(patchNames_.size(), vector::zero);

    vectorField viscousForces(patchNames_.size(), vector::zero);
    vectorField viscousMoments(patchNames_.size(), vector::zero);

    forAll (patchNames_, patchI)
    {
        const label curPatch =
            mesh.boundaryMesh().findPatchID(patchNames_[patchI]);

        if (curPatch > -1)
        {
            // Patch found

            // Pressure component
            pressureForces[patchI] =
                gSum
                (
                    alpha.boundaryField()[curPatch]*
                    p.boundaryField()[curPatch]*Sfb[curPatch]
                );

            pressureMoments[patchI] =
                gSum
                (
                    (mesh.Cf().boundaryField()[curPatch] - origin_)
                  ^ (
                      alpha.boundaryField()[curPatch]*
                        p.boundaryField()[curPatch]*Sfb[curPatch]
                    )
                );

            // Viscous component
            viscousForces[patchI] =
                gSum
                (
                    alpha.boundaryField()[curPatch]*
                    (Sfb[curPatch] & devRhoReffb[curPatch])
                );

            viscousMoments[patchI] =
                gSum
                (
                    (mesh.Cf().boundaryField()[curPatch] - origin_)
                  ^ (
                        alpha.boundaryField()[curPatch]*
                        (Sfb[curPatch] & devRhoReffb[curPatch])
                    )
                );
        }
        else
        {
            FatalErrorIn
            (
                "bool rasVofForceAndTorqueFunctionObject::execute()"
            )   << "Patch named " << patchNames_[patchI] << " not found.  "
                << "Available patches are: "
                << mesh.boundaryMesh().names()
                << abort(FatalError);
        }
    }

    vectorField totalForces = pressureForces + viscousForces;
    vectorField totalMoments = pressureMoments + viscousMoments;

    Info<< "Pressure forces = " << pressureForces << nl
        << "Pressure moments = " << pressureMoments << nl
        << "Viscous forces = " << viscousForces << nl
        << "Viscous moments = " << viscousMoments << nl
        << "Total forces = " << totalForces << nl
        << "Total moments = " << totalMoments << endl;

    of_ << time_.value() << tab
        << sum(totalForces).x() << tab
        << sum(totalForces).y() << tab
        << sum(totalForces).z() << tab
        << sum(totalMoments).x() << tab
        << sum(totalMoments).y() << tab
        << sum(totalMoments).z() << endl;

    if (patchNames_.size() > 1)
    {
        Info<< "sum pressure forces = " << sum(pressureForces) << nl
            << "sum pressure moments = " << sum(pressureMoments) << nl
            << "sum viscous forces = " << sum(viscousForces) << nl
            << "sum viscous moments = " << sum(viscousMoments) << nl
            << "sum total forces = " << sum(totalForces) << nl
            << "sum total moments = " << sum(totalMoments) << nl << endl;
    }
    else
    {
        Info << nl << endl;
    }

    return true;
}


bool Foam::rasVofForceAndTorqueFunctionObject::read(const dictionary& dict)
{
    patchNames_ = wordList(dict.lookup("patches"));
    origin_ = point(dict.lookup("origin"));

    return false;
}

// ************************************************************************* //
