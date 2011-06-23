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

#include "sixDofMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "cellSet.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "fixedValueTetPolyPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDofMotion, 0);

    addToRunTimeSelectionTable(dynamicFvMesh, sixDofMotion, IOobject);
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDofMotion::sixDofMotion(const IOobject& io)
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
    bodies_(*this, dict_.lookup("bodies")),
    motionPtr_(motionSolver::New(*this))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDofMotion::~sixDofMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sixDofMotion::update()
{
    // Set the motion onto the motion patch
    tetPointVectorField& motionU =
        const_cast<tetPointVectorField&>
        (
            this->objectRegistry::lookupObject<tetPointVectorField>("motionU")
        );

    // Set motion from bodies
    bodies_.setMotion(motionU);

    // Calculate the motion and move points
    fvMesh::movePoints(motionPtr_->newPoints());

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //
