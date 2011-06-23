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

#include "TaitBryanCoordinateRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TaitBryanCoordinateRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        TaitBryanCoordinateRotation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::TaitBryanCoordinateRotation::calcTransform
(
    const scalar phiAngle,
    const scalar thetaAngle,
    const scalar psiAngle,
    const bool inDegrees
)
{
    scalar alpha = phiAngle;  // Yaw
    scalar beta = thetaAngle; // Pitch
    scalar gamma = psiAngle;  // Roll

    if (inDegrees)
    {
        alpha *= mathematicalConstant::pi/180.0;
        beta *= mathematicalConstant::pi/180.0;
        gamma *= mathematicalConstant::pi/180.0;
    }

    tensor::operator=
    (
        tensor
        (
            cos(alpha)*cos(beta),
            cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma),
            cos(alpha)*sin(beta)*cos(gamma) - sin(alpha)*sin(gamma),

            sin(alpha)*cos(beta),
            sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),
            sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma),

            -sin(beta),
            cos(beta)*sin(gamma),
            cos(beta)*cos(gamma)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TaitBryanCoordinateRotation::TaitBryanCoordinateRotation()
:
    coordinateRotation()
{}


Foam::TaitBryanCoordinateRotation::TaitBryanCoordinateRotation
(
    const vector& phiThetaPsi,
    const bool inDegrees
)
:
    coordinateRotation()
{
    calcTransform
    (
        phiThetaPsi.component(vector::X),
        phiThetaPsi.component(vector::Y),
        phiThetaPsi.component(vector::Z),
        inDegrees
    );
}


Foam::TaitBryanCoordinateRotation::TaitBryanCoordinateRotation
(
    const scalar phiAngle,
    const scalar thetaAngle,
    const scalar psiAngle,
    const bool inDegrees
)
:
    coordinateRotation()
{
    calcTransform(phiAngle, thetaAngle, psiAngle, inDegrees);
}


Foam::TaitBryanCoordinateRotation::TaitBryanCoordinateRotation
(
    const dictionary& dict
)
:
    coordinateRotation()
{
    vector rotation(dict.lookup("rotation"));

    bool inDegrees = true;
    if (dict.found("degrees"))
    {
        inDegrees = Switch(dict.lookup("degrees"));
    }

    calcTransform
    (
        rotation.component(vector::X),
        rotation.component(vector::Y),
        rotation.component(vector::Z),
        inDegrees
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
