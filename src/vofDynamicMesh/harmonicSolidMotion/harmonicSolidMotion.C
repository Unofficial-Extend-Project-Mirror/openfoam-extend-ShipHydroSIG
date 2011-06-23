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

#include "harmonicSolidMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "axisCoordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(harmonicSolidMotion, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, harmonicSolidMotion, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::harmonicSolidMotion::periodToFreq(const vector period)
{
    return
        vector
        (
            2.0*mathematicalConstant::pi/(period.x() + SMALL),
            2.0*mathematicalConstant::pi/(period.y() + SMALL),
            2.0*mathematicalConstant::pi/(period.z() + SMALL)
        );
}


Foam::vector Foam::harmonicSolidMotion::centre() const
{
    const scalar t = time().value();

    return
        motionCentre_
      + cmptMultiply
        (
            transL_,
            vector
            (
                sin(transOmega_.x()*t + transPsi_.x()),
                sin(transOmega_.y()*t + transPsi_.y()),
                sin(transOmega_.z()*t + transPsi_.z())
            )
        );
}


Foam::vector Foam::harmonicSolidMotion::rotation() const
{
    const scalar t = time().value();

    return
        cmptMultiply
        (
            rotL_, 
            vector
            (
                sin(rotOmega_.x()*t + rotPsi_.x()),
                sin(rotOmega_.y()*t + rotPsi_.y()),
                sin(rotOmega_.z()*t + rotPsi_.z())
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::harmonicSolidMotion::harmonicSolidMotion(const IOobject& io)
:
    dynamicFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    refPoints_
    (
        IOobject
        (
            "points",
            time().constant(),
            polyMesh::meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    motionCentre_(dict_.lookup("motionCentre")),
    transL_(dict_.lookup("translationAmplitude")),
    transOmega_(periodToFreq(dict_.lookup("translationPeriod"))),
    transPsi_(dict_.lookup("translationPhase")),
    rotL_(dict_.lookup("rotationAmplitude")),
    rotOmega_(periodToFreq(dict_.lookup("rotationPeriod"))),
    rotPsi_(dict_.lookup("rotationPhase"))
{
    refPoints_.rename("refPoints");

    Switch inDegrees(dict_.lookup("inDegrees"));

    if (inDegrees)
    {
        transPsi_ *= mathematicalConstant::pi/180.0;
        rotL_ *= mathematicalConstant::pi/180.0;
        rotPsi_ *= mathematicalConstant::pi/180.0;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::harmonicSolidMotion::~harmonicSolidMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::harmonicSolidMotion::update()
{
    Info << "Updating point position" << endl;

    vector c = centre();
    vector r = rotation();
    tensor R = axisCoordinateRotation(r.x(), r.y(), r.z(), false);

    fvMesh::movePoints((R & (refPoints_ - motionCentre_)) + c);      

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //
