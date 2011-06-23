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

#include "graphSolidMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "axisCoordinateRotation.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(graphSolidMotion, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, graphSolidMotion, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::graphSolidMotion::centre() const
{
    const scalar t = time().value();

    return motionCentre_
      + vector
        (
            interpolateXY(t, x_.x(), x_.y()),
            interpolateXY(t, y_.x(), y_.y()),
            interpolateXY(t, z_.x(), z_.y())
        );
}


Foam::vector Foam::graphSolidMotion::rotation() const
{
    const scalar t = time().value();

    scalar xAngle = interpolateXY(t, xRot_.x(), xRot_.y());
    scalar yAngle = interpolateXY(t, yRot_.x(), yRot_.y());
    scalar zAngle = interpolateXY(t, zRot_.x(), zRot_.y());

    if (inDegrees_)
    {
        xAngle *= mathematicalConstant::pi/180.0;
        yAngle *= mathematicalConstant::pi/180.0;
        zAngle *= mathematicalConstant::pi/180.0;
    }

    return vector(xAngle, yAngle, zAngle);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::graphSolidMotion::graphSolidMotion(const IOobject& io)
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
    x_
    (
        "t",
        "x-translation",
        "x",
        IFstream(time().path()/time().caseConstant()/word(dict_.lookup("x")))()
    ),
    y_
    (
        "t",
        "y-translation",
        "y",
        IFstream(time().path()/time().caseConstant()/word(dict_.lookup("y")))()
    ),
    z_
    (
        "t",
        "z-translation",
        "z",
        IFstream(time().path()/time().caseConstant()/word(dict_.lookup("z")))()
    ),
    xRot_
    (
        "t",
        "x-rotation",
        "x",
        IFstream(time().path()/time().caseConstant()/word(dict_.lookup("xRot")))()
    ),
    yRot_
    (
        "t",
        "y-rotation",
        "yRot",
        IFstream(time().path()/time().caseConstant()/word(dict_.lookup("yRot")))()
    ),
    zRot_
    (
        "t",
        "z-rotation",
        "zRot",
        IFstream(time().path()/time().caseConstant()/word(dict_.lookup("zRot")))()
    ),
    inDegrees_(dict_.lookup("inDegrees"))
{
    refPoints_.rename("refPoints");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::graphSolidMotion::~graphSolidMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::graphSolidMotion::update()
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
