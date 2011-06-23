/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

#include "excentreRotationFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(excentreRotationFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, excentreRotationFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::excentreRotationFvMesh::excentreRotationFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
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
    rpm_
    (
        readScalar(dynamicMeshCoeffs_.lookup("rpm"))
    ),
    excentre_
    (
        readScalar(dynamicMeshCoeffs_.lookup("excentre"))
    ),
    rotAxis_
    (
       readScalar(dynamicMeshCoeffs_.lookup("rotAxis"))
    ),
    time0_
    (
        readScalar(dynamicMeshCoeffs_.lookup("time0"))
    ),
    time1_
    (
        readScalar(dynamicMeshCoeffs_.lookup("time1"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::excentreRotationFvMesh::~excentreRotationFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::excentreRotationFvMesh::update()
{
    scalar omega = 0.0;
      
    scalar localTime = time().value();

    scalar deltaT1 = time1_ - time0_ ;
      
    if (localTime >= 0.0 && localTime < time0_)
    {
        omega = 0.0;
    }
    else if (localTime < time1_)
    {
        scalar tempTime = localTime - time0_;

        omega = rpm_ * 
            (
                3.0*sqr(tempTime)/deltaT1/deltaT1
              - 2.0*pow3(tempTime)/pow3(deltaT1)
            );
    }
    else if (localTime > time1_)
    {
        omega = rpm_;           
    }

    Info << "Eccentricity = " << excentre_ << endl;
    Info << "  Local Time = " << localTime << endl;
    Info << "         rpm = " << omega << nl << endl;

    scalar theta = (omega*360.0*time().value()/60.0)*
        mathematicalConstant::pi/180.0;

    scalar deltaTheta = (omega*360.0*time().deltaT().value()/60.0)*
        mathematicalConstant::pi/180.0;
      
    const pointField& motionPoints = allPoints();

    pointField newPoints(motionPoints.size(), vector::zero);

    if (rotAxis_ == 1)
    {
        newPoints.replace
        (
            vector::X,
            motionPoints.component(vector::X)
          + excentre_*(cos(theta) - cos(theta + deltaTheta)) 
        ); 
      
        newPoints.replace
        (
            vector::Y,
            motionPoints.component(vector::Y)
        );

        newPoints.replace
        (
            vector::Z,
            motionPoints.component(vector::Z)
          + excentre_*(-sin(theta) + sin(theta + deltaTheta)) 
        );
    }
    else if (rotAxis_ == 2)
    {
        newPoints.replace
        (
            vector::X,
            motionPoints.component(vector::X)
          + excentre_*(cos(theta) - cos(theta + deltaTheta))
        );

        newPoints.replace
        (
            vector::Y,
            motionPoints.component(vector::Y)
          + excentre_*(-sin(theta) + sin(theta + deltaTheta))
        );

        newPoints.replace
        (
            vector::Z,
            motionPoints.component(vector::Z) 
        );
    }
    else if (rotAxis_ == 3)
    {
        newPoints.replace
        (
            vector::X,
            motionPoints.component(vector::X) 
        );

        newPoints.replace
        (
            vector::Y,
            motionPoints.component(vector::Y)
          + excentre_*(cos(theta) - cos(theta + deltaTheta))
        );

        newPoints.replace
        (
            vector::Z,
            motionPoints.component(vector::Z)
          + excentre_*(-sin(theta) + sin(theta + deltaTheta))
        );
    }

    fvMesh::movePoints(newPoints);      

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //
