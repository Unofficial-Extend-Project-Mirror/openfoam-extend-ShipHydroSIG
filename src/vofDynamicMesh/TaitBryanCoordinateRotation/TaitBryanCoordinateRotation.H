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

Class
    Foam::TaitBryanCoordinateRotation

Description
    A coordinateRotation defined in the z-y'-x''  convention.

    The 3 rotations are defined in the TaitBryan convention
    (around Z, around Y' and around X') (Yaw, Pitch, Roll).
    Note that it is the reverse transformation
    (local->global) that is defined here.

    - the rotation angles are in degrees, unless otherwise explictly specified:

    @verbatim
        coordinateRotation
        {
            type        TaitBryanRotation
            degrees     false;
            rotation    (0 0 3.141592654);
        }
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef TaitBryanCoordinateRotation_H
#define TaitBryanCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class TaitBryanCoordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class TaitBryanCoordinateRotation
:
    public coordinateRotation
{
    // Private member functions

        //- Calculate transformation tensor
        void calcTransform
        (
            const scalar phiAngle,
            const scalar thetaAngle,
            const scalar psiAngle,
            const bool inDegrees = true
        );


public:

    //- Runtime type information
    TypeName("TaitBryanRotation");


    // Constructors

        //- Construct null
        TaitBryanCoordinateRotation();

        //- Construct from rotation vector
        TaitBryanCoordinateRotation
        (
            const vector& phiThetaPsi,
            const bool inDegrees = true
        );

        //- Construct from components of rotation vector
        TaitBryanCoordinateRotation
        (
            const scalar phiAngle,
            const scalar thetaAngle,
            const scalar psiAngle,
            const bool inDegrees = true
        );

        //- Construct from dictionary
        TaitBryanCoordinateRotation(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
