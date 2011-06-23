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
    harmonicSolidMotion

Description
    Harmonic mesh motion: 6-DOF translation and rotation
    with independent amplitudes, periods and phase lags
 
SourceFiles
    harmonicSolidMotion.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef harmonicSolidMotion_H
#define harmonicSolidMotion_H

#include "dynamicFvMesh.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class harmonicSolidMotion Declaration
\*---------------------------------------------------------------------------*/

class harmonicSolidMotion
:
    public dynamicFvMesh
{
    // Private data

        //- Moving mesh dictionary
        dictionary dict_;

        //- Reference points from constant directory
        pointIOField refPoints_;

        //- Centre of motion
        vector motionCentre_;


        // Translation

            //- Translation amplitude
            vector transL_;

            //- Translation period
            vector transOmega_;

            //- Translation phase
            vector transPsi_;


        // Rotation

            //- Rotation amplitude
            vector rotL_;

            //- Rotation period
            vector rotOmega_;

            //- Rotation phase
            vector rotPsi_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        harmonicSolidMotion(const harmonicSolidMotion&);

        //- Disallow default bitwise assignment
        void operator=(const harmonicSolidMotion&);


        // Return current position of centroid
        vector centre() const;

        //- Return current local-to-global transform
        vector rotation() const;

        //- Convert period to frequency
        static vector periodToFreq(const vector deg);


public:

    //- Runtime type information
    TypeName("harmonicSolidMotion");


    // Constructors

        //- Construct from IOobject
        harmonicSolidMotion(const IOobject& io);


    // Destructor

        ~harmonicSolidMotion();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
