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
    graphSolidMotion

Description
    Automatic motion of the mesh given graphs for
    6-DOF translation and 3-axis rotation
 
SourceFiles
    graphSolidMotion.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef graphSolidMotion_H
#define graphSolidMotion_H

#include "dynamicFvMesh.H"
#include "dimensionedTypes.H"
#include "Switch.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class graphSolidMotion Declaration
\*---------------------------------------------------------------------------*/

class graphSolidMotion
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

            //- x-translation
            graph x_;

            //- y-translation
            graph y_;

            //- x-translation
            graph z_;


        // Rotation
    
            //- x-rotation
            graph xRot_;

            //- y-rotation
            graph yRot_;

            //- z-rotation
            graph zRot_;

            //- Is rotation given in degrees?
            Switch inDegrees_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        graphSolidMotion(const graphSolidMotion&);

        //- Disallow default bitwise assignment
        void operator=(const graphSolidMotion&);


        // Return current position of centroid
        vector centre() const;

        //- Return current local-to-global transform
        vector rotation() const;


public:

    //- Runtime type information
    TypeName("graphSolidMotion");


    // Constructors

        //- Construct from IOobject
        graphSolidMotion(const IOobject& io);


    // Destructor

        ~graphSolidMotion();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
