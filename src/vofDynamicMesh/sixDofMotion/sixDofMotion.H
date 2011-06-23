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
    sixDofMotion

Description

SourceFiles
    sixDofMotion.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef sixDofMotion_H
#define sixDofMotion_H

#include "dynamicFvMesh.H"
#include "floatingBodies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                        Class sixDofMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDofMotion
:
    public dynamicFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Floating bodies
        floatingBodies bodies_;

        //- Mesh motion solver
        autoPtr<motionSolver> motionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDofMotion(const sixDofMotion&);

        //- Disallow default bitwise assignment
        void operator=(const sixDofMotion&);


public:

    //- Runtime type information
    TypeName("sixDofMotion");


    // Constructors

        //- Construct from IOobject
        sixDofMotion(const IOobject& io);


    // Destructor

        ~sixDofMotion();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
