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
    floatingBodyFvMesh

Description
    Multi-body motion with GGI interfaces.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    floatingBodyFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef floatingBodyFvMesh_H
#define floatingBodyFvMesh_H

#include "topoChangerFvMesh.H"
#include "floatingBodies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                        Class floatingBodyFvMesh Declaration
\*---------------------------------------------------------------------------*/

class floatingBodyFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- List of bodies
        floatingBodies bodies_;

        //- Number of steps between mesh motion updates
        label nSteps_;

        //- Mesh motion solver
        autoPtr<motionSolver> motionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        floatingBodyFvMesh(const floatingBodyFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const floatingBodyFvMesh&);


public:

    //- Runtime type information
    TypeName("floatingBodyFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit floatingBodyFvMesh(const IOobject& io);


    // Destructor

        ~floatingBodyFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
