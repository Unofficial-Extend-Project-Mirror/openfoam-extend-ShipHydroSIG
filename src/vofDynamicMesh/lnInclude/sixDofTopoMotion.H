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
    sixDofTopoMotion

Description
    Multi-body motion with sliding interfaces.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    sixDofTopoMotion.C

\*---------------------------------------------------------------------------*/

#ifndef sixDofTopoMotion_H
#define sixDofTopoMotion_H

#include "topoChangerFvMesh.H"
#include "floatingBodies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                        Class sixDofTopoMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDofTopoMotion
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Floating bodies
        floatingBodies bodies_;

        //- Use topological sliding.  Other option is GGI
        const Switch useTopoSliding_;

        //- Mesh motion solver
        autoPtr<motionSolver> motionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDofTopoMotion(const sixDofTopoMotion&);

        //- Disallow default bitwise assignment
        void operator=(const sixDofTopoMotion&);


        // Topology controls

            //- Add sliding interface zones and modifiers
            void addZonesAndModifiers();

            //- Return true if sliders are attached
            bool attached() const;


public:

    //- Runtime type information
    TypeName("sixDofTopoMotion");


    // Constructors

        //- Construct from IOobject
        explicit sixDofTopoMotion(const IOobject& io);


    // Destructor

        ~sixDofTopoMotion();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
