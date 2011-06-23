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
    threeDofMotion

Description
    Dynamic mesh class accomodating solid body motion based on a
    3-DOF translational solver.

SourceFiles
    threeDofMotion.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef threeDofMotion_H
#define threeDofMotion_H

#include "dynamicFvMesh.H"
#include "translationODE.H"
#include "polyPatchID.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                        Class threeDofMotion Declaration
\*---------------------------------------------------------------------------*/

class threeDofMotion
:
    public dynamicFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- 3-DOF equation
        translationODE equation_;

        //- ODE solver
        autoPtr<ODESolver> solver_;

        //- ODE solver accuracy
        scalar epsilon_;

        //- Mesh motion solver
        autoPtr<motionSolver> motionPtr_;

        //- Identifier of a hull patch
        polyPatchID hullPatch_;

        //- Current time index
        label curTimeIndex_;

        //- Force at the beginning of time-step
        vector oldForce_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        threeDofMotion(const threeDofMotion&);

        //- Disallow default bitwise assignment
        void operator=(const threeDofMotion&);


        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return force on the hull
        vector hullForce() const;


public:

    //- Runtime type information
    TypeName("threeDofMotion");


    // Constructors

        //- Construct from IOobject
        threeDofMotion(const IOobject& io);


    // Destructor

        virtual ~threeDofMotion();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
