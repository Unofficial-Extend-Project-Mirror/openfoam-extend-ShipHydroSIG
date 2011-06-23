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
    threeDofMotionSubset

SourceFiles
    threeDofMotionSubset.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef threeDofMotionSubset_H
#define threeDofMotionSubset_H

#include "dynamicFvMesh.H"
#include "translationODE.H"
#include "polyPatchID.H"
#include "ODESolver.H"
#include "fvMeshSubset.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class threeDofMotionSubset Declaration
\*---------------------------------------------------------------------------*/

class threeDofMotionSubset
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


        //- Subset mesh
        fvMeshSubset subsetMesh_;

        //- Mesh motion solver
        autoPtr<motionSolver> motionPtr_;

        //- Identifier of a hull patch
        polyPatchID hullPatch_;

        //- Identifier of a motion patch on a subset
        word motionPatch_;

        //- Centre of gravity
        vector centroid_;

        //- Current time index
        label curTimeIndex_;

        //- Force at the beginning of time-step
        vector oldForce_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        threeDofMotionSubset(const threeDofMotionSubset&);

        //- Disallow default bitwise assignment
        void operator=(const threeDofMotionSubset&);


        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return force on the hull
        vector hullForce() const;


public:

    //- Runtime type information
    TypeName("threeDofMotionSubset");


    // Constructors

        //- Construct from IOobject
        threeDofMotionSubset(const IOobject& io);


    // Destructor

        virtual ~threeDofMotionSubset();


    // Member Functions

        //- Return reference to subset mesh
        const fvMeshSubset& subsetMesh() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
