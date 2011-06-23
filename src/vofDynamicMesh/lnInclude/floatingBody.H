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
    floatingBody

Description
    Floating body is a representation of the mesh and 6-DOF forces with the
    following characteristics:

    1) Mesh representation is given as a list of patches representing the hull.
    Patches can be used to capture individual components (eg. keel,
    bulb, rudder, hull etc.)

    2) Mesh is structured in the following way
    Hull is created with symmetry around the y-z plane
    Incoming velocity is in -x direction
    Mesh can be of one of following types
    2.1) Single zone mesh, with or without prismatic layers
    2.2) Two-component GGI mesh with a spherical interface

    2) Flow force calculation is based on turbulent single-phase or
    free surface flow for all patches.  Only the force on wetted surfaces
    (gamma = 1) is considered

    3) External propulsion force is given as direction and height above
     centre of gravity.  Force magnitude equals to total drag force
    in opposing direction

    4) 6-DOF solver
    Holds inertial parameters of the boat
    Force balance includes hydrodynamic forces and external propulsion force
    Motion can be constrained in each of the following:
    Forward motion (x)
    Sideways motion (y)
    Vertical motion (z)
    Roll angle (x-rot)
    Pitch angle (y-rot)
    Yaw angle (z-rot)

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    floatingBody.C

\*---------------------------------------------------------------------------*/

#ifndef floatingBody_H
#define floatingBody_H

#include "word.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "fvMesh.H"
#include "sixDOFqODE.H"
#include "ODESolver.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "uniformDimensionedFields.H"

#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class floatingBody Declaration
\*---------------------------------------------------------------------------*/

class floatingBody
{
    // Private data

        //- Name of floating body
        word name_;

        //- Reference to mesh
        const fvMesh& mesh_;


        // Mesh variables

            //- Hull patches
            PtrList<polyPatchID> hullPatches_;

            //- Hull slider patch: moving with hull translation and rotation
            polyPatchID hullSlider_;

            //- Hull slider zone: moving with hull translation and rotation
            faceZoneID hullSliderZone_;

            //- Fixed slider patch: moving with hull translation only
            polyPatchID fixedSlider_;

            //- Fixed slider zone: moving with hull translation only
            faceZoneID fixedSliderZone_;

            //- Centre of slider rotation relative to centre of gravity
            vector centreOfSlider_;


        // External force variables

            //- Add propulsion force
            Switch addPropulsionForce_;

            //- Propulsion direction
            vector propulsionDirection_;

            //- Propulsion force point relative to centre of gravity
            vector propulsionForceCentroid_;


        // ODE parameters

            //- 6-DOF motion equation
            sixDOFqODE sixDOF_;

            //- ODE solver
            autoPtr<ODESolver> solver_;

            //- ODE solver accuracy
            scalar epsilon_;

            //- Fixed motion components

                //- Fixed surge (x-translation)
                Switch fixedSurge_;

                //- Fixed sway (y-translation)
                Switch fixedSway_;

                //- Fixed heave (z-translation)
                Switch fixedHeave_;

                //- Fixed roll (rotation around x)
                Switch fixedRoll_;

                //- Fixed pitch (rotation around y)
                Switch fixedPitch_;

                //- Fixed yaw (rotation around z)
                Switch fixedYaw_;


        //- Steady-state 6-DOF simulation
        Switch steadyState_;

        // Helper variables

            //- Current time index
            label curTimeIndex_;

            //- Current centre of rotation
            vector rotCentre_;

            //- Force at the beginning of time-step
            vector oldForce_;

            //- Moment at the beginning of time-step
            vector oldMoment_;


        //- Output file for translation and rotation
        OFstream of_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        floatingBody(const floatingBody&);

        //- Disallow default bitwise assignment
        void operator=(const floatingBody&);


        //- Check definition
        void check() const;

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return drag force on the hull
        vector dragForce() const;

        //- Return drag moment on the hull
        vector dragMoment() const;

        //- Fix translation
        void fixTranslation(vector& v) const;

        //- Fix rotation
        void fixRotation(vector& rot) const;

        //- Return velocity of the centre of rotation
        vector rotCentreVelocity() const;

        //- Calculate boundary velocity
        tmp<vectorField> boundaryVelocity
        (
            const vectorField& boundaryPoints
        ) const;


public:


    // Constructors

        //- Construct from dictionary
        floatingBody
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor - default


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return ID of hull slider patch
        const polyPatchID& hullSlider() const
        {
            return hullSlider_;
        }

        //- Return ID of fixed slider patch
        const polyPatchID& fixedSlider() const
        {
            return fixedSlider_;
        }

        //- Set boundary motion
        void setMotion(tetPointVectorField& motionU);

        //- Correct motion.  This will move zones in parallel
        void correctMotion(vectorField& newAllPoints) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
