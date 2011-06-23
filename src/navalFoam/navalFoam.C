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

Application
    interDyMFoam

Description
    Solver for 2 incompressible fluids capturing the interface with
    dynamic mesh motion.

    Turbulence modelling is generic, i.e.  laminar, RAS or LES may be selected.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "dynamicFvMesh.H"
#include "interfaceProperties.H"
#include "subCycle.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "numericalBeach.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"
#   include "readControls.H"
#   include "correctPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Execute function objects
    runTime.functionObjects().execute();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (runTime.value() > startMotionTime)
        {
            bool meshChanged = mesh.update();

#           include "volContinuity.H"

            // Mesh motion update
            if (correctPhi && meshChanged)
            {
#               include "correctPhi.H"
            }

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

            if (meshChanged)
            {
#               include "CourantNo.H"
            }
        }

        // Pressure-velocity corrector
        int oCorr = 0;
        do
        {
            gh = (g & mesh.C());
            ghf = (g & mesh.Cf());

            twoPhaseProperties.correct();

#           include "UEqn.H"

            // --- PISO loop
            for (int corr=0; corr < nCorr; corr++)
            {
#               include "pdEqn.H"
            }

#           include "limitU.H"

            if (mesh.moving())
            {
                // Make the fluxes relative
                phi -= fvc::meshPhi(rho, U);
                rhoPhi = phi*fvc::interpolate(rho);
            }

#           include "movingMeshRhoUContinuityErrs.H"

            // Calculate static pressure
            if (oCorr == nOuterCorr - 1)
            {
                // Calculate static pressure
                Info << "Calculating static pressure" << endl;
#               include "pEqn.H"
            }

#           include "alphaEqn.H"

            turbulence->correct();
        } while (++oCorr < nOuterCorr);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
