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

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "tetFem.H"
#include "tetFec.H"
#include "coordinateSystem.H"
#include "TaitBryanCoordinateRotation.H"
#include "dynamicFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    // Create displacement field
    tetPointVectorField& motionU =
        const_cast<tetPointVectorField&>
        (
            mesh.lookupObject<tetPointVectorField>("motionU")
        );

    coordinateSystem cs
    (
        "LocalCS",
        vector(0, 0, 0),
//         TaitBryanCoordinateRotation(-5, 3.5, 5)
        TaitBryanCoordinateRotation(-5, 3.5, 0)
    );

    pointField lp = motionU.mesh().boundary()[0].localPoints();

    pointField targetPoints = cs.globalPosition(lp) + vector(0, 0, -0.036939);

//     for (runTime++; !runTime.end(); runTime++)

    runTime++;
    {
        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "Setting motion boundary condition on patch "
            << motionU.mesh()().boundaryMesh()[0].name() << endl;
        motionU.boundaryField()[0] ==
            (targetPoints - motionU.mesh().boundary()[0].localPoints())/
//             (runTime.endTime().value() - runTime.value());
            runTime.deltaT().value();

        mesh.update();

        mesh.checkMesh(true);

        mesh.write();
        runTime.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
