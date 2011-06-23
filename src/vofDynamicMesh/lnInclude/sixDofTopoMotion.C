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

\*---------------------------------------------------------------------------*/

#include "sixDofTopoMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "cellSet.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "fixedValueTetPolyPatchFields.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"

#include "slidingInterface.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDofTopoMotion, 0);

    addToRunTimeSelectionTable(dynamicFvMesh, sixDofTopoMotion, IOobject);
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDofTopoMotion::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if (useTopoSliding_)
    {
        if
        (
            pointZones().size() > 0
         || faceZones().size() > 0
         || cellZones().size() > 0
        )
        {
            Info<< "void sixDofTopoMotion::addZonesAndModifiers() : "
                << "Zones and modifiers already present.  Skipping."
                << endl;

            if (topoChanger_.size() == 0)
            {
                FatalErrorIn
                (
                    "void sixDofTopoMotion::addZonesAndModifiers()"
                )   << "Mesh modifiers not read properly"
                    << abort(FatalError);
            }

            return;
        }

        Info<< "Time = " << time().timeName() << endl
            << "Adding zones and modifiers to the mesh" << endl;

        // Add zones
        List<pointZone*> pz(3*bodies_.size());
        List<faceZone*> fz(3*bodies_.size());
        List<cellZone*> cz(0);

        label npz = 0;
        label nfz = 0;
        label nSliders = 0;

        forAll (bodies_, bodyI)
        {
            const floatingBody& curBody = bodies_[bodyI];

            if
            (
                curBody.hullSlider().active()
             && curBody.fixedSlider().active()
            )
            {
                nSliders++;

                // Add an empty zone for cut points
                pz[npz] = new pointZone
                (
                    curBody.name() + "CutPointZone",
                    labelList(0),
                    npz,
                    pointZones()
                );
                npz++;

                // Do face zones for slider

                // Inner slider
                const polyPatch& innerSlider =
                    boundaryMesh()[curBody.hullSlider().index()];

                labelList isf(innerSlider.size());

                forAll (isf, i)
                {
                    isf[i] = innerSlider.start() + i;
                }

                fz[nfz] = new faceZone
                (
                    curBody.name() + "InsideSliderZone",
                    isf,
                    boolList(innerSlider.size(), false),
                    nfz,
                    faceZones()
                );
                nfz++;

                // Outer slider
                const polyPatch& outerSlider =
                    boundaryMesh()[curBody.fixedSlider().index()];

                labelList osf(outerSlider.size());

                forAll (osf, i)
                {
                    osf[i] = outerSlider.start() + i;
                }

                fz[nfz] = new faceZone
                (
                    curBody.name() + "OutsideSliderZone",
                    osf,
                    boolList(outerSlider.size(), false),
                    nfz,
                    faceZones()
                );
                nfz++;

                // Add empty zone for cut faces
                fz[nfz] = new faceZone
                (
                    curBody.name() + "CutFaceZone",
                    labelList(0),
                    boolList(0, false),
                    nfz,
                    faceZones()
                );
                nfz++;
            }
        }

        pz.setSize(npz);
        fz.setSize(nfz);


        Info << "Adding point and face zones" << endl;
        addZones(pz, fz, cz);

        // Add topology modifiers
        topoChanger_.setSize(nSliders);
        label nTopos = 0;

        forAll (bodies_, bodyI)
        {
            const floatingBody& curBody = bodies_[bodyI];

            if
            (
                curBody.hullSlider().active()
             && curBody.fixedSlider().active()
            )
            {
                topoChanger_.set
                (
                    nTopos,
                    new slidingInterface
                    (
                        curBody.name() + "Slider",
                        nTopos,
                        topoChanger_,
                        curBody.name() + "OutsideSliderZone",
                        curBody.name() + "InsideSliderZone",
                        curBody.name() + "CutPointZone",
                        curBody.name() + "CutFaceZone",
                        curBody.fixedSlider().name(),
                        curBody.hullSlider().name(),
                        slidingInterface::INTEGRAL,   // Edge matching algorithm
                        true,                         // Attach-detach action
                        intersection::VISIBLE         // Projection algorithm
                    )
                );
            }
        }

        topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
        topoChanger_.write();
        write();
    }
    else if (topoChanger_.size() > 0)
    {
        if (topoChanger_.size() != 0)
        {
            FatalErrorIn
            (
                "void sixDofTopoMotion::addZonesAndModifiers()"
            )   << "Mesh modifiers are present and topo sliding is set to off"
                << abort(FatalError);
        }
    }
}


bool Foam::sixDofTopoMotion::attached() const
{
    const polyTopoChanger& topoChanges = topoChanger_;

    bool result = false;

    forAll (topoChanges, modI)
    {
        if (typeid(topoChanges[modI]) == typeid(slidingInterface))
        {
            result =
                result
             || refCast<const slidingInterface>(topoChanges[modI]).attached();
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDofTopoMotion::sixDofTopoMotion(const IOobject& io)
:
    topoChangerFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    bodies_(*this, dict_.lookup("bodies")),
    useTopoSliding_(dict_.lookup("useTopoSliding")),
    motionPtr_(motionSolver::New(*this))
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDofTopoMotion::~sixDofTopoMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sixDofTopoMotion::update()
{
    // Detaching the interface
    if (useTopoSliding_)
    {
        if (attached())
        {
            // Changing topology by hand
            autoPtr<mapPolyMesh> topoChangeMap1 = topoChanger_.changeMesh();

            if (topoChangeMap1.valid())
            {
                Info << "Decoupling sliding interfaces" << endl;
                motionPtr_->updateMesh(topoChangeMap1());
            }
        }
        else
        {
            Info << "Sliding interfaces decoupled" << endl;
        }
    }

    // Set the motion onto the motion patch
    tetPointVectorField& motionU =
        const_cast<tetPointVectorField&>
        (
            this->objectRegistry::lookupObject<tetPointVectorField>("motionU")
        );

    // Set motion from bodies
    bodies_.setMotion(motionU);

    // Save old points
    pointField oldPointsNew = allPoints();

    // Calculate the motion and move points
    movePoints(motionPtr_->newPoints());

    // Attach the interface
    // Changing topology by hand
    if (useTopoSliding_)
    {
        autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

        bool meshChanged2 = topoChangeMap2.valid();

        if (meshChanged2)
        {
            Info << "Coupling sliding interfaces" << endl;
            motionPtr_->updateMesh(topoChangeMap2());

            if (debug)
            {
                Info << "Moving points post slider attach" << endl;
            }

            pointField mappedOldPointsNew(allPoints().size());
            mappedOldPointsNew.map(oldPointsNew, topoChangeMap2->pointMap());

            movePoints(mappedOldPointsNew);
            resetMotion();
            setV0();

            movePoints(topoChangeMap2->preMotionPoints());
        }

        // Topological change: return true
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
