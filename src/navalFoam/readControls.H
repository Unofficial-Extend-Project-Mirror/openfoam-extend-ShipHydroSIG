#   include "readTimeControls.H"
#   include "readPIMPLEControls.H"

    bool implicitAlpha = true;
    if (pimple.found("implicitAlpha"))
    {
        implicitAlpha = Switch(pimple.lookup("implicitAlpha"));
    }

    bool implicitP = true;
    if (pimple.found("implicitP"))
    {
        implicitP = Switch(pimple.lookup("implicitP"));
    }

    bool correctPhi = false;
    if (pimple.found("correctPhi"))
    {
        correctPhi = Switch(pimple.lookup("correctPhi"));
    }

    bool checkMeshCourantNo = false;
    if (pimple.found("checkMeshCourantNo"))
    {
        checkMeshCourantNo = Switch(pimple.lookup("checkMeshCourantNo"));
    }

    scalar startMotionTime(readScalar(pimple.lookup("startMotionTime")));
