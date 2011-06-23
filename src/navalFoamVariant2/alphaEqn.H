{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    label nAlphaSubCycles
    (
        readLabel(pimple.lookup("nAlphaSubCycles"))
    );

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            if (implicitAlpha)
            {
#               include "alphaImplicit.H"
            }
            else
            {
#               include "alphaMules.H"
            }
        }
    }
    else
    {
        if (implicitAlpha)
        {
#           include "alphaImplicit.H"
        }
        else
        {
#           include "alphaMules.H"
        }
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;

    interface.correct();

    // Calculate density using limited alpha1
    rho = twoPhaseProperties.rho();
}
