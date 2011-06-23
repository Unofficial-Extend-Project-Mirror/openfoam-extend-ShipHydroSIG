    surfaceScalarField nuEff
    (
        "nuEff",
        twoPhaseProperties.nuf()
      + fvc::interpolate(turbulence->nut())
    );

    // Convection term to be added separately, depending on singlePhase switch
    // HJ, 4/Aug/2008
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      - fvm::laplacian(nuEff, U)
      - (fvc::grad(U) & fvc::grad(nuEff))
      ==
      - alpha1*fvm::Sp(beach.damping(), U)
    );


    Switch singlePhase(pimple.lookup("singlePhase"));

    if (singlePhase)
    {
        UEqn += pos(alpha1 - 0.1)*fvm::div(phi, U);
//         UEqn += alpha1*fvm::div(phi, U);
    }
    else
    {
        UEqn += fvm::div(phi, U);
    }

    if (oCorr == nOuterCorr - 1)
    {
        UEqn.relax(1);
    }
    else
    {
        UEqn.relax();
    }

    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(pd)
            )*mesh.magSf()
        )/rho
    );

    if (singlePhase)
    {
        U *= pos(alpha1 - 0.1);
        U.correctBoundaryConditions();
    }
