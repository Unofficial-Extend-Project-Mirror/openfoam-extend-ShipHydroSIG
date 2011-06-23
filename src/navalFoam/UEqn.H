    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    // Convection term to be added separately, depending on singlePhase switch
    // HJ, 4/Aug/2008
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      - fvm::laplacian(muEff, U, "laplacian(mut,U)")
      - (fvc::grad(U) & fvc::grad(muEff))
      ==
      - alpha1*rho*fvm::Sp(beach.damping(), U)
    );


    Switch singlePhase(pimple.lookup("singlePhase"));

    if (singlePhase)
    {
        UEqn += pos(alpha1 - 0.1)*fvm::div(rhoPhi, U);
//         UEqn += alpha1*fvm::div(rhoPhi, U);
    }
    else
    {
        UEqn += fvm::div(rhoPhi, U);
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
        )
    );

    if (singlePhase)
    {
        U *= pos(alpha1 - 0.1);
        U.correctBoundaryConditions();
    }
