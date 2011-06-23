{
    if (nOuterCorr != 1)
    {
        pd.storePrevIter();
    }

    rAU = 1.0/UEqn.A();
    surfaceScalarField rAUf = fvc::interpolate(rAU/rho);

    U = rAU*UEqn.H();
    U.correctBoundaryConditions();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
        // Restart accuracy lost with ddtCorr.  Removed, HJ, 16/Nov/2010
//       + fvc::ddtPhiCorr(rAU, U, phi)
    );

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf();

    if (!mesh.moving())
    {
        adjustPhi(phi, U, pd);
    }

    for(int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAUf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pRefCell, pRefValue);

        // Solve pressure equation to a tight tolerance in the last corrector
        if
        (
            oCorr == nOuterCorr - 1
         && corr == nCorr - 1
         && nonOrth == nNonOrthCorr
        )
        {
            pdEqn.solve(mesh.solver("pdFinal"));
        }
        else
        {
            pdEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pdEqn.flux();
        }
    }

    // Explicitly relax pressure except for last corrector
    if (oCorr != nOuterCorr - 1)
    {
        pd.relax();
    }

    U += rAU/rho*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();
}
