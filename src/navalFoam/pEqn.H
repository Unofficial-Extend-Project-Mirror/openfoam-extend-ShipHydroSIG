if (implicitP)
{
    // Add gravity term, re-using the momentum equation
    UEqn += interface.sigmaK()*fvc::grad(alpha1) - rho*g;

    volScalarField rUAPres = 1.0/UEqn.A();

    // Pick up velocity boundary conditions
    volVectorField Upres
    (
        "Upres",
        U
    );
    Upres = rUAPres*UEqn.H();

    surfaceScalarField phiPres = fvc::interpolate(Upres) & mesh.Sf();

    adjustPhi(phiPres, Upres, p);

    p.storePrevIter();

    for(int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAPres, p) == fvc::div(phiPres)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve(mesh.solver(pd.name()));
    }

    Info << "Relaxation: p = " << mesh.relaxationFactor(p.name()) << endl;
    p.relax();
}
else
{
    p = pd + rho*gh;
}
