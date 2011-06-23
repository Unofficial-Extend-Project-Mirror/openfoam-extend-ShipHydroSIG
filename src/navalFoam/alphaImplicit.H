{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1, alphaScheme)
      + fvm::div
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        )
    );

    alpha1Eqn.solve();
}
