{
    surfaceScalarField phialpha1 =
        fvc::flux
        (
            phi,
            alpha1,
            alphaScheme
        )
      + fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

    MULES::explicitSolve(alpha1, phi, phialpha1, 1, 0);
}
