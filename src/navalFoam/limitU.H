{
    scalar limitMagU = readScalar(pimple.lookup("limitMagU"));

    volScalarField magU(mag(U));

    scalar maxMagU = max(magU).value();

    Info<< "mag(U): max: " << maxMagU
        << " avg: " << magU.weightedAverage(mesh.V()).value();

    if (maxMagU > limitMagU)
    {
        U.internalField() *=
            1.0 - pos(magU.internalField() - limitMagU)*
            limitMagU/magU.internalField();
        U.correctBoundaryConditions();
        Info << " ...limiting" << endl;
    }
    else
    {
        Info << endl;
    }
}
