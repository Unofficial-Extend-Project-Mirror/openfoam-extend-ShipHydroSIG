{
#   include "continuityErrs.H"

    wordList pcorrTypes(pd.boundaryField().types());

    for (label i=0; i < pd.boundaryField().size(); i++)
    {
        if (pd.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", pd.dimensions(), 0.0),
        pcorrTypes
    );

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve(mesh.solver(pd.name()));

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pcorrEqn.flux();
        }

        // Fluxes are corrected to absolute velocity and further corrected
        // later.  HJ, 6/Feb/2009
    }

#   include "continuityErrs.H"
}
