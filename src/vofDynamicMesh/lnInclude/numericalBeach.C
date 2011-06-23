/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "numericalBeach.H"
#include "numericalBeachFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

numericalBeach::numericalBeach(volVectorField& U)
:
    U_(U)
{}


numericalBeach::~numericalBeach()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> numericalBeach::damping() const
{
    tmp<volScalarField> tDamping
    (
        new volScalarField
        (
            IOobject
            (
                "damping",
                U_.db().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0)
        )
    );

    volScalarField& damping = tDamping();

    const volVectorField::GeometricBoundaryField& bvf = U_.boundaryField();
    forAll(bvf, patchi)
    {
        if (isType<numericalBeachFvPatchField>(bvf[patchi]))
        {
            const numericalBeachFvPatchField& beach =
                dynamic_cast<const numericalBeachFvPatchField&> (bvf[patchi]);

            damping = max(damping, beach.internalDamping());
        }
    }

    return tDamping;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
