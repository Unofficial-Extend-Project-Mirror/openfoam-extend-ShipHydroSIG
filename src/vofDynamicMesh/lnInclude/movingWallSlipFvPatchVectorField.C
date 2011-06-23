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

#include "movingWallSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingWallSlipFvPatchVectorField::movingWallSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    basicSymmetryFvPatchVectorField(p, iF)
{}


movingWallSlipFvPatchVectorField::movingWallSlipFvPatchVectorField
(
    const movingWallSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basicSymmetryFvPatchVectorField(ptf, p, iF, mapper)
{}


movingWallSlipFvPatchVectorField::movingWallSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


movingWallSlipFvPatchVectorField::movingWallSlipFvPatchVectorField
(
    const movingWallSlipFvPatchVectorField& pivpvf
)
:
    basicSymmetryFvPatchVectorField(pivpvf)
{}


movingWallSlipFvPatchVectorField::movingWallSlipFvPatchVectorField
(
    const movingWallSlipFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    basicSymmetryFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingWallSlipFvPatchVectorField::evaluate()
{
    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();
    const fvMesh& mesh = dimensionedInternalField().mesh();

    const pointField& oldAllPoints = mesh.oldAllPoints();

    vectorField oldFc(pp.size());

    forAll(oldFc, i)
    {
        oldFc[i] = pp[i].centre(oldAllPoints);
    }

    vectorField Up = (pp.faceCentres() - oldFc)/mesh.time().deltaT().value();

    const volVectorField& U =
        db().lookupObject<volVectorField>(dimensionedInternalField().name());
    scalarField phip = 
        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

    vectorField nHat = p.nf();
    const scalarField& magSf = p.magSf();
    scalarField Un = phip/(magSf + VSMALL);

    vectorField tangentialPart =
    (
        this->patchInternalField()
      + transform(I - 2.0*sqr(nHat), this->patchInternalField())
    )/2.0;

    vectorField normalPart = Up + nHat*(Un - (nHat & Up));

    vectorField::operator=(normalPart + tangentialPart);
}


void movingWallSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingWallSlipFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
