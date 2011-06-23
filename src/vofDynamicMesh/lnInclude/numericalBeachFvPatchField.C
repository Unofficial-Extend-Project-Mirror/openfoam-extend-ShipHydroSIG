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

#include "numericalBeachFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "patchWave.H"
#include "mathematicalConstants.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::mathematicalConstant;

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

numericalBeachFvPatchField::numericalBeachFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchField<vector>(p, iF),
    thickness_(0.0),
    offset_(0.0),
    nud_(0.0)
{}


numericalBeachFvPatchField::numericalBeachFvPatchField
(
    const numericalBeachFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchField<vector>(ptf, p, iF, mapper),
    thickness_(ptf.thickness_),
    offset_(ptf.offset_),
    nud_(ptf.nud_)
{}


numericalBeachFvPatchField::numericalBeachFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchField<vector>(p, iF, dict),
    thickness_(readScalar(dict.lookup("thickness"))),
    nud_(readScalar(dict.lookup("nud")))
{
    if (dict.found("offset"))
    {
        offset_ = readScalar(dict.lookup("offset", 0.0));

        if (offset_ > thickness_)
        {
            FatalErrorIn
            (
                "numericalBeachFvPatchField::numericalBeachFvPatchField("
                "const fvPatch& p,"
                "const DimensionedField<vector, volMesh>& iF,"
                "const dictionary& dict)"
            )
            << "Offset must be smaller than thickness"
            << exit(FatalError);
        }
    }
}


numericalBeachFvPatchField::numericalBeachFvPatchField
(
    const numericalBeachFvPatchField& pivpvf
)
:
    inletOutletFvPatchField<vector>(pivpvf),
    thickness_(pivpvf.thickness_),
    offset_(pivpvf.offset_),
    nud_(pivpvf.nud_)
{}


numericalBeachFvPatchField::numericalBeachFvPatchField
(
    const numericalBeachFvPatchField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchField<vector>(pivpvf, iF),
    thickness_(pivpvf.thickness_),
    offset_(pivpvf.offset_),
    nud_(pivpvf.nud_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void numericalBeachFvPatchField::updateCoeffs()
{
    inletOutletFvPatchField<vector>::updateCoeffs();
}


tmp<volScalarField> numericalBeachFvPatchField::internalDamping() const
{
    // Get patchids of walls
    labelHashSet wallPatchIDs(1);
    wallPatchIDs.insert(this->patch().index());

    // Calculate distance starting from wallPatch faces.
    patchWave wave(patch().boundaryMesh().mesh(), wallPatchIDs, false);

    const scalarField& distance = wave.distance();

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    tmp<volScalarField> tDamping
    (
        new volScalarField
        (
            IOobject
            (
                "damping",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0)
        )
    );

    volScalarField& damping = tDamping();

    forAll(distance, cellI)
    {
        if ( distance[cellI] < offset_ )
        {
            damping[cellI] = nud_;
        }
        else if ( distance[cellI] < thickness_ )
        {
            scalar nonDimDist =
                (distance[cellI] - offset_)/(thickness_ - offset_);

            damping[cellI] = nud_*(0.5 + 0.5*cos(pi*nonDimDist));
        }
    }

    return tDamping;
}


void
numericalBeachFvPatchField::write(Ostream& os) const
{
    inletOutletFvPatchField<vector>::write(os);

    os.writeKeyword("thickness") << thickness_ << token::END_STATEMENT << nl;
    if ( offset_ > 0 )
    {
        os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("nud") << nud_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    numericalBeachFvPatchField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
