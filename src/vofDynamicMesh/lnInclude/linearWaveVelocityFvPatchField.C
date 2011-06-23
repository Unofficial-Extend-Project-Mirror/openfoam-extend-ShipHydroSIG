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

#include "linearWaveVelocityFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearWaveVelocityFvPatchField::linearWaveVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    amplitude_(0),
    period_(0),
    valueAbove_(vector::zero),
    curTimeIndex_(-1)
{}


linearWaveVelocityFvPatchField::linearWaveVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    period_(readScalar(dict.lookup("period"))),
    valueAbove_(dict.lookup("valueAbove")),
    curTimeIndex_(-1)
{
    evaluate();
}


linearWaveVelocityFvPatchField::linearWaveVelocityFvPatchField
(
    const linearWaveVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    period_(ptf.period_),
    valueAbove_(ptf.valueAbove_),
    curTimeIndex_(-1)
{}


linearWaveVelocityFvPatchField::linearWaveVelocityFvPatchField
(
    const linearWaveVelocityFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    amplitude_(ptf.amplitude_),
    period_(ptf.period_),
    valueAbove_(ptf.valueAbove_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void linearWaveVelocityFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        scalarField coord = this->patch().Cf().component(vector::Y);

        scalar omega = twoPi/period_;
        scalar lambda = 9.81*period_*period_/twoPi;
        scalar k = twoPi/lambda;

        scalar time = this->db().time().value();

        scalar h = amplitude_*sin(twoPi*time/period_);

        Info << "t = " << this->db().time().value()
            << " h = " << h
            << " U(0) = " <<  omega*cos(omega*time) << endl;

        forAll (coord, faceI)
        {
            vector newU = valueAbove_;

            scalar eky = exp(k*coord[faceI]);

            if(coord[faceI] > h)
            {
                eky = exp((h - coord[faceI])*3.0*k);
                newU.y() = 0.0;
            }
            else
            {
                newU.y() = omega*amplitude_*eky*cos(omega*time);
            }

            newU.x() = omega*amplitude_*eky*sin(omega*time);
            newU.z() = 0.0;

            this->operator[](faceI) = newU;
        };

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void linearWaveVelocityFvPatchField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("period")
        << period_ << token::END_STATEMENT << nl;
    os.writeKeyword("valueAbove")
        << valueAbove_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, linearWaveVelocityFvPatchField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
