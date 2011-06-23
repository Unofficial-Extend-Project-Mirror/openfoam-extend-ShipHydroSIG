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

#include "linearWaveFvPatchField.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    valueAbove_(pTraits<Type>::zero),
    valueBelow_(pTraits<Type>::zero),
    amplitude_(0),
    period_(1),
    axis_("x"),
    curTimeIndex_(-1),
    phiName_("none")
{}


template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    valueAbove_(pTraits<Type>(dict.lookup("valueAbove"))),
    valueBelow_(pTraits<Type>(dict.lookup("valueBelow"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    period_(readScalar(dict.lookup("period"))),
    axis_(dict.lookup("axis")),
    curTimeIndex_(-1),
    phiName_(dict.lookupOrDefault<word>("phiName","none"))
{
    this->refValue() = (valueBelow_);

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(valueBelow_);
    }

    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1.0;
}


template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const linearWaveFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    valueAbove_(ptf.valueAbove_),
    valueBelow_(ptf.valueBelow_),
    amplitude_(ptf.amplitude_),
    period_(ptf.period_),
    axis_(ptf.axis_),
    curTimeIndex_(-1),
    phiName_(ptf.phiName_)
{}


template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const linearWaveFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    valueAbove_(ptf.valueAbove_),
    valueBelow_(ptf.valueBelow_),
    amplitude_(ptf.amplitude_),
    period_(ptf.period_),
    axis_(ptf.axis_),
    curTimeIndex_(-1),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
template<class Type>
void linearWaveFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        //Field<Type>& patchField = *this;

        scalar twoPi = 2*mathematicalConstant::pi;

        scalarField coord;

        if (axis_ == "x")
        {
            coord = this->patch().Cf().component(vector::X);
        }
        else if (axis_ == "y")
        {
            coord = this->patch().Cf().component(vector::Y);
        }
        else if (axis_ == "z")
        {
            coord = this->patch().Cf().component(vector::Z);
        }
        else
        {
            FatalErrorIn("void linearWaveFvPatchField<Type>::updateCoeffs()")
                << "Unknown axis: " << axis_ << ".  Should be x, y or z"
                << abort(FatalError);
        }

        scalar h = amplitude_*
            Foam::sin(twoPi*this->db().time().value()/period_);

        Info<< "t = " << this->db().time().value()
            << " h = " << h << endl;

        this->refValue() =
            pos(coord - h)*valueAbove_ + neg(coord - h)*valueBelow_;

        curTimeIndex_ = this->db().time().timeIndex();
    }

    if (phiName_ != "none")
    {
        const Field<scalar>& phip = this->patch().lookupPatchField
        (
            phiName_,
            reinterpret_cast<const surfaceScalarField*>(NULL),
            reinterpret_cast<const scalar*>(NULL)
        );

        this->valueFraction() = 1.0 - pos(phip);
    }
    else
    {
        this->valueFraction() = 1.0;
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


// Write
template<class Type>
void linearWaveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("valueAbove")
        << valueAbove_ << token::END_STATEMENT << nl;
    os.writeKeyword("valueBelow")
        << valueBelow_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("period")
        << period_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis")
        << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("phiName")
        << phiName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
