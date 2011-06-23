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

Class
    Foam::numericalBeachFvPatchField

Description
    Foam::numericalBeachFvPatchField

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved.

SourceFiles
    numericalBeachFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef numericalBeachFvPatchField_H
#define numericalBeachFvPatchField_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class numericalBeachFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class numericalBeachFvPatchField
:
    public inletOutletFvPatchField<vector>
{
    // Private data

        //- Thickness measured from the patch
        scalar thickness_;

        //- Offset to allow some constant damping towards the walls
        scalar offset_;

        //- Damping factor
        scalar nud_;


public:

    //- Runtime type information
    TypeName("numericalBeach");


    // Constructors

        //- Construct from patch and internal field
        numericalBeachFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        numericalBeachFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  numericalBeachFvPatchField
        //  onto a new patch
        numericalBeachFvPatchField
        (
            const numericalBeachFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        numericalBeachFvPatchField
        (
            const numericalBeachFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new numericalBeachFvPatchField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        numericalBeachFvPatchField
        (
            const numericalBeachFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new numericalBeachFvPatchField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

        //- Return the damping in the internal field
        tmp<volScalarField> internalDamping() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
