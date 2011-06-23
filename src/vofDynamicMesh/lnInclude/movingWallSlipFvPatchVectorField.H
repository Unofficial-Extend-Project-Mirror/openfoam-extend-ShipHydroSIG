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
    movingWallSlipFvPatchVectorField

Description
    Moving slip wall

SourceFiles
    movingWallSlipFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingWallSlipFvPatchVectorField_H
#define movingWallSlipFvPatchVectorField_H

#include "fvPatchFields.H"
#include "basicSymmetryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class movingWallSlipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class movingWallSlipFvPatchVectorField
:
    public basicSymmetryFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("movingWallSlip");


    // Constructors

        //- Construct from patch and internal field
        movingWallSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingWallSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingWallSlipFvPatchVectorField
        //  onto a new patch
        movingWallSlipFvPatchVectorField
        (
            const movingWallSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingWallSlipFvPatchVectorField
        (
            const movingWallSlipFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingWallSlipFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingWallSlipFvPatchVectorField
        (
            const movingWallSlipFvPatchVectorField&,
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
                new movingWallSlipFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
