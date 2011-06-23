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
    floatingBodies

Description
    A list of floating bodies.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef floatingBodies_H
#define floatingBodies_H

#include "PtrList.H"
#include "floatingBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class floatingBodies Declaration
\*---------------------------------------------------------------------------*/

class floatingBodies
:
    public PtrList<floatingBody>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        floatingBodies(const floatingBodies&);

        //- Disallow default bitwise assignment
        void operator=(const floatingBodies&);


public:

    // Constructors

        //- Construct from Istream
        floatingBodies
        (
            const fvMesh& mesh,
            Istream& is
        )
        {
            PtrList<entry> entries(is);

            setSize(entries.size());

            forAll (entries, i)
            {
                this->set
                (
                    i,
                    new floatingBody
                    (
                        entries[i].keyword(),
                        mesh,
                        entries[i].dict()
                    )
                );
            }
        }


    // Destructor - default


    // Member Functions

        //- Set boundary motion for all bodies
        void setMotion(tetPointVectorField& motionU)
        {
            PtrList<floatingBody>& bodies = *this;

            forAll (bodies, bodyI)
            {
                bodies[bodyI].setMotion(motionU);
            }
        }

        //- Correct motion for all bodies
        void correctMotion(vectorField& newAllPoints) const
        {
            const PtrList<floatingBody>& bodies = *this;

            forAll (bodies, bodyI)
            {
                bodies[bodyI].correctMotion(newAllPoints);
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
