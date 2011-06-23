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
    rasVofForceAndTorqueFunctionObject

Description
    FunctionObject calculating force and torque on walls in turbulent
    VOF free surface flows.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    rasVofForceAndTorqueFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef rasVofForceAndTorqueFunctionObject_H
#define rasVofForceAndTorqueFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class rasVofForceAndTorqueFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class rasVofForceAndTorqueFunctionObject
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- List of patches
        wordList patchNames_;

        //- Origin
        point origin_;

        //- Output file
        OFstream of_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        rasVofForceAndTorqueFunctionObject
        (
            const rasVofForceAndTorqueFunctionObject&
        );

        //- Disallow default bitwise assignment
        void operator=(const rasVofForceAndTorqueFunctionObject&);

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

public:

    //- Runtime type information
    TypeName("rasVofForceAndTorque");


    // Constructors

        //- Construct from components
        rasVofForceAndTorqueFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
