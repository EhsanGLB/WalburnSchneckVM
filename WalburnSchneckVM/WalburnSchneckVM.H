/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::viscosityModels::WalburnSchneckVM

Description
    An incompressible Cross-Power law non-Newtonian viscosity model.

SourceFiles
    WalburnSchneckVM.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef WalburnSchneckVM_H
#define WalburnSchneckVM_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class WalburnSchneckVM Declaration
\*---------------------------------------------------------------------------*/

class WalburnSchneckVM
:
    public viscosityModel
{
    // Private data

        dictionary WalburnSchneckVMCoeffs_;

        dimensionedScalar a1_;
        dimensionedScalar a2_;
        dimensionedScalar a3_;
        dimensionedScalar a4_;
        dimensionedScalar TPMA_;
        dimensionedScalar rho_;
        dimensionedScalar hct_;

        volScalarField nu_;

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("WalburnSchneckVM");


    // Constructors

        //- Construct from components
        WalburnSchneckVM
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~WalburnSchneckVM()
        {}


    // Member Functions

        //- Return the laminar viscosity
        virtual const volScalarField& nu() const
        {
            return nu_;
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
