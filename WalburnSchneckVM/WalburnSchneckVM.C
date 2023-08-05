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

\*---------------------------------------------------------------------------*/

#include "WalburnSchneckVM.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(WalburnSchneckVM, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        WalburnSchneckVM,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::WalburnSchneckVM::calcNu() const
{
    dimensionedScalar dimMu("dimMu", dimMass/dimLength/dimTime, 1.0);
    dimensionedScalar dimSR("dimMu", dimless/dimTime, 1.0);
    dimensionedScalar dimLambda("dimLambda", dimless/dimTime, 1);

    return (dimMu/rho_) * a1_ * exp(a2_*hct_*100) * exp(a4_*(TPMA_/pow(hct_*100,2))) * pow((dimLambda+strainRate())/dimSR, -1*a3_*hct_*100);
}

//*(1 + 0.0*strainRate()/dimLambda);//
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::WalburnSchneckVM::WalburnSchneckVM
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    WalburnSchneckVMCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    a1_(WalburnSchneckVMCoeffs_.lookup("a1")),
    a2_(WalburnSchneckVMCoeffs_.lookup("a2")),
    a3_(WalburnSchneckVMCoeffs_.lookup("a3")),
    a4_(WalburnSchneckVMCoeffs_.lookup("a4")),
    TPMA_(WalburnSchneckVMCoeffs_.lookup("TPMA")),
    rho_(WalburnSchneckVMCoeffs_.lookup("rho")),
    hct_(WalburnSchneckVMCoeffs_.lookup("hct")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::WalburnSchneckVM::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    WalburnSchneckVMCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    WalburnSchneckVMCoeffs_.lookup("a1") >> a1_;
    WalburnSchneckVMCoeffs_.lookup("a2") >> a2_;
    WalburnSchneckVMCoeffs_.lookup("a3") >> a3_;
    WalburnSchneckVMCoeffs_.lookup("a4") >> a4_;
    WalburnSchneckVMCoeffs_.lookup("TPMA") >> TPMA_;
    WalburnSchneckVMCoeffs_.lookup("rho") >> rho_;
    WalburnSchneckVMCoeffs_.lookup("hct") >> hct_;

    return true;
}


// ************************************************************************* //
