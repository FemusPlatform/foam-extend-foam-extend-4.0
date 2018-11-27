/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "externalHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

externalHeatFluxTemperatureFvPatchScalarField::
externalHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    alphaEffName_("undefinedAlphaEff"),
    CpName_("undefinedCp"),
    h_ex_(p.size(), 1.),
    Tref_(p.size(), 0.0)
{}


externalHeatFluxTemperatureFvPatchScalarField::
externalHeatFluxTemperatureFvPatchScalarField
(
    const externalHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(ptf.q_, mapper),
    alphaEffName_(ptf.alphaEffName_),
    CpName_(ptf.CpName_),
    h_ex_(ptf.h_ex_),
    Tref_(ptf.Tref_, mapper)
{}


externalHeatFluxTemperatureFvPatchScalarField::
externalHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_("q", dict, p.size()),
    Tref_("Tref", dict, p.size()),
    h_ex_("h_ex", dict, p.size()),   
    alphaEffName_(dict.lookup("alphaEff")),
    CpName_(dict.lookup("Cp"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


externalHeatFluxTemperatureFvPatchScalarField::
externalHeatFluxTemperatureFvPatchScalarField
(
    const externalHeatFluxTemperatureFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_),
    CpName_(thftpsf.CpName_),
    h_ex_(thftpsf.h_ex_),
    Tref_(thftpsf.Tref_)
{}


externalHeatFluxTemperatureFvPatchScalarField::
externalHeatFluxTemperatureFvPatchScalarField
(
    const externalHeatFluxTemperatureFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_),
    CpName_(thftpsf.CpName_),
    h_ex_(thftpsf.h_ex_),
    Tref_(thftpsf.Tref_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void externalHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    q_.autoMap(m);
    Tref_.autoMap(m);
    h_ex_.autoMap(m);
}


void externalHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const externalHeatFluxTemperatureFvPatchScalarField& thftptf =
        refCast<const externalHeatFluxTemperatureFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
    Tref_.rmap(thftptf.Tref_, addr);
    h_ex_.rmap(thftptf.h_ex_, addr);
}


void externalHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tbound =
        lookupPatchField<volScalarField, scalar>("T");
    
    const scalarField& alphaEffp =
        lookupPatchField<volScalarField, scalar>(alphaEffName_);

    const scalarField& Cpp =
        lookupPatchField<volScalarField, scalar>(CpName_);

    q_ = h_ex_*(Tbound - Tref_);
    
    gradient() = q_/(Cpp*alphaEffp);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void externalHeatFluxTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    q_.writeEntry("q", os);
    os.writeKeyword("alphaEff") << alphaEffName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cp") << CpName_ << token::END_STATEMENT << nl;
    h_ex_.writeEntry("h_ex", os);
    Tref_.writeEntry("Tref", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    externalHeatFluxTemperatureFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam


// ************************************************************************* //

