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

#include "exchangeWithSolidFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

exchangeWithSolidFvPatchScalarField::
exchangeWithSolidFvPatchScalarField
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
    Tref_(p.size(), 0.0),
    solidTemperature_(p.size(), 0.0),
    solidCp_(p.size(), 0.0),
    solidMass_(p.size(), 0.0)
{

}


exchangeWithSolidFvPatchScalarField::
exchangeWithSolidFvPatchScalarField
(
    const exchangeWithSolidFvPatchScalarField& ptf,
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
    Tref_(ptf.Tref_, mapper),
    solidTemperature_(ptf.solidTemperature_, mapper),
    solidCp_(ptf.solidCp_, mapper),
    solidMass_(ptf.solidMass_, mapper)
{

}


exchangeWithSolidFvPatchScalarField::
exchangeWithSolidFvPatchScalarField
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
    solidTemperature_("solidTemperature", dict, p.size()),
    solidCp_("solidCp", dict, p.size()),
    solidMass_("solidMass", dict, p.size()),
    alphaEffName_(dict.lookup("alphaEff")),
    CpName_(dict.lookup("Cp"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;

}


exchangeWithSolidFvPatchScalarField::
exchangeWithSolidFvPatchScalarField
(
    const exchangeWithSolidFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_),
    CpName_(thftpsf.CpName_),
    h_ex_(thftpsf.h_ex_),
    Tref_(thftpsf.Tref_),
    solidTemperature_(thftpsf.solidTemperature_),
    solidCp_(thftpsf.solidCp_),
    solidMass_(thftpsf.solidMass_)
{

}


exchangeWithSolidFvPatchScalarField::
exchangeWithSolidFvPatchScalarField
(
    const exchangeWithSolidFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_),
    CpName_(thftpsf.CpName_),
    h_ex_(thftpsf.h_ex_),
    Tref_(thftpsf.Tref_),
    solidTemperature_(thftpsf.solidTemperature_),
    solidCp_(thftpsf.solidCp_),
    solidMass_(thftpsf.solidMass_)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void exchangeWithSolidFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    q_.autoMap(m);
    Tref_.autoMap(m);
    h_ex_.autoMap(m);
    solidTemperature_.autoMap(m);
    solidCp_.autoMap(m);
    solidMass_.autoMap(m);

}


void exchangeWithSolidFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const exchangeWithSolidFvPatchScalarField& thftptf =
        refCast<const exchangeWithSolidFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
    Tref_.rmap(thftptf.Tref_, addr);
    h_ex_.rmap(thftptf.h_ex_, addr);
    solidTemperature_.rmap(thftptf.solidTemperature_, addr);
    solidCp_.rmap(thftptf.solidCp_, addr);
    solidMass_.rmap(thftptf.solidMass_, addr);

}


void exchangeWithSolidFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tbound    = lookupPatchField<volScalarField, scalar>("T");
    const scalarField& alphaEffp = lookupPatchField<volScalarField, scalar>(alphaEffName_);
    const scalarField& Cpp       = lookupPatchField<volScalarField, scalar>(CpName_);
    const volScalarField& Tvolume = db().lookupObject<volScalarField>("T");
        
    scalar dt = db().time().deltaTValue();

    scalar TotalFlux = 0.;
    scalar Area      = 0.;
    
    // calculate solid temperature and heat flux
    forAll(q_, faceI)
    {
       TotalFlux += q_[faceI] * patch().magSf()[faceI];
       Area += patch().magSf()[faceI];     
    }
    
    Info<<" Tsolid before "<< solidTemperature_[0];
    
    solidTemperature_ -= TotalFlux * dt / (solidMass_ * solidCp_);
    
    Info<<" Tsolid after "<< solidTemperature_[0] << " Total Flux "<<TotalFlux<<" Area "<<Area<<"\n";
    
    forAll(q_, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];
        q_[faceI] = -(Tvolume[faceCellI] - solidTemperature_[faceI]) * Cpp[faceI]*alphaEffp[faceI] * (this->patch().deltaCoeffs()[faceI]);
    }
        
    gradient() = q_/(Cpp*alphaEffp);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void exchangeWithSolidFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    q_.writeEntry("q", os);
    os.writeKeyword("alphaEff") << alphaEffName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cp") << CpName_ << token::END_STATEMENT << nl;
    h_ex_.writeEntry("h_ex", os);
    Tref_.writeEntry("Tref", os);
    solidTemperature_.writeEntry("solidTemperature", os);
    solidCp_.writeEntry("solidCp", os);
    solidMass_.writeEntry("solidMass", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    exchangeWithSolidFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam


// ************************************************************************* //

