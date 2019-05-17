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

#include "exchangeWithSolidWithFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOReferencer.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

exchangeWithSolidWithFluxFvPatchScalarField::
exchangeWithSolidWithFluxFvPatchScalarField
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
    solidMass_(p.size(), 0.0),
    extFlux_(p.size(), 0.0)
{

}


exchangeWithSolidWithFluxFvPatchScalarField::
exchangeWithSolidWithFluxFvPatchScalarField
(
    const exchangeWithSolidWithFluxFvPatchScalarField& ptf,
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
    solidMass_(ptf.solidMass_, mapper),
    extFlux_(ptf.solidMass_, mapper)
{

}


exchangeWithSolidWithFluxFvPatchScalarField::
exchangeWithSolidWithFluxFvPatchScalarField
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
    CpName_(dict.lookup("Cp")),
    extFlux_("extFlux", dict, p.size())
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;

}


exchangeWithSolidWithFluxFvPatchScalarField::
exchangeWithSolidWithFluxFvPatchScalarField
(
    const exchangeWithSolidWithFluxFvPatchScalarField& thftpsf
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
    solidMass_(thftpsf.solidMass_),
    extFlux_(thftpsf.extFlux_)
{

}


exchangeWithSolidWithFluxFvPatchScalarField::
exchangeWithSolidWithFluxFvPatchScalarField
(
    const exchangeWithSolidWithFluxFvPatchScalarField& thftpsf,
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
    solidMass_(thftpsf.solidMass_),
    extFlux_(thftpsf.extFlux_)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void exchangeWithSolidWithFluxFvPatchScalarField::autoMap
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
    extFlux_.autoMap(m);
}


void exchangeWithSolidWithFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const exchangeWithSolidWithFluxFvPatchScalarField& thftptf =
        refCast<const exchangeWithSolidWithFluxFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
    Tref_.rmap(thftptf.Tref_, addr);
    h_ex_.rmap(thftptf.h_ex_, addr);
    solidTemperature_.rmap(thftptf.solidTemperature_, addr);
    solidCp_.rmap(thftptf.solidCp_, addr);
    solidMass_.rmap(thftptf.solidMass_, addr);
    extFlux_.rmap(thftptf.extFlux_, addr);

}


void exchangeWithSolidWithFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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

    solidCp_ =   450. 
               + 0.280*solidTemperature_ 
               - 2.91*1.e-4*solidTemperature_*solidTemperature_
               + 1.34*1.e-7*solidTemperature_*solidTemperature_*solidTemperature_;
    
    solidTemperature_ -= (TotalFlux - extFlux_*Area) * dt / (solidMass_ * solidCp_);
        
    Info<<"Patch "<<patch().name();
    
    int counter = 0;
    forAll(q_, faceI)
    {
        if(counter == 0) {// this allows to print info when running in parallel
            Info<<" Tsolid "<< solidTemperature_[faceI] << " Sol-Fluid Flux "<<TotalFlux/Area<<" Sol Flux "<<extFlux_[faceI]<<"\n";
        }
        counter ++;
        label faceCellI = patch().faceCells()[faceI];
        q_[faceI] = -(Tvolume[faceCellI] - solidTemperature_[faceI]) * Cpp[faceI]*alphaEffp[faceI] * (this->patch().deltaCoeffs()[faceI]);
    }
        
    gradient() = q_/(Cpp*alphaEffp);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void exchangeWithSolidWithFluxFvPatchScalarField::write(Ostream& os) const
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
    extFlux_.writeEntry("extFlux", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    exchangeWithSolidWithFluxFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam


// ************************************************************************* //

