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

Class
    Foam::exchangeWithSolidWithFluxFvPatchScalarField

Description
    Fixed heat flux boundary condition for temperature.

SourceFiles
    exchangeWithSolidWithFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef exchangeWithSolidWithFluxFvPatchScalarFields_H
#define exchangeWithSolidWithFluxFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
       Class exchangeWithSolidWithFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class exchangeWithSolidWithFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
// Private data

    //- Heat flux [W/m2]
    scalarField q_;
    scalarField Tref_;
    scalarField h_ex_;
    scalarField solidCp_;
    scalarField solidMass_;
    scalarField solidTemperature_;
    scalarField extFlux_;
    
    //- Name of effective thermal diffusivity field
    word alphaEffName_;
    
    //- Name of specific heat capacity field
    word CpName_;
public:
    
    //- Runtime type information
    TypeName("exchangeWithSolidWithFlux");


    // Constructors

        //- Construct from patch and internal field
        exchangeWithSolidWithFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        exchangeWithSolidWithFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  exchangeWithSolidWithFluxFvPatchScalarField onto
        //  a new patch
        exchangeWithSolidWithFluxFvPatchScalarField
        (
            const exchangeWithSolidWithFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        exchangeWithSolidWithFluxFvPatchScalarField
        (
            const exchangeWithSolidWithFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new exchangeWithSolidWithFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        exchangeWithSolidWithFluxFvPatchScalarField
        (
            const exchangeWithSolidWithFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new exchangeWithSolidWithFluxFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
