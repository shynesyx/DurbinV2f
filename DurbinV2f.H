/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    DurbinV2f

Description
    Durbin V2f turbulence model.

SourceFiles
    DurbinV2f.C
Author
    Yanxiang Shi

\*---------------------------------------------------------------------------*/

#ifndef DurbinV2f_H
#define DurbinV2f_H

//#include "turbulenceModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  //namespace turbulenceModels
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class DurbinV2f Declaration
\*---------------------------------------------------------------------------*/

class DurbinV2f
:
    public turbulenceModel
{
    // Private data

        // Model coefficients
            dimensionedScalar Cmu_;
            dimensionedScalar CmuKE_;
            dimensionedScalar Ceps10_;
            dimensionedScalar Ceps11_;
            dimensionedScalar Ceps2_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar oneOnSigmaK_;
            dimensionedScalar oneOnSigmaEps_;
            dimensionedScalar yStarLim_;

            dimensionedScalar f0_;

            wallDist yw_;

        // Fields
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField v2_;
            volScalarField f_;

            volScalarField nut_;


    // Private member functions

        tmp<volScalarField> T_() const;
        tmp<volScalarField> L_() const;

public:

    //- Runtime type information
    TypeName("DurbinV2f");

    // Constructors

        //- from components
        DurbinV2f
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~DurbinV2f(){}


    // Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DkEff", oneOnSigmaK*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DepsilonEff", oneOnSigmaEps*nut_ + nu())
            );
        }

        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Note that epsilonTilda is returned as epsilon.
        //  This is the appropriate variable for most purposes.
        tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        tmp<volScalarField> v2() const
        {
            return v2_;
        }

        tmp<volScalarField> f() const
        {
            return f_;
        }

        tmp<volSymmTensorField> R() const;

        tmp<fvVectorMatrix> divR(volVectorField& U) const;

        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
  //} // End namespace turbulenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
