/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::RASModels::DurbinV2f

Description
    Durbin V2f turbulence model for incompressible flows

    The default model coefficients correspondg to the following:
    \verbatim
        DurbinV2fCoeffs
	{
	    Cmu
	}
    \endverbatim

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
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class DurbinV2f Declaration
\*---------------------------------------------------------------------------*/

class DurbinV2f
:
    public RASModel
{

protected:

    // Protected data

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
            dimensionedScalar sigmaEps_;
            dimensionedScalar yStarLim_;

            dimensionedScalar fMin_;

            wallDist yw_;

        // Fields
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField v2_;
            volScalarField f_;

            volScalarField nut_;


    // Protected member functions

        tmp<volScalarField> T() const;
        tmp<volScalarField> L() const;

public:

    //- Runtime type information
    TypeName("DurbinV2f");

    // Constructors

        //- from components
        DurbinV2f
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    // Destructor

    virtual ~DurbinV2f()
    {}


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
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Note that epsilonTilda is returned as epsilon.
        //  This is the appropriate variable for most purposes.
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return turbulent velocity scale v^2
        virtual tmp<volScalarField> v2() const
        {
            return v2_;
        }

        //- Return 
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;
  
        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;
    
        // virtual tmp<fvVectorMatrix> divR(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
  //} // End namespace turbulenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
