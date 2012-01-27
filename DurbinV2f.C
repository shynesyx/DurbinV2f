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

\*---------------------------------------------------------------------------*/

#include "DurbinV2f.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DurbinV2f, 0);
addToRunTimeSelectionTable(RASModel, DurbinV2f, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> DurbinV2f::T() const
{
    volScalarField yStar_=pow(CmuKE_,0.25)*sqrt(k_)*yw_/nu();
    return max(k_/(epsilon_ + epsilonMin_),
               pos(yStarLim_-yStar_)*6.0*sqrt(nu()/(epsilon_ + epsilonMin_)));
}

tmp<volScalarField> DurbinV2f::L() const
{
    volScalarField yStar_=pow(CmuKE_,0.25)*sqrt(k_)*yw_/nu();
    return
        CL_*max(pow(k_,1.5)/(epsilon_ + epsilonMin_),
               pos(yStarLim_-yStar_)
               *Ceta_*pow(pow(nu(),3.0)/(epsilon_ + epsilonMin_),0.25));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
DurbinV2f::DurbinV2f
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    // Coefficients are obtained from Bad Nauheim, Deutschland 2005
    Cmu_    // C_mu of v2f
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),

    CmuKE_    // C_mu of k-epsilon
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKE",
            coeffDict_,
            0.09
        )
    ),
    Ceps10_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps10",
            coeffDict_,
            0.09
        )
     ),
    Ceps11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps11",
            coeffDict_,
            0.05
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.90
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.40
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.30
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.23
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceta",
            coeffDict_,
            70.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    yStarLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "yStarLim",
            coeffDict_,
            30.0
        )
    ),

    fMin_("fMinSmall", dimless/dimTime, SMALL),

    yw_(mesh_),

    // Field scalars
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    // Calculate viscosity (with Davidson correction - 2003)
    nut_(min(CmuKE_*sqr(k_)/(epsilon_ + epsilonMin_),
             Cmu_*v2_*T()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> DurbinV2f::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

tmp<volSymmTensorField> DurbinV2f::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

tmp<fvVectorMatrix> DurbinV2f::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

bool DurbinV2f::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        CmuKE_.readIfPresent(coeffDict());
	Ceps10_.readIfPresent(coeffDict());
        Ceps11_.readIfPresent(coeffDict());
	Ceps2_.readIfPresent(coeffDict());
	C1_.readIfPresent(coeffDict());
	C2_.readIfPresent(coeffDict());
	CL_.readIfPresent(coeffDict());
	Ceta_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
	yStarLim_.readIfPresent(coeffDict());
	
        return true;
    }
    else
    {
        return false;
    }
}

void DurbinV2f::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    // Production of TKE
    volScalarField G("RASModel::G", nut_*2*magSqr(symm(fvc::grad(U_))));

    // Coefficient
    const volScalarField Ceps1 = Ceps10_*(scalar(1.0)+Ceps11_*min(sqrt(k_/v2_),scalar(10.0)));

    // Dissipation rate equation

//#   include "epsilonV2fWallI.H"
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Ceps1*G/T()
      - fvm::Sp(Ceps2_/T(), epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

//#   include "wallDissipationV2fI.H"

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(1.0/T(), k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

    // f equation
    
    tmp<fvScalarMatrix> fEqn
    (
//      - fvm::laplacian(scalar(1.0),f_)
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L()),f_)
      - ((C1_-scalar(6.0))*v2_/k_ - 2.0/3.0*(C1_-scalar(1.0)))/(sqr(L())*T())
      + C2_*G/(k_*sqr(L()))
/*      - fvm::Sp(1.0/sqr(L()),f_)
      - ((C1_-scalar(1.0))*(v2_/k_-scalar(2.0/3.0))/T()
       - C2_*G/k_
       - 5.0*epsilon_*v2_/sqr(k_))/sqr(L())*/
    );

    fEqn().relax();
    solve(fEqn);
    bound(f_, fMin_);

    // v2 equation

    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      - fvm::laplacian(DkEff(), v2_)
     ==
    //    k_*f_ 
    // Davidson correction - 2003
        min(k_*f_, 
          -((C1_-scalar(6.0))*v2_ - 2.0/3.0*k_*(C1_-scalar(1.0)))/T() + C2_*G)
      - fvm::Sp(6.0*epsilon_/k_, v2_)
    );

    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, kMin_);

    // Re-calculate viscosity (with Davidson correction - 2003)
    nut_ = min(CmuKE_*sqr(k_)/(epsilon_ + epsilonMin_),
               Cmu_*v2_*T());

    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
