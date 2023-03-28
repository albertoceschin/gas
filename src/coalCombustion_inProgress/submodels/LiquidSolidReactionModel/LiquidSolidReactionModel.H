/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::LiquidSolidReactionModel

Group
    grpLagrangianIntermediateLiquidSolidReactionSubModels

Description
    Templated surface reaction model class

SourceFiles
    LiquidSolidReactionModel.C
    LiquidSolidReactionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef LiquidSolidReactionModel_H
#define LiquidSolidReactionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class LiquidSolidReactionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LiquidSolidReactionModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected Data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //-Runtime type information
    TypeName("surfaceReactionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        LiquidSolidReactionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        LiquidSolidReactionModel(CloudType& owner);

        //- Construct from dictionary
        LiquidSolidReactionModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

        //- Construct copy
        LiquidSolidReactionModel(const LiquidSolidReactionModel<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<LiquidSolidReactionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~LiquidSolidReactionModel() = default;


    //- Selector
    static autoPtr<LiquidSolidReactionModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& cloud
    );


    // Member Functions

        //- Update surface reactions
        //  Returns the heat of reaction
        virtual scalar calculate
        (
            const scalar dt,
            const scalar Re,
            const scalar nu,
            const label celli,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const = 0;

        //- Add to devolatilisation mass
        void addToLiquidSolidReactionMass(const scalar dMass);


        // I-O

            //- Write injection info
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeLiquidSolidReactionModel(CloudType)                                    \
                                                                               \
    typedef Foam::CloudType::reactingMultiphaseCloudType                       \
        reactingMultiphaseCloudType;                                           \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::LiquidSolidReactionModel<reactingMultiphaseCloudType>,               \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            LiquidSolidReactionModel<reactingMultiphaseCloudType>,                 \
            dictionary                                                         \
        );                                                                     \
    }


#define makeLiquidSolidReactionModelType(SS, CloudType)                            \
                                                                               \
    typedef Foam::CloudType::reactingMultiphaseCloudType                       \
        reactingMultiphaseCloudType;                                           \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<reactingMultiphaseCloudType>, 0);                            \
                                                                               \
    Foam::LiquidSolidReactionModel<reactingMultiphaseCloudType>::                  \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<reactingMultiphaseCloudType>>                                \
        add##SS##CloudType##reactingMultiphaseCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LiquidSolidReactionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
