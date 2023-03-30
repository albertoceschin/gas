## Strategy

* 27 March: copy COxidationKineticDiffusionLimitedRate and extend to Pyrolysis.
*src/lagrangian/intermediate/clouds/Templates/ReactingMultiphaseCloud/ReactingMultiphaseCloud.C* has most of the capabilities needed. 
It has also *src/lagrangian/intermediate/submodels/HeterogeneousReactingModel/MUCSheterogeneousRate* from D. Papanastassiou and G. Bitsianes, Modelling of Heterogeneous Gas-Solid Reactions.

Derived *src/lagrangian/coalCombustion/submodels/coalCloud* has following */surfaceReactionModel*:
* COxidationDiffusionLimitedRate
* COxidationHurtMitchell
* COxidationIntrinsicRate
* COxidationKineticDiffusionLimitedRate
* COxidationMurphyShaddix

Hence coalCloud can be copied and missing feature (ie. liquid to solid reactions) added.

## Notes

In the tutorial simplifiedSiwek, *limestoneParcels* is a *basicThermoCloud* defined in *applications/solvers/lagrangian/coalChemistryFoam/createClouds.H*.
On the opposite, *coalParcels* is a *coalCloud* with (ie) its singleMixtureFractionCoeffs for mass fractions and compositions of phases.

It seems a patch is needed. Explore if also gas, not only liquid, need fix. And if the fix is right.
