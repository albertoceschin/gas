## Description

* From tutorial: OpenFOAM-v2212/tutorials/lagrangian/coalChemistryFoam/simplifiedSiwek. Goal is to adapt it as much as possible to KAUST gasifier, in order to understand what is missing and has to be implemented.

## Tasks 

- [x] Change names of patches in /save and /constant in order to adapt tutorial to the pipe mesh. 
- [x] Adapt geometry in blockMeshDict to something meaningful. 
- [x] Changed injection to cone injection.
- [x] Test if reaction on, if conditions similar to the tutorial.
- [ ] Test multiphase reaction.
- [ ] Update species and flowrates (coal -> VRO). Idea was to use one liquid, one solid and multiple gas species to implement missing reactions. For now, realistic but no real.

## Notes

* In tutorial cell set for spark takes only one cell.
* Combustion allowed after once all volatile components evolved. In detail, *residualCoeff* constant in *devolatilisationModel*.
