## Strategy

Use *applications/solvers/lagrangian/coalChemistryFoam* as a starting point. 
We need to collect the features we need to add/modify.

## Notes

Features to implement - we should find the correct order to implement them:
* Secondary breakup model: Junjun used the Stochastic Secondary Droplet (SSD) Model, Paolo wanted to use the one he published in JFM.
* Char combustion: surface reaction implemented in coalChemistryFoam already, this should suffice for the time being, once we know more on the char morphology we can work on diffusion-limited models or refine the existing surface model
* Pyrolysis model: Elia's kinetics - not really sure how to add this in a Lagrangian framework, requires more detailed model of the liquid phase I guess. We have to account for the four SARA fractions, pyrolysis products are char and volatiles that are vaporized.
* Gas-phase combustion: pyrolysis products vaporize from the droplets and then burn.
