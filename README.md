# Remote sensing the gas exchange velocity and carbon efflux using the SWOT Satellite

This project deals with devolping an empirical model for the normalized gas exchange velocity $k_{600}$, and a subsequent algorithm which uses Bayesian inference to remotely sense $k_{600}$ in preperation for the launch of the NASA/CNES/UKSA/CSA SWOT mission. This project is a work in progress.

## To run
*On local machine*
1) Open RSK600.Rproj in RStudio
2) Run main.R inside of this project.
3) That's it!

*on HPC cluster*
1) set working directory to '~/RSK600/'
2) run 'run.sh'
3) That's it!

Note that this bash script isn't running the MC analysis because it requires an external dataset too big for git (and an R package I can't get installed)...

# Manuscript
The manuscript for this project is under development and written in RMarkdown. The most recent version is stored in the project directory, while old versions are hosted in the manuscripts folder.

The manuscript file draws all figures and results from the analysis outputs in order to have 100% reproducibility. Just use main.R!! Note that it can't be generated on an HPC currently.

## Notes
1) The actual algorithm is written as an R package and available at https://github.com/craigbrinkerhoff/BIKER. The script needed to run these tests is included here (runBIKER.R), though it is not run in the main.R script and is set up specifically for implementation on a high performance computing cluster.

2) The stan model implemented in the algorithm is saved here purely for posterity's sake (along with past versions as this is all a work in progress). The most recent version is the one that is implemented in the BIKER R package.

3) Note that the hydraulics dataset needed to run the MC simulations is not included here as it is too large to store on Git. It is availble in dropbox from its paper's original repo (https://github.com/craigbrinkerhoff/2019_AMHG) but you will need to store it in the following location for all code to run smoothly: ~inputs/MonteCarlo/field_measurements.csv

## To do
1) Decide on slope errors
2) Write!