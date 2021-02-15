# Remote sensing the gas exchange velocity and carbon efflux using the SWOT Satellite

This project deals with developing an empirical model for the normalized gas exchange velocity k600, and a subsequent algorithm which uses Bayesian inference to remotely sense k600 in preparation for the launch of the NASA/CNES/UKSA/CSA SWOT mission. This project is a work in progress but is completely reproducible and will generate all results (and the manuscript) by following the directions below.

## To run
I suggest running most of the analysis on an HPC to take advantage of parallelization and memory issues with BIKER, i.e both the runBIKER.R and FCO2_implications.R scripts are written to be ran on an HPC and likely will crash many local machines in their current setups. However, both the Monte Carlo analysis and manuscript knitting must be performed on a local machine.

*On local machine*
1) Set number of cores within '~/scripts/inputs.R' to the number you want to run the parallel processes on.
2) Open RSK600.Rproj in RStudio
3) Run main.R inside of this project.

*On HPC cluster*
1) Set number of cores within '~/scripts/inputs.R' to the number you want to run the parallel processes on.
2) set working directory to '~/RSK600/' (make sure your R installation is activated!)
3) run 'run.sh' as a bash script

## Manuscript
The manuscript for this project is under development and written in RMarkdown. The most recent version is stored in the project directory, while old versions are hosted in the 'manuscripts_vc' folder.

The manuscript file draws all figures and results from the analysis outputs in order to have 100% reproducibility. Just use main.R!! Note that it can't be generated on an HPC currently because of issues with pandoc.

## Notes
1) The actual algorithm is written as an R package and available at https://github.com/craigbrinkerhoff/BIKER and is called directly in the validation script in this analysis.

3) Note that the hydraulics dataset needed to run the MC simulations is not included here as it is too large to store on Git. It is available in dropbox from its paper's original repo (https://github.com/craigbrinkerhoff/2019_AMHG) but you will need to store it in the following location for all code to run smoothly: ~inputs/MonteCarlo/field_measurements.csv

## To do
1) Decide on slope errors
2) Write!
