# Remote sensing the gas exchange velocity and carbon efflux using the SWOT Satellite

This project deals with developing an algorithm which uses Bayesian inference to remotely sense k600 in preparation for the launch of the NASA/CNES/UKSA/CSA SWOT mission. This project is a work in progress but is completely reproducible and will generate all results (and the manuscript) by following the directions below.

## To run
I suggest running most of the analysis on an HPC to take advantage of all the parallelized processes in this code (as well as big memory needs for BIKER). Some of these scripts will crash an underpowered local machine in their current setups. Further, the manuscript knitting must be performed on a local machine unless pandoc is set up on your HPC.

*On local machine*
1) Set number of cores within '~/scripts/inputs.R' to the number you want to run the parallel processes on.
2) Open RSK600.Rproj in RStudio
3) Run main.R inside of this project. This will run all analyses, generate all results, and knit the manuscript. But, it will likely crash/take forever on 'normal' local machines.

*On HPC cluster*
1) Set number of cores within '~/scripts/inputs.R' to the number you want to run the parallel processes on.
2) set working directory to '~/RSK600/' (make sure your R installation is activated!)
3) run 'run.sh' as a bash script. This will run all analyses and generate all results. This will run easily and time will be a function of your HPC session setup, but it will take much longer if generating new Monte Carlo results (by default munge2 == 0 and this is not performed as I have included these results in the repo). But, you must clone the repo locally to knit the manuscript.

## Manuscript
The manuscript for this project is under development and written in RMarkdown. The most recent version is stored in the project directory, while old versions are hosted in the 'manuscripts_vc' folder.

The manuscript file draws all figures and results from the analysis outputs in order to have 100% reproducibility. Just use main.R!! Note that it can only be generated on an HPC if pandoc is installed.

## Notes
1) The actual algorithm is written as an R package and available at https://github.com/craigbrinkerhoff/BIKER and is called directly in the validation script in this analysis.

3) Note that the hydraulics dataset needed to run the MC simulations is not included here as it is too large to store on Git. I have included the amended version that the MC simulations are ran on. Set munge == 1 if you desire to regenerate this amended version, but you will need to store it in the following location for all code to run smoothly: ~inputs/MonteCarlo/field_measurements.csv

## To do
1) rerun all analysis and make sure it's working smoothly
2) Write!
