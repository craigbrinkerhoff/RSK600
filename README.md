# Developing and validating the BIKER algorithm

This repo is for developing and validating a model to remotely sense the gas exchange velocity via satellite observations made by the NASA/CNES/UKSA/CSA SWOT satellite. The model is called BIKER, or the 'Bayesian Inference/Inversion of the $k_{600}$ Rate' algorithm. This project is completely reproducible and will generate all results (and the manuscript) by following the directions below.

Note that the actual algorithm is written as an R package and available at <https://github.com/craigbrinkerhoff/BIKER>. It is called directly in this project's validation workflow.

## Notes

- This workflow was built to run in parallel on an HPC cluster. If you run into problems on a local machine, don't say I didn't warn you! <br>
- You need to have the R library `ProjectTemplate` already installed to run this code. Also check out the other required libraries at `/config/global.dcf` <br>
- User specific inputs are at `/lib/inputs.R`. Here is where you can change the number of cores requested to run the parallel processes on.

## To run validation

- The gas exchange model is developed and validated via the '\~src/swot_k\_model.R'. This script is run independently of the algorithm validation (detailed next) and is not set up for an HPC cluster. To run the algorithm validation, set working directory to `/RSK600/`. Then, execute 'run.sh' as a bash script. Options for how this job is submitted are stored in that file as well and are up to the user.

## To generate manuscript

-   The manuscript for this project is written in RMarkdown. The most recent version is stored in the main directory, while old versions are stored in `/docs/`.
-   To generate the manuscript as a `.docx` file (pre-formatted correctly for submission to an AGU journal), open `RSK600.Rproj` just run `generateMS.R`. This will pull all figures and results from the analysis outputs in order to have 100% reproducibility.

## Questions

Feel free to reach out at cbrinkerhoff[at]umass[dot]edu!
