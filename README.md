# Developing and validating the BIKER algorithm

This repo is for developing and validating an algorithm to remotely sense the gas exchange velocity via satellite observations made by the NASA/CNES/UKSA/CSA SWOT satellite. The model is called BIKER, or the 'Bayesian Inference/Inversion of the $k_{600}$ Rate' algorithm. Note that the actual algorithm is written as an R package and available at <https://github.com/craigbrinkerhoff/BIKER>. It is called directly in this project's validation workflow.

## Notes

-   This workflow was built to run in parallel on an HPC cluster. If you run into problems on a local machine, don't say I didn't warn you! <br>
-   User specific inputs are at `/lib/inputs.R`. Here is where you can change the number of cores requested to run the parallel processes on.
-   The gas exchange model is developed and validated via the `~src/swot_k_model_vX.R`. This script is run independently of the algorithm validation (detailed next) and ran on a local machine.
-   The virtual environment used to run these analyses on a Linux HPC is reproducible using `environment.yml`. Note that, within that conda environment, R packages were generally (though not exclusively) installed into a private library via the `renv` package. See `config/global.dcf` for the full list of required packages and check that they are all installed either via `conda-forge` or `renv`.

## To run
- Submit a job to your HPC using the `run.sh` and `src/main.R` scripts.

## To generate manuscript

-   The manuscript for this project is written in RMarkdown. The most recent version is stored in the main directory, while old versions are stored in `/docs/`.
-   To generate the manuscript as a `.docx` file (pre-formatted correctly for submission to an AGU journal), open `RSK600.Rproj` just run `generateMS.R`. This will pull all figures and results from the analysis outputs in order to have 100% reproducibility.

## Questions

Feel free to reach out at cbrinkerhoff[at]umass[dot]edu!
