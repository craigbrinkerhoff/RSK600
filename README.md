# Remote sensing the gas transfer velocity using SWOT observations

See the 'outline.html' RMarkdown file, which explains this project in detail. In short, this project deals with devolping an empirical model for the normalized gas transfer velocity k600, and a subsequent algorithm which uses Bayesian inference to remotely sense k600 in preperation for the launch of the NASA/CNES/UKSA/CSA SWOT mission. This project is a work in progress.

##To run
1) run main.R
That's it!


##Notes
1) Included in this repo are validation results after running the algorithm on simulated SWOT data. The actual algorithm is at https://github.com/craigbrinkerhoff/BIGEER and the script needed to run these tests is included here (though is not run in the main.R script and is set up specifally for implementation on a high performance computing cluster).

2) Note that the hydraulics dataset needed to run the MC simulations is not included here as it is too large to store on Git. It is availble in dropbox from its paper's original repo (https://github.com/craigbrinkerhoff/2019_AMHG) but you will need to store it in the following location for all code to run smoothly: ~inputs/MonteCarlo/field_measurements.csv

##Things to do
1) Rerun entire validation with correct posterior uncertainity --> 0.25
2) Figure out a 3rd research goal.
