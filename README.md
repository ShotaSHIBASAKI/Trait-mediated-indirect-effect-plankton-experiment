# Trait-mediated-indirect-effect-plankton-experiment
TThis repository contains the data and codes used in the manuscript "Quantifying a trait-mediated indirect effect of an environmental stressor on population dynamics" by Shota Shibasaki and Masato Yamamichi (2025, bioRxiv), https://doi.org/10.1101/2025.09.04.674171.

# File discription
## ./DataAnalysis
DataAnalysis. R replicates the analyses in the manuscript ecept for the parameter estimation.
The  first chunk "Clump size comparaison" shows the distribution of the clump sizes of the three algae species under four conditions, using the data in the folder "./ClumpSize", which continas the clump sizes in csv files.
There are four csv files unsed in this analysis, corresponding to the clump sizes of the three species without salt or rotifers (None), with rotifers (Predator), with salt (Salt), and with the two stressors (PredatorSlat).
This cunk also perform teh two-way aligned rank transformation ANOVA (i.e., a non-parametric two-way ANOVA) to see how the culture conditions affect the clump sizes of each green alga.

The second chunk compared the distributioned of the clump sizes of Chlamydomonas sphaeroides after 6 hour with or withour rotifers, using  "./ClumpSize/N2242_ClumpSize_Salt_Predator.csv".
The codes also peform LMM to see how the log2 of clump sizes changed due to the time and the presence of the rotifers.

The third chunk plot show the three-month dynamics of the mean clump diameters and show that the phenotypic plasticity causes the clump formation (see Supporting information), using the data in "ExperimentalDataN2242 - CellCulture.csv" and "ExperimentalDataN2242 - Plasticity.csv".

## ./Algae-rotifer-Rstan/ParameterEstimation.Rmd
This file performs the analysis of parameter estimation using the Bayesian statistics with Rstan. All data and stan files are contined in the "./Algae-rotifer-Rstan/".

