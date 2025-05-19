# MetabolicNetworks
This repository is in supplementation to my MSc Thesis Project "Developing Enzyme-Constrained Genome-Scale Metabolic Model Pipeline to Decipher Microbiome Metabolic Alterations". Here you can find relevant scripts and results.

fva.py:
This python script automatically determines which fva results are still not created and calls fva.m for corresponding samples. It also performs a rudimentary differential abundance analysis (not included in the Thesis report), and creates several csv files containing p-values from DA, taxonomy csv, and 2 abundance dictionaries.

fva.m:
This MatLab function script is what is called by fva.py to perform fva. The script does all the necessary fva calculations and saves the results. Raven Toolbox is required for usage, Gurobi Optimization software is highly recommended. 

DA.R:
This R script performs differential abundance analysis for 2 cohorts with edgeR and DESeq2.
