# PARSEC_ABC_FAR

## Codes for "Rationalised experiment design for parameter estimation with sensitivity clustering".

Problem statement (copied from the article): "We explore the dynamics of a three-gene repressilator system, a synthetic genetic circuit characterized by the sequential repression of gene expression in a cyclic manner (Elowitz, M.B. and Leibler, S., Nature, 403(6767), 2000). Using an ordinary differential equation (ODE) model detailed in Supplementary SI S3, we monitor the protein levels produced by these genes. We analyze a set of six simultaneous measurements of protein B and C levels over 72 hours to obtain an accurate estimate of the parameters using PARSEC. Briefly, PARSEC calculates the parameter sensitivity indices (PSI) at a specific set of time points for proteins B and C. We combine the indices for proteins B and C to create a PARSEC-PSI vector for the candidate measurements. We then use k-means clustering to group the measurements into six clusters, with one representative candidate randomly chosen from each cluster as part of the final design."

We compare the designs using ABC-FAR to identify the most informative ones.
Labels for the designs: SSD, SRD and W1SD denote the PARSEC, random and anti-PARSEC designs respectively.

## functions_ABC_FAR
This folder contains the codes for parameter estimation

## PARSEC_eg_Represilator
This folder has a solved/implemented example.
## PARSEC_template
This folder has the template.

## Contents of the last two folders:
### Files/folder in the template
1. #### designing_experiments.mlx - A live script to be executed. (Requires changes according to the problem statement)
2. model_system.m - The file where the model is to be written. (Requires changes based on model equations)
3. TSA - Folder containing codes to perform Sensitivity analysis.
4. PARSEC - Folder containing codes to perform clustering and design selection.
5. design_evaluation - Folder containing codes to evaluate designs.

### Files and folders created during execution (Default names)
1. GT_val.mat - Parameter values used for training and testing
2. TSA_full_*.mat - Data files storing results of sensititivity analysis.
3. data_set_* - Folder containing the excel sheet(s) describing the designs used in the analysis.
4. 
