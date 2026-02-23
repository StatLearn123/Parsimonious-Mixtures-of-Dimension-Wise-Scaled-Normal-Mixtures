# Parsimonious-Mixtures-of-Dimension-Wise-Scaled-Normal-Mixtures

This repository contains the code for fitting parsimonious mixtures of Dimension-Wise Scaled Shifted Exponential Normal distributions. 
Below is a description of the functions (and their arguments) provided in Main.R. 
To run the models, both Main.R and Utils.R must be sourced in the current R session.

## DSNM_M.fit ##

### Description ###

Fits, by using EM-based algorithms, parsimonious mixtures of Dimension-Wise Scaled Shifted Exponential Normal Mixtures. Parallel computing is implemented and highly recommended for faster calculation.

### Usage ###

DSNM_M.fit (X, k, corr = "all", scale = "all", tailedness = "all", rel.tol = 0.001, iter.max = 1000, verbose = TRUE, nThreads = 1)

### Arguments ###

* X: An n x p matrix, where n is the number of observations and p is the number of variables.
* k: A vector (or a number) containing the groups to be tried. If k is a vector, it must be 1:k. 
* corr: A character indicating the parsimonious structure for the reference correlation matrix P. Possible values are: "Identity", "Free", "E.comp", or "all", where "E.comp" stands for equality across components. When "all" is used, all the parsimonious structures are considered.
* scale: A character indicating the parsimonious structure for the reference marginal standard deviations matrix T. Possible values are: "Free", "E.dim", "E.comp", "E.comp.dim", or "all", where "E.dim" stands for equality across dimensions, and "E.comp.dim" stands for equality across components and dimensions. When "all" is used, all the parsimonious structures are considered.
* tailedness: A character indicating the parsimonious structure for the tailedness parameter \theta. Possible values are: "Free", "E.dim", "E.comp", "E.comp.dim", "Asympt", or "all", where "Asympt" stands for asymptotic value. When "all" is used, all the parsimonious structures are considered.
* rel.tol: A numeric indicating the stopping criterion for the algorithm.
* iter.max: A numeric indicating the maximum number of iterations of the algorithm.
* verbose: A logical specifying whether to display fitting information.
* nThreads: A positive integer indicating the number of cores used for running in parallel.  

## extract.bestM ##

### Description ###

This function extracts the best-fitting model according to the Bayesian information criterion (BIC) or Integrated classification likelihood (ICL).

### Usage ###

extract.bestM (results, sp.th = 0, criterion = "BIC")

### Arguments ###

* results: The output of the DSNM_M.fit() function.
* sp.th: Numeric threshold defining the minimum allowable mixture weight for each cluster, used to avoid potential spurious solutions. When set to 0 (default), the constraint is disabled.
* criterion: A character indicating the information criterion to be used. Possible values are "BIC" or "ICL".

# Example files

The "Sim 1_Example File.R" provides an example for the first simulation study of the paper. The "Sim 2_File.R" and "Real Data_File.R" provide the code for the second simulation study and the real-data analysis of the paper.
