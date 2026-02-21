# Parsimonious-Mixtures-of-Dimension-Wise-Scaled-Normal-Mixtures

This repository contains the code for fitting parsimonious mixtures of Dimension-Wise Scaled Shifted Exponential Normal Mixtures.
In the following, you find a description of the main functions (and their arguments) contained in the Main.R and Utils.R files.

## DSNM_M.fit (Main.R) ##

### Description ###

Fits, by using EM-based algorithms, parsimonious mixtures of Dimension-Wise Scaled Shifted Exponential Normal Mixtures. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

DSNM_M.fit (X, k, corr = "all", scale = "all", tailedness = "all", rel.tol = 0.001, iter.max = 1000, Asympt = 100, verbose = TRUE, nThreads = 1)

### Arguments ###


