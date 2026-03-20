README.txt
===========

This supplementary folder contains R scripts used to compute and evaluate
confidence and credible intervals for prevalence estimation under binomial
and hypergeometric sampling frameworks, as reported in the associated article.

File descriptions
-----------------

prevCI_hypgeom_exact.r
    R script implementing an exact confidence interval calculation based on
    the hypergeometric distribution. The function computes exact prevalence
    confidence intervals without relying on asymptotic approximations.

prevCI_hypgeom_simul.r
    R script implementing an exact-type hypergeometric confidence interval
    using Monte Carlo simulation to approximate the interval limits.

prevCI_hypgeom_Bayes.r
    R script implementing a simple analytic Bayesian highest density
    credible interval for prevalence under a hypergeometric model.
    The implementation is fully analytic and does not rely on MCMC.

prevCI_Ge2024.R
    R script implementing the Bayesian prevalence interval proposed by
    Ge et al. (2024), reproduced here for comparison.

prevCI_Ge2024_mod.R
    R script implementing the Bayesian prevalence interval proposed by
    Ge et al. (2024), modified to truncate the limits of the interval,
    as per the suggestion of an anonymous reviewer.

BinomialConfIntSeSp.R
    R script implementing an exact binomial confidence interval that
    accounts for imperfect sensitivity and specificity, following
    Reiczigel et al. (2010)

prevCI_hypgeom_prLik.r
    R script implementing a profile likelihood–based confidence interval
    for prevalence under a hypergeometric model, incorporating estimates
    of sensitivity (Se) and specificity (Sp).

benchmark.R
    R script for benchmarking and comparing the runtime performance of
    the prevalence interval functions included in this folder.

simulation_1_exact.R
    R script performing simulation studies to evaluate coverage
    properties of the exact, simulation-based, and Bayesian interval
    methods.

simulation_2_prL.R
    R script performing simulation studies to evaluate coverage properties
    of the profile likelihood–based confidence interval.

Requirements
------------
All scripts are written for the R statistical computing environment.
Unless otherwise indicated within the scripts, only base R functions
are required.

Usage
-----
Scripts can be loaded into an R session using, for example:

    source("prevCI_hypgeom_exact.r")

Simulation and benchmarking scripts are intended to be run as standalone
files. Please refer to comments within individual scripts for details on
function arguments, outputs, and example usage.

Notes
-----
These scripts are provided to support transparency and reproducibility of
the results presented in the associated article. The code is supplied
as-is, without warranty.
