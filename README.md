# A unified approach to penalized likelihood estimation of covariance matrices in high dimensions

## How to install the package:
```
library(devtools)
install_github("luca-cibinel/gicf", build_vignette = F)

library(gicf)
```


## Overview
This package implements the Generalised Iterative Conditional Fitting for the optimisation of the penalised Gaussian loglikelihood

$$-\log{|\Sigma|} - \text{tr}(\Sigma^{-1}S) - \lambda\|\Sigma - \text{diag}(\Sigma)\|_1 - \kappa\|\text{diag}(\Sigma^{-1})\|_1,$$

under the constraint that $\Sigma$ satisfies a given pattern of zeros.

The package also implements some helper functions which allow to compute the maximum value of the parameters $\kappa$ and $\lambda$ for which the solution is not trivial.

## Reproducibility

Inside the folder `Simulations and analysis` there are the files used to perform both the simulation studies and the analysis on the sonar data. These uses a local implementation of GICF (equivalent to the `gicf` package), contained in the folder `Simulations and analysis/gicf`. When executing the R scripts, the working directory should be set to the directory of the script.

### Simulation studies
The folder `Simulations and analysis/simulations` contains the simulated data and the R scripts of the simulation studies:
- `simulation_mle.R` compares the MLE estimate versus the ridge-regularised estimate of $\Sigma$ in under the specification of an adjacency matrix.
- `simulation_time.R` compares the computational time required to GICF and covglasso to estimate a sparse covariance matrix under a known sparsity pattern.
- `simulation_lasso.R` compares the covglasso estimate versus the ridge-regularised covglasso estimate of $\Sigma$.

The data contained in the `Simulations and analysis/simulations/data` are simulated datasets sampled from a multivariate normal distribution with mean $0$. Each dataset is described by two files:
```
sigma_mod_RB_d_[D]_p_[P]_n_[N].dat
simul_mod_RB_d_[D]_p_[P]_n_[N].dat
```
where the prefix `sigma` indicates the file which contains the true covariance matrix and the prefix `simul` contains the simulated data. In the name of each file, `D` indicates $10$ times the density of the covariance matrix, `P` indicates the number of covariates and `N` indicates the number of observations.

The folder `Simulations and analysis/simulations/environments` contains one R environment for each simulation study. If these environment are loaded, the output can be recovered by running the section "OUTPUT".

### Sonar data analysis
Inside the folder `Simulations and analysis/sonar data analysis` there is the R script which performs the analysis. The data is downloaded directly by the script.

Together with the script there are two R environments, for the banded and non-banded estimators, which contains the computed values of the cross validation objective function, used to perform model selection. If those enviornments are loaded, the output can be recoverd directly by running the section "OUTPUT".
