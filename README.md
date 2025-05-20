# A unified approach to penalized likelihood estimation of covariance matrices in high dimensions

This repository contains the code developed in the context of the paper [Cibinel et al. (2024)](https://arxiv.org/abs/2410.02403), which includes both the implementation of the Generalised Iterative Conditional Fitting (see below), the simulation studies and the practical analyses that have been conducted.

The algorithm is provided both in the form of the R package `gicf` and as the original R/C++ scripts created during the development of the methodology. The structure of the repository mirrors that of the R package submitted to CRAN, with all the information needed for the reproducibility of the results presented in Cibinel et al. (2024) stored in the folder `Simulations and analysis`.

## How to install the package:
To install from CRAN (**recommended**):
```
install.packages("gicf")

library(gicf)
```

To install directly from this repository:
```
library(devtools)
install_github("luca-cibinel/gicf", build_vignette = F)

library(gicf)
```

## Overview
The Generalised Iterative Conditional Fitting optimises the penalised Gaussian loglikelihood

$$-\log{|\Sigma|} - \text{tr}(\Sigma^{-1}S) - \lambda\|\Sigma - \text{diag}(\Sigma)\|_1 - \kappa\|\text{diag}(\Sigma^{-1})\|_1,$$

under the constraint that $\Sigma$ satisfies a given pattern of zeros.

The package also implements some helper functions which allow to compute the maximum value of the parameters $\kappa$ and $\lambda$ for which the solution is not trivial.

## Reproducibility

Inside the folder `Simulations and analysis` there are the files used to perform both the simulation studies and the analysis on the sonar data. These use a local implementation of the GICF algorithm (equivalent to the `gicf` package), contained in the folder `Simulations and analysis/gicf`. When executing the R scripts, the working directory should be set to the directory of the script.

### Simulation studies
The folder `Simulations and analysis/simulations` contains the simulated data and the R scripts of the simulation studies:
- `simulation_mle.R` compares the MLE estimate versus the ridge-regularised estimate of $\Sigma$ in under the specification of an adjacency matrix.
- `simulation_time.R` compares the computational time required to the GICF and the covglasso algorithms to estimate a sparse covariance matrix under a known sparsity pattern.
- `simulation_lasso.R` compares the covglasso estimate versus the ridge-regularised covglasso estimate of $\Sigma$.

The data contained in the `Simulations and analysis/simulations/data` are simulated datasets sampled from a multivariate normal distribution with mean $0$. Each dataset is described by two files:
```
sigma_mod_RB_d_[D]_p_[P]_n_[N].dat
simul_mod_RB_d_[D]_p_[P]_n_[N].dat
```
where the prefix `sigma` indicates the file which contains the true covariance matrix and the prefix `simul` contains the simulated data. In the name of each file, `D` indicates $10$ times the density of the covariance matrix, `P` indicates the number of covariates and `N` indicates the number of observations. The simulations regarding computational time are an exception to this format, due to the large amount of data used in this study. Instead of relying on a local copy of the data, new data are sampled each time, using the true covariance matrix stored in the appropriate `.dat` file. The original datasets can be made available upon request.

The folder `Simulations and analysis/simulations/environments` contains one R environment for each simulation study. If these environment are loaded, the output can be recovered by running the section "OUTPUT".

### Sonar data analysis
Inside the folder `Simulations and analysis/sonar data analysis` there is the R script which performs the analysis. The data is downloaded directly by the script.

Together with the script there are two R environments, for the banded and non-banded estimators, which contain the computed values of the cross validation objective function, used to perform model selection. If those enviornments are loaded, the output can be recoverd directly by running the section "OUTPUT".

## References
- Cibinel, L., Roverato, A., & Vinciotti, V. (2024). "A unified approach to penalized likelihood estimation of covariance matrices in high dimensions", **arXiv preprint**, [https://arxiv.org/abs/2410.02403](https://arxiv.org/abs/2410.02403)
