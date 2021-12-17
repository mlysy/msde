# msde: Bayesian Inference for Multivariate Stochastic Differential Equations

*Martin Lysy, Feiyu Zhu, JunYong Tong, Trevor Kitt, Nigel Delaney*

---

### Description

Implements an MCMC sampler for the posterior distribution of arbitrary time-homogeneous multivariate stochastic differential equation (SDE) models with possibly latent components.  The package provides a simple entry point to integrate user-defined models directly with the sampler's C++ code, and parallelizes large portions of the calculations when compiled with OpenMP.

### Installation

To install the latest R release:
```r
install.packages("msde")
```
To install the latest development version, first install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/msde")
```

### Usage

Please see tutorial `vignette("msde-quicktut")` and several provided SDE models `vignette("msde-exmodels")`.
