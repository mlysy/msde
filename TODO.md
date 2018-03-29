All tasks have been done.

Keep track of TODO lists

1.  Finish testthat test for `sde.pf`.  See why there is a discrepancy, perhaps due to invalid proposals?  See `validx(x, theta)` in e.g., `test-hest_sd.R`.  Also, randomize inputs as much as possible (can do for other tests as well).

2.  Make `sde.pf` usable.  First of all, no need to return *all* `X/lwgt` output.  Last observation is sufficient.  Can potentially add a logical argument to output if desired.  Second, don't always provide `Z` matrix.  In fact, `sdeFilter` constructor can overload without `Z`, which requires significantly less memory (i.e., only enough for one obs), i.e., much faster if we rerun `sde.pf` over and over, realllocating memory every time.  Also, provide **R**-level argument to `smc::sampler` RESAMPLE type (no resampling + prespecified `Z` is *mainly* for debugging, though does have other uses).  See RESAMPLE argument in `RcppSMC`, as it explains how/when to resample.

3.  Example in documentation.

4.  augment this to sde.pf, i.e, any sdeModel.
    - make sdePF <= particleEval  a virtual function of sdeCobj/sdeRobj
    - every call to sde.pf allocates/deallocates full memory.
    - sdeSMC.cpp should become (possible multiple) header files.

5.  add smc debug to each model using msde-test_debug.R

6.  testthat check deterministic Z input, history = TRUE/FALSE.

7.  Test a PMCMC written in R.
    - You can check against a so-called multivariate Ornstein-Uhlenbeck model, for which Kalman filter analytically does filter calculation.  This is coded in \code{mou.loglik}.
    - Compare speed/accuracy against \code{sde.post} using true posterior as
      benchmark.
    - Write up a report in the form of a vignette, which also shows how to
      use \code{sde.pf}.

8. finish debugging

9. Z should be passed as 3d array at R level, and data output should also be a 3d array.
