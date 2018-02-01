1.  Finish testthat test for `sde.pf`.  See why there is a discrepancy, perhaps due to invalid proposals?  See `validx(x, theta)` in e.g., `test-hest_sd.R`.  Also, randomize inputs as much as possible (can do for other tests as well).

2.  Make `sde.pf` usable.  First of all, no need to return *all* `X/lwgt` output.  Last observation is sufficient.  Can potentially add a logical argument to output if desired.  Second, don't always provide `Z` matrix.  In fact, `sdeFilter` constructor can overload without `Z`, which requires significantly less memory (i.e., only enough for one obs), i.e., much faster if we rerun `sde.pf` over and over, realllocating memory every time.  Also, provide **R**-level argument to `smc::sampler` RESAMPLE type (no resampling + prespecified `Z` is *mainly* for debugging, though does have other uses).  See RESAMPLE argument in `RcppSMC`, as it explains how/when to resample.

3.  Example in documentation.
