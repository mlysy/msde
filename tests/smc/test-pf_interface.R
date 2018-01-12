# SMC for multivariate SDEs.

# A basic PF for msde models has now been implemented.

# * arbitrary euler discretization is supported
# * bad proposals are _partially_ supported, in the sense that the log-weight
#   is set to -Inf but the particle is updated with an illegal value.  Thus,
#   the next observation could be e.g., full of NAs.
#   TODO: have some way of knowing that particle will be thrown out eventually
#         so no sense in updating proposal.
# * Only filtering is currently implemented (forward pass).
#   TODO: implement smoothing (backward pass).

# the next step is to interface this with the msde package.
# do this in 2 steps:
# 1.  sde.make.pf returns an XPtr to an sdeParFilt object, which has a public
#     method accepting parameters (and maybe data?)
#     Do this by adding a method to sdeCObject which returns an object.
#     Then wrap this object as XPtr and return.
# 2.  sde.pf takes the pointer and data/parameter inputs, and returns the particle cloud + log weights at the last step.
#     NOTE: because memory for calculations gets "deep copied" into SMC sampler,
#     observed data / normal draws can only be passed in once.  Otherwise, might
#     as well re-allocate all memory.


#---------------------------------------------------------------------------

# TODO: implement the function eou.pf directly into msde and formally test
# it using testthat interface.

# 1. convert pf_eval to eou.pf, and add it directly to the msde package.  Dont' need to document, but function should be available to users.  Essentially this means adding sdeSMC.cpp to msde/src, and eouModel.h to msde/inst/include.

# 2. include eouModel into sde.examples (see how hestModel is done).  you will also need to modify the files sdeExports and sde.examples.

# 3. incorporate the "with SMTC" into testthat.  to do this:

# 3.1 put smc-functions.R into msde/tests/testthat  (it has some redundancy with msde-testfunctions so you can get rid of those).

# 3.2 create a file tests/testthat/test-eou.pf.R which contains the testthat code.  this code will create emod object, then for many choices of theta, nObs, dT, etc., will (i) simulate data (ii) calculate pf in R (iii) calculate pf in C++ (iv) check that output is identical.
