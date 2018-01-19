#' particle filter prototype for eou model
#' 
#' @export
eou.pf <- function(init, npart, Z) {
    ndims <- 2 # for eou model
    nparams <- 5
    # initial values
    if(class(init) != "sde.init") {
        # all argument checking done at construction time
        stop("init must be an sde.init object.")
    }
    # dt <- init$dt.m
    # par.index <- init$nvar.obs.m
    # init.data <- init$data
    # init.params <- init$params
    ncomp <- nrow(init$data)
    # normal draws
    if(missing(Z)) Z <- matrix(rnorm(npart*ndims*(ncomp-1)), ncomp-1, npart*ndims)
    ans <- .pf_eval(initParams = init$params, initData = t(init$data),
                    dT = init$dt.m, nDimsPerObs = init$nvar.obs.m,
                    NormalDraws = t(Z))
    ans$X <- t(ans$X)
    ans$lwgt <- t(ans$lwgt)
    return(ans)
}