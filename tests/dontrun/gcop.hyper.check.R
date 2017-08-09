gcop.hyper.check <- function(hyper, param.names, data.names) {
  nparams <- length(param.names)
  ndims <- length(data.names)
  var.names <- c(param.names, data.names)
  if(!is.null(hyper)) {
    # error checking
    if(class(hyper) != "gcop") {
      stop("hyper must be NULL or a gcop object.")
    }
    # check variable names
    vnames <- names(hyper$XDens)
    if(!msde:::is.valid.vars(vnames, var.names)) {
      stop("Variable names in hyper must be a subest of param.names and data.names.")
    }
    # indices of the variables
    var.id <- sapply(vnames, function(x) which(x == var.names))
    # order the variables
    var.ord <- order(var.id)
    hyper$XDens <- hyper$XDens[var.ord]
    hyper$Rho <- hyper$Rho[var.ord,var.ord]
    var.id <- sort(var.id)
    # separate into theta and x components
    theta.id <- var.id[var.id <= nparams]
    x.id <- var.id[var.id > nparams] - nparams
    # format arguments
    # order:
    # nBreaks, range, dx, pdf, logPdf, cdf,
    # mean, sd, RhoCholSd
    # paramId, dataId
    prior.args <- sapply(names(hyper$XDens[[1]]), function(nm) {
      lapply(hyper$XDens, function(xd) xd[[nm]])
    }, simplify = FALSE)
    #dx <- apply(matrix(prior.args$xrng,nrow = 2), 2, diff)
    #dx <- dx/prior.args$ndens
    dx <- sapply(prior.args$xrng, diff)/unlist(prior.args$ndens)
    prior.args <- c(prior.args[1:2], list(dx = dx),
                    prior.args[3:7],
                    list(Rho = chol(hyper$Rho),
                         thetaId = theta.id, xId = x.id))
    names(prior.args)[c(1,2,4,5,6,9)] <- c("nBreaks", "range", "pdf",
                                         "logPdf", "cdf", "RhoCholSd")
    prior.args <- lapply(prior.args, function(x) {
      if(length(x) == 0) NULL else as.double(unlist(x))
    })
  } else {
    prior.args <- list(NULL)
  }
  prior.args
}
