#--- somehow Xptr doesn't get deleted properly ---------------------------------

require(msde)
require(GaussCop)
setwd("~/Documents/proj/msde/msdeHeaders/tests/dontrun")

# compile model
replicate(10, {
  param.names <- c("alpha", "gamma", "beta", "sigma", "rho", "lambda")
  data.names <- c("X", "Z")
  #hmod <- NULL
  #if(exists("hmod")) {
    #rm(hmod)
    #Sys.sleep(1)
    #gc()
  #}
  hmod <- sde.make.model(ModelFile = "ghestModel.h",
                         param.names = param.names,
                         data.names = data.names)
  #Hmod <- list(ptr = .sde_MakeModel())
  nRV <- 7
  nsamples <- 1e5
  X <- matrix(rnorm(nsamples*nRV), nsamples, nRV)
  gcop <- gcopFit(X = X, fitXD = "normal")
})

hmod <- NULL
hmod <- list(ptr = .sde_MakeModel())

nRV <- 7
nsamples <- 1e5
X <- matrix(rnorm(nsamples*nRV), nsamples, nRV)
gcop <- gcopFit(X = X, fitXD = "normal")
