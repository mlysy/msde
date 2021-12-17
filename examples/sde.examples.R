# Heston's model
hmod <- sde.examples("hest") # load pre-compiled model

# inspect model's C++ code
hfile <- sde.examples("hest", file.only = TRUE)
cat(readLines(hfile), sep = "\n")

\dontrun{
# compile it from scratch
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = hfile,
                       param.names = param.names,
                       data.names = data.names)
}
