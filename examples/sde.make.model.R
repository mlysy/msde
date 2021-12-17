# header (C++) file for Heston's model
hfile <- sde.examples("hest", file.only = TRUE)
cat(readLines(hfile), sep = "\n")

\donttest{
# compile the model
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = hfile,
                       param.names = param.names,
                       data.names = data.names)

hmod
}
