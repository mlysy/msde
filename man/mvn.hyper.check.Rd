% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvn.hyper.check.R
\name{mvn.hyper.check}
\alias{mvn.hyper.check}
\title{Argument checking for the default multivariate normal prior.}
\usage{
mvn.hyper.check(hyper, param.names, data.names)
}
\arguments{
\item{hyper}{The normal prior's hyperparameters: \code{NULL}, or a list with elements \code{mu} and \code{Sigma}, corresponding to a named mean vector and variance matrix (see Details).}

\item{param.names}{Vector of parameter names (see Details).}

\item{data.names}{Vector of data names (see Details).}
}
\value{
A list with the following elements:
\describe{
\item{\code{mean}}{The mean vector.}
\item{\code{cholSd}}{The upper upper Cholesky factor of the variance matrix.}
\item{\code{thetaId}}{The index of the corresponding variables in \code{theta}.}
\item{\code{xId}}{The index of the corresponding variables in \code{x0}.}
}
}
\description{
Argument checking for the default multivariate normal prior.
}
\details{
This function is not meant to be called directly by the user, but rather to parse the hyper-parameters of a default multivariate normal prior distribution to be passed to the C++ code in \code{\link[=sde.prior]{sde.prior()}} and \code{\link[=sde.post]{sde.post()}}.  This default prior is multivariate normal on the elements of \verb{(theta, x0)} specified by each of \code{names(mu)}, \code{rownames(Sigma)}, and \code{colnames(Sigma)}.  The remaining components are given Lebesgue priors, or a full Lebesgue prior if \code{hyper == NULL}.  If the names of \code{mu} and \code{Sigma} are inconsistent an error is thrown.
}
