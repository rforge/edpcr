\name{conv.mcmc}
\alias{conv.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert MCMC chain between parameterisations
}
\description{
Applies \code{conv} to each row of a MCMC chain.
}
\usage{
conv.mcmc(chain, from = "m", to = "p")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain}{
mcmc object with each line parameterisation as given by \code{from}.
}
  \item{from}{
"m", "p" or "c" represents initial parameterisation.
}
  \item{to}{
"m", "p" or "c" represents final parameterisation.
}
}
\details{
"c" and cpars represent calculation pars (includes \eqn{\lambda} used in CMP distribution).\cr
"p" and pars represent original pars (includes \eqn{\mu}, parameter of interest).\cr
"m" and mpars represent pars transformed to find mode using optim and to run MCMC.
}
\value{
mcmc object.
}
\author{
Philip Wilson
}
\seealso{
\code{\link{conv}}
}