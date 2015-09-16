\name{mu2lam}
\alias{mu2lam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate lambda of Conway Maxwell Poisson distribution from mean and dispersion
}
\description{
Uses \code{\link{uniroot}} to calculate lambda of the Conway Maxwell Poisson distribution with given mean and dispersion parameters.
}
\usage{
mu2lam(mu, disp, l0 = NULL, l1 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
     mean parameter (positive)
}
  \item{disp}{
     dispersion parameter (positive)
}
  \item{l0, l1}{
     initial end points of interval serached by uniroot, \code{l0 < l1}. The root does not need to lie between these end points.
}
}
\details{
For the Conway Maxwell Poiosson ditribution the probability of \eqn{x} successes is given by
\deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{(x!)^\nu}}{p(x) = \lambda^x exp(-\lambda)/(x!)^\nu}
for \eqn{x = 0, 1, 2, \ldots}.
 
\code{mu2lam} gives \eqn{\lambda}{\lambda} for the Conway Maxwell Poisson distribtution with mean \eqn{mu} and dispersion \eqn{\nu}{\nu}.
Uses \code{uniroot} to calculate \eqn{\lambda} of the Conway Maxwell Poisson distribution with given mean and dispersion.
}
\value{
The value of \eqn{\lambda}.
}
\references{
Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., "A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution," J. Royal Statist. Soc., v54, pp. 127-142, 2005.
}
\author{
Philip Wilson
}
\seealso{
\code{\link{p2lam}}
}
\examples{
mu2lam(1.5,1.2)
}