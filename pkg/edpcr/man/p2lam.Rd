\name{p2lam}
\alias{p2lam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate \eqn{\lambda} of Conway Maxwell Poisson distribution from probability of 0 sucesses
}
\description{
Uses \code{\link{uniroot}} to calculate lambda of the Conway Maxwell Poisson distribution with given probability of 0 successes and dispersion.
}
\usage{
p2lam(p0, disp, n = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p0}{
probability of 0 successes.
}
  \item{disp}{
dispersion parameter (positive).
}
  \item{n}{
truncation integer.
}
}
\details{
For the Conway Maxwell Poiosson ditribution the probability of \eqn{x} successes is given by
\deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{(x!)^\nu}}{p(x) = \lambda^x exp(-\lambda)/(x!)^\nu}
for \eqn{x = 0, 1, 2, \ldots}.
 
\code{p2lam} gives \eqn{\lambda}{\lambda} for the Conway Maxwell Poisson distribtution with \eqn{p(0)} and \eqn{\nu}{\nu}.
 
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
\code{\link{mu2lam}}, \code{\link{p2mu}}
}
\examples{
p2lam(.05,1.2)
}