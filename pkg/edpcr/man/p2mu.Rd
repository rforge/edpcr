\name{p2mu}
\alias{p2mu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate mean of Conway Maxwell Poisson distribution from probability of 0 sucesses
}
\description{
Uses \code{\link{p2lam}} to calculate \eqn{\lambda}{lambda} of Conway Maxwell Poisson distribution, then uses that to calculate the mean.
}
\usage{
p2mu(p0, disp, n = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p0}{
probability of 0 successes.
}
  \item{disp}{
dispersion parameter \eqn{\nu} (positive).
}
  \item{n}{
truncation integer.
}
}
\details{
For the Conway Maxwell Poiosson ditribution the probability of \eqn{x} successes is given by
\deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{(x!)^\nu}}{%
p(x) = \lambda^x exp(-\lambda)/(x!)^\nu}
for \eqn{x = 0, 1, 2, \ldots}.
 
\code{p2mu} gives the mean for the Conway Maxwell Poisson distribtution with \eqn{p(0)} and \eqn{\nu}{\nu} 
}
\value{
The value of \eqn{\mu}.
}
\references{
Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., "A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution," J. Royal Statist. Soc., v54, pp. 127-142, 2005.}
\author{
Philip Wilson
}

\seealso{
\code{\link{mu2lam}}
}
\examples{
p2mu(.05,1.2)
}