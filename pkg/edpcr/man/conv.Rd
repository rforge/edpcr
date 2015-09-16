\name{conv}
\alias{conv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert between parameterisations
}
\description{
This function is used by \code{cqmc} to convert between parameterisations. Optionally includes determinant required by \code{cqmc}.
}
\usage{
conv(pars, from = "m", to = "c", det = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
vector of parameters with names 
}
  \item{from}{
"m", "p" or "c" represents initial parameterisation
}
  \item{to}{
"m", "p" or "c" represents final parameterisation
}
  \item{det}{
logical, if TRUE then output is list with determinant of transformation
}
}
\details{
"c" and cpars represent calculation pars (includes lambda used in CMP distribution).\cr
"p" and pars represent original pars (includes mu, parameter of interest).\cr
"m" and mpars represent pars transformed to find mode using optim and to run MCMC.\cr
}
\value{
If det is FALSE then transformed vector.

If det is TRUE then a list with:
\item{new.pars}{The transformed vector}
\item{dfac}{The absolute value of the determinant of the transformation if from is "m" and to is "c" and NA otherwise}
}
\author{
Philip Wilson
}
\seealso{
\code{\link{conv.mcmc}}
}
\examples{
pars<-c(mu=1.5,E=0.95,A=10^-15,trendx=10^-3,
 trendy=-10^-3,E1=0.8,disp=1.5)

mpars<-conv(pars,"p","m")

cpars<-conv(mpars,"m","c",det=TRUE)

conv(cpars[[1]],"c","p")
}