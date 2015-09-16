\name{get.pars0}
\alias{get.pars0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get initial Cq-dPCR parameters
}
\description{
Calculates initial parameter estimates from data. Used by \code{\link{cqmc}}.
}
\usage{
get.pars0(data, extra = c("trendx", "trendy", "disp", "E1"),
  h = 0.1, n0 = NULL, n1 = NULL, 
  nt0 = 0, nt1 = 0, n.lo = 0, n.hi = 0, 
  maxn0 = 7, E.init = 0.9, E1.init = E.init, d.init = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
object or data data frame. Contains x-locations, y-locations and Cq values including \code{NA}'s for negative partitions.
}
  \item{extra}{
vector of names of parameters to include in model in additon to mu, E, and A. Defaults to all possibitilies: \code{c("E1","trendx","trendy","disp")}.
}
  \item{h}{
threshold value. Not needed if data is cqdat object.
}
  \item{n0, n1, nt0, nt1, n.hi, n.lo}{
counts of the number of negative partitions (\code{n0}), positive partitions (\code{n1}), trimmed negative partitions (\code{nt0}), trimmed positive partitions (\code{nt1}), low outliers (\code{n.lo}) and high outliers (\code{n.hi}). Not needed if data is \code{cqdat} object.
}
  \item{maxn0}{
positive truncation of Poisson distribution in calculating initial value for A.
}
  \item{E.init, E1.init}{
initial parameter values (probabilities) for \eqn{E} and \eqn{E_1}.
}
  \item{d.init}{
initial parameter value (positive) for dispersion (\eqn{\nu}).
}
}
\details{
Initial values for \eqn{nu}, \eqn{E} and \eqn{E_1} are given by the variables \code{d.init}, \code{E.init} and \code{E1.init} respectively. Initial values for \code{trendx} and \code{trendy} come robust linear regression using \code{\link{lmRob}}. The initial value for \eqn{\mu} is derived from the proportion of negative partitions, and if \code{d.init} is 1 this will be the same as for the standard count-based method. 

The initial value of \eqn{A} is based on \eqn{A=h(1+E1)^-1*(1+E)^(1-m1)} where \eqn{m1} is an estimate of the mean Cq for 1 initial molecule and is based in turn on the sample mean of the Cq values.
}
\value{
vector of initial parameter values.
}
\author{
Philip Wilson
}
\examples{
dat<-fetch(Exp37a,1)
get.pars0(dat)
}