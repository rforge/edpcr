\name{llcq}
\alias{llcq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-posterior for Cq-dPCR data
}
\description{
Calculates the log-posterior or log-likelihood for the data and the parameters.
}
\usage{
llcq(mpars, data = fetch(2, 1), maxn0 = 7, c0 = 6, 
  extra = c("trendx", "trendy", "disp", "E1"), 
  h = 0.1, n0 = 300, n1 = 465,
  nt0 = 0, nt1 = 0, n.lo = 0, n.hi = 0, 
  nx = 0, ny = 0, 
  full = TRUE, delta = 1e-05, sum. = TRUE, 
  report = 0, probreport = 0, prior = FALSE, 
  mu.fun = function(mu) 1, A.fun = function(A) 1, 
  E.fun = function(E) 1, E1.fun = function(E1) 1, 
  trendx.fun = function(x) 1, trendy.fun = function(y) 1, 
  disp.fun = function(d) 1, prior.funs = NULL, 
  invisible = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpars}{
parameter vector in "m" parameterisation. Use \code{conv} if necessary to reparameterise.
}
  \item{data}{
\code{cqdat} object or data frame.
}
  \item{maxn0}{
the maximum number of initial molecules used in computation.
}
  \item{c0}{
number of cycles for which exact probabilities are caclulated. This can have a significant impact on speed of computation.
}
  \item{extra}{
vector of names of parameters to include in model in additon to \code{"mu"}, \code{"E"}, and \code{"A"}. Defaults to all possibitilies: \code{c("E1","trendx","trendy","disp")}.
}
  \item{h}{
threshold value  Not needed if data is a cqdat object.
}
  \item{n0, n1, nt0, nt1, n.lo, n.hi}{
counts of the number of negative partitions (\code{n0}), positive partitions (\code{n1}), trimmed negative partitions (\code{nt0}), trimmed positive partitions (\code{nt1}), low outliers (\code{n.lo}) and high outliers (\code{n.hi}). Not needed if data is \code{cqdat} object.
}
  \item{nx, ny}{
numbers of rows, columns respectively.
}
  \item{full}{
logical, if \code{TRUE} then more robust method used for small log-likelihood contributions.
}
  \item{delta}{
positive number used in finite difference method.
}
  \item{sum.}{
logical, if \code{TRUE} calculate log-likelihood for all data, otherwise calculate vector of log-densities for each Cq value.
}
  \item{report}{
if postive then tracing information is produced. High values may produce more information.
}
  \item{probreport}{
if \code{report > 0}, the probability that "'" is printed when a valid log-likelihood is calculated.
}
  \item{prior}{
logical, if \code{FALSE} use constant priors, which is equivalent to returning the log-likelihood. Otherwise use priors specified by \code{mu.fun}, \code{A.fun}, \code{E.fun}, \code{E1.fun}, \code{trendx.fun}, \code{trendy.fun} and \code{disp.fun}.
}
  \item{mu.fun, A.fun, E.fun, E1.fun, trendx.fun, trendy.fun, disp.fun}{
priors as single parameter functions.
}
  \item{prior.funs}{
list of named prior functions.
}
  \item{invisible}{
logical. If \code{TRUE} then the result is printed.
  }
}
\details{
Due to computation limits very small log-posteriors cannot be calculated and \code{-10^10} is returned instead. If \code{report>0} "." will also be printed. For other log-posteriors if \code{report>0} then "'" is printed with probability \code{probreport}.
}
\value{
If \code{sum.} is \code{TRUE} the log-posterior, otherwise a vector of log-densities for the Cq values.
}
\author{
Philip Wilson
}
\seealso{
\code{\link{cqmc}},\code{\link{plotfun}}
}
\examples{
dat<-fetch(Exp37a,1)
pars0<-get.pars0(dat)
llcq(conv(pars0,"p","m"),dat)
}