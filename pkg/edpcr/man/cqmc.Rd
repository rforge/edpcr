\name{cqmc}
\alias{cqmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mode and/or MCMC sample of Cq-dPCR posterior
}
\description{
Estimates the mode and produces an MCMC sample of the posterior distribution of the parameters of the Cq-dPCR model.
}
\usage{
cqmc(data, mc.rep = 10^4, h = NULL, 
  n0 = NULL, n1 = NULL, 
  nt0 = 0, nt1 = 0, n.lo = 0, n.hi = 0, 
  pars0 = NULL, report = 1, probreport = 0.1, 
  extra = c("trendx", "trendy", "disp", "E1"), 
  c0 = 6, maxn0 = 7, 
  mod.method = "Nelder-Mead", mod.rep = 2000, 
  burnin = 0, nreport = 10, cq.xlim = NA, tune = 1, 
  E1.init = 0.8, E.init = 0.9, d.init = 1, 
  prior = TRUE, 
  mu.fun = function(x) dgamma(x, 1.5, 1.5), 
  A.fun = function(x) x^-1, 
  E.fun = function(x) dbeta(x, 60, 5), 
  E1.fun = E.fun, 
  trendx.fun = function(x) 1, trendy.fun = function(x) 1, 
  disp.fun = function(x) dgamma(x, 10, 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
cqdat object or data frame
}
  \item{mc.rep}{
number of MCMC samples. If 0 then MCMC not performed.
}
  \item{h}{
Threshold value
}
  \item{n0, n1, nt0, nt1, n.lo, n.hi}{
integer counts of the negative partitions (\code{n0}), positive partitions (\code{n1}), number of trimmed negative partitions (\code{nt0}), trimmed positive partitions (\code{nt1}), low outliers (\code{n.lo}) and high outliers (\code{n.hi}). Not needed if data is \code{cqdat} object.
}
  \item{pars0}{
initial model parameter values
}
  \item{report}{
if postive tracing information is produced. High values may produce more information.
}
  \item{probreport}{
if \code{report>0}, the probability that "." is printed when a valid log-likelihood is calculated. Useful for showing progress in finding mode or simulating MCMC chain.
}
  \item{extra}{
vector of names of parameters to include in model in additon to \code{"mu"}, \code{"E"}, and \code{"A"}. Defaults to all possibitilies: \code{c("E1","trendx","trendy","disp")}.
}
  \item{c0}{
number of cycles for which exact probabilities are caclulated. This can have a significant impact on speed of computation.
}
  \item{maxn0}{
the maximum number of initial molecules used in computation.
}
  \item{mod.method}{
character string, name of optimsation method used by \code{optim} to find posterior mode.
}
  \item{mod.rep}{
maximum number of iterations used by \code{optim} to find mode, which is used to initialise MCMC process if \code{mc.rep>0}. If 0 then mode not found.
}
  \item{burnin}{
proportion of MCMC samples discarded as burn-in.
}
  \item{nreport}{
number of times the function value and acceptance rate of the MCMC sample are printed.
}
  \item{cq.xlim}{
if \code{NA} then no rows trimmed, if vector then \code{cq.xlim} is the start and end rows for which Cq data is retained, with positive and negative counts only kept for remaining trimmed rows.
}
  \item{tune}{
the tuning parameter for the Metropolis sampling. If a vector, then the same length as the parameter vector.
}
  \item{E.init, E1.init}{
vectors of initial parameter values (probabilities) for \eqn{E} and \eqn{E_1}.
}
  \item{d.init}{
initial parameter value (positive) for dispersion (\eqn{\nu}{\nu}).
}
  \item{prior}{
logical. if \code{FALSE} use constant priors, in which case the mode is equivalent to the MLE. Otherwise use priors specified by \code{mu.fun}, \code{A.fun}, \code{E.fun}, \code{E1.fun}, \code{trendx.fun}, \code{trendy.fun} and \code{disp.fun}.
}
  \item{mu.fun, A.fun, E.fun, E1.fun, trendx.fun, trendy.fun, disp.fun}{
priors as single parameter functions.
}
}
\details{
This function can be used to find the posterior mode using \code{link{optim}} and/or simulate a posterior MCMC sample of the Enhanced dPCR model using \code{\link{MCMCmetrop1R}}. If \code{prior} is \code{FALSE} then the former is equivalent to finding the Maximum Likelihood estimate (MLE).

The mode found for a particular set of initial paramemeters may not be the global maximum. For this reason it is recommended to search from a number combinations of \eqn{E} and \eqn{E_1}{E1} set through the vectors \code{E.init} and \code{E1.init}.

If \code{report>0} then "." is printed when the posterior is calculated as 0 (due to computational limitations). If \code{report>0} and \code{probreport>0} then "'" is printed with probability \code{probreport} when a non-zero posterior is calculated. If \code{optim} fails to find non-zero posteriors, and thus fails to work then a better initial value may be required.
}
\value{
A list with components:
\item{cqdata}{the cqdata object containing the data.}
\item{counts}{a vector of the count data.}
\item{pars0}{initial parameters.}
\item{logval0}{log-posterior at \code{pars0}.}
\item{h}{threshold value.}
\item{nx, ny}{column and row numbers.}
Also included if \code{mod} is \code{TRUE}:
\item{pars.mod}{parameters at posterior mode.}
\item{logval.mod}{log-posterior at \code{pars.mod}.}
\item{mod.res}{matrix of parameters and log-posteriors found from running optim at different starting values of \code{E} and \code{E1}.}
Also included if \code{mc} is \code{TRUE}:
\item{pars.mc}{sample mean of MCMC sample of posterior.}
\item{logval.mc}{log-posterior at \code{pars.mc}}
\item{pars.sum}{summary of MCMC sample of posterior.}
\item{pars.chain}{MCMC sample of posterior.}
\item{mc.vals}{log-posteriors associated with MCMC sample of posterior}
\item{mc.acc}{acceptance rate.}
\item{mcmod.pars}{posterior mode from MCMC posterior sample.}
\item{mcmod.val}{log-posterior at \code{mcmod.pars}.}
\item{prior.funs}{list of prior distributions.}
}
\author{
Philip Wilson
}
\examples{

\dontrun{
# These examples take some time to run.

dat<-fetch(Exp37a,panel=1)
# Estimate posterior mode for single intitial parameter vector
res1<-cqmc(dat,mc.rep=0,E.init=.9,E1.init=.8,probreport=.1)

res2<-cqmc(dat,mod.rep=10^3,mc.rep=10^3,E.init=.9,E1.init=c(.5,.7,.9),probreport=.1)
# Estimate posterior mode from 9 initial parameter vectors. 
# Simulate MCMC sample from posterior distribution 
#  using estimated mode as initial value.
res3<-cqmc(dat,mc.rep=10^4,
  E.init=rep(c(.9,.7),each=3),
  E1.init=rep(c(.9,.7,.5),2),probreport=.1) 
}
}