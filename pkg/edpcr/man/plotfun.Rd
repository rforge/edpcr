\name{plotfun}
\alias{plotfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of Cq distribution
}
\description{
Produces plot of the Cq distribution given the provided parameters. Can include a line representing the proprtion of negative partitions.
}
\usage{
plotfun(pars = NULL, h = 0.1, xlim = NULL, 
  dxlim = c(-3, 2.5), new = TRUE, npoints = 10^3, 
  c0 = 6, maxn0 = 7, inc0 = TRUE, type = 0, 
  xlab = expression("C"["q"]), ylab = "density", 
  col = "blue", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
parameters or list, such as \code{cqres} object produced by \code{cqmc}, that includes the parameters (see \code{type}).
}
  \item{h}{
threshold value
}
  \item{xlim}{
x-axis range
}
  \item{dxlim}{
x-axis range relative to approximate mean Cq for 1 initial molecule (used if \code{xlim} is \code{NULL}).
}
  \item{new}{
logical, if \code{FALSE}, adds plot to existing plot
}
  \item{npoints}{
number of points calculated for plotting
}
  \item{c0}{
c0 number of cycles for which exact probabilities are caclulated. This can have a significant impact on speed of computation.
}
  \item{maxn0}{
the maximum number of initial molecules used in computation.
}
  \item{inc0}{
logical, if TRUE include line representing fraction of negative partitions.
}
  \item{type}{
character specifying which vector of parameters to use if pars is a list ("i": \code{pars0}, "m": \code{pars.mod}, "p": \code{pars.mc}).
}
  \item{xlab}{
x-axis label
}
  \item{ylab}{
y-axis lable
}
  \item{col}{
colour of lines
}
  \item{\dots}{
plotting parameters
}
}
\details{
Plots the density of Cq given the specified parameters. The parameters are given by pars and if pars is a list also by type.

If the parameters include "trendx" and/or "trendy" they are ignored.
}
\value{
A list with components:
\item{x}{x values plotted.}
\item{y}{y values plotted.}
\item{p0}{probability of 0 molecules if \code{inc0} is \code{TRUE} and \code{NA} otherwise.}
%% ...
}
\author{
Philip Wilson
}
\seealso{
\code{\link{plotall}}
}
\examples{

pars0<-get.pars0(fetch(Exp37a,1))
plotfun(pars0)

# Plot Cq density for posterior mean parameters 
# using result of cqmc.
cqresult<-cqmc(fetch(Exp37a,1),E.init=.9,E1.init=.5,mc.rep=0,mod.rep=200)
plotfun(cqresult,type="m",dxlim=c(-4,3))
}
\keyword{hplot}