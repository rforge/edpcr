\name{plotall}
\alias{plotall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of Cq density/densities and/or density plot of data.
}
\description{
Produces plot of any combination of density plot of data (using \code{\link{plotden}}) and Cq density for initial, mode and mean parameters (using \code{\link{plotfun}}). Can include lines representing proprtion of negative partitions.
}
\usage{
plotall(res, types = "dimp", detpars = NULL, leg = TRUE, 
  main = NULL, lwds = c(1, 1, 1, 1), bw = 0.025, 
  c0 = 6, maxn0 = 7, 
  cols = c("blue", "orange", "purple", "red"), 
  xlim = NULL, dxlim = c(-3.5, 2), inc0 = TRUE, 
  legx = "topright", legy = NULL, inset = 0.025, detr = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
\code{cqres} object
}
  \item{types}{
charcters included in this string indicate which plots to include - data density plot (d), Cq densities for initial parameters (i), mode parameters (m) and mean parameters (p). 
}
  \item{detpars}{
vector determining detrending if any of data for density plot. If it contains "trendx" and/or "trendy" then the trend(s) is/are removed from the data.
}
  \item{leg}{
logical, if \code{TRUE} then legend is included.
}
  \item{main}{
the main title.
}
  \item{lwds}{
vector line widths associated with data, initial, mode and mean plots respectively (see \code{inc} argument).
}
  \item{bw}{
bandwidth for density plot.
}
  \item{c0}{
number of cycles for which exact probabilities are caclulated. This can have a significant impact on speed of computation.
}
  \item{maxn0}{
the maximum number of initial molecules used in computation.
}
  \item{cols}{
colours of lines associated with data, initial, mode and mean plots respectively (see \code{inc} argument).
}
  \item{xlim}{
x-axis range (if density plot of data not included).
}
  \item{dxlim}{
x-axis range relative to approximate mean Cq for 1 initial molecule (if density plot of data not included and \code{xlim} is \code{NULL}).
}
  \item{inc0}{
logical, whether to include line representing fraction of negative partitions.
}
  \item{legx, legy, inset}{
x, y and inset arguments of legend function.
}
  \item{detr}{
a character string. If x and/or y are included then the data is detrended for that/those variables using \code{rlm} for the density plot.
}
  \item{\dots}{
plotting parameters.
}
}
\details{
Trends are excluded from parameter density plots, so to be comparable to a data density plot the trend(s) must be removed from the data density plot (see examples).
}
\value{
\code{NULL}
}
\author{
Philip Wilson
}
\seealso{
\code{\link{density}}
}
\examples{
# Plot Cq density for posterior mean parameters using result of cqmc.
# Remove trend from data density plot to make comparable.
cqresult<-cqmc(fetch(Exp37a,1),E.init=.9,E1.init=.5,mc.rep=0,mod.rep=200)
plotall(cqresult,types="dim",detpars=cqresult[["pars.mc"]])
}
\keyword{hplot}