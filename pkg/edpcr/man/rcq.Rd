\name{rcq}
\alias{rcq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Cq data
}
\description{
Simulates Cq data by simulating molecule growth and calculating Cq values through simple algorithm.
}
\usage{
rcq(E=.95,E1=.95,mu=1.5,disp=1.2,N=10^4,nx=70,ny=11,h=.1,A=10^-10,c.mean=NULL,trendx=NULL,trendy=NULL,incxy=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
efficiency (probability).
}
  \item{E1}{
cycle 1 efficiency (probability).
}
  \item{mu}{
average number of molecules per partition (positive).
}
  \item{disp}{
dispersion parameter \eqn{\nu} (positive).
}
  \item{N}{
number of Cq trials (each results in either a Cq value or contributes to the count of negative partitions).
}
  \item{nx, ny}{
numbers of columns and rows.
}
  \item{h}{
threshold value.
}
  \item{A}{
fluorescence per molecule.
}
  \item{c.mean}{
mean Cq (alternative to specifying \code{A}).
}
  \item{trendx,trendy}{
x- and y-trends added to Cq data.
}
  \item{incxy}{
if \code{TRUE} then output includes x- and y-locations.
}
}
\details{
Simulates molecule growth where each molecule amplifies with probability E (E_1 in cycle 1). Cq values are calculated based on exponential growth between the cycles immediately before and after the threshold is crossed.
}
\value{
If \code{incxy} is \code{TRUE} a matrix including Cq values, x- and y-locations, otherwise a vector of Cq values. For negative partitions, Cq values are given as \code{NA}. 
}
\author{
Philip Wilson
}
\examples{
sim1<-rcq(.95,.65,N=1000,c.mean=25.5)
sim2<-rcq(.9,.9,nx=70,ny=11,c.mean=25.5,trendx=10^-3,trendy=10^-4)
}