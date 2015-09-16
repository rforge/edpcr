\name{detrend}
\alias{detrend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detrend Cq data
}
\description{
Removes x-trend and/or y-trend from Cq data using \code{\link{lmRob}}.
}
\usage{
detrend(data = NULL, detr = "xy", 
  nx = NULL, ny = NULL, inc.bs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
\code{cqdat} object or data frame
}
  \item{detr}{
character string: if "x", "y" or "xy" remove x-trend, y-trend or both respectively.
}
  \item{nx, ny}{
numbers of rows and columns. If NULL, maxima of the x- and y-locations from \code{data} are used.
}
  \item{inc.bs}{
logical, if \code{TRUE} include removed trends in result.
  }
}
\value{
If \code{inc.bs} is \code{FALSE}, the matrix of detrended data.\cr

If \code{inc.bs} is \code{TRUE} then a list containing:
\item{data}{the detrended data}
\item{b.x}{the removed trend for x}
\item{b.y}{the removed trend for y}
}
\author{
Philip Wilson
}
\examples{
dat<-fetch(Exp37a,1)

dat2<-detrend(dat,inc.bs=TRUE)
dat2[c("b.x","b.y")]
}