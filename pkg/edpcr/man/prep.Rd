\name{prep}
\alias{prep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create cqdat object
}
\description{
Creates cqdat object from data frame containing x-location, y-location and Cq values, or from a vector just containing Cq values
}
\usage{
prep(data, failcode = NA, cq.xlim = NA,
  out = c(3.5, 5), h = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data frame of Cq data with locations containing colnames "cq","x" and "y", or a vector of Cq data.
}
  \item{failcode}{
either NA, number, character string: designation of negative partitions within cq column
}
  \item{cq.xlim}{
if NA then no columns trimmed, if vector then \code{cq.xlim} is the start and end columns for which Cq data is retained, with positive and negative counts only kept for remaining trimmed columns.
}
  \item{out}{
vector of 2 positive numbers defining distance from mean (after detrending) for excluding low and high outliers respectively
}
  \item{h}{
threshold value
}
}
\details{
Where there are spurious trends in the Cq data at the edges of the panel they may be considered unreliable. The \code{cq.xlim} argument permits columns to be censored where the partitions with Cq values are treated as having at least 1 molecule. The number of such partitions is code{nt1}, with the negative partitions from the censored columns being code{nt0}.

There may also be values that are considered outliers for being too low or too high. Through the \code{out} these outliers can be censored, with the associated counts \code{n.lo} and code{n.hi} respectively. They are considered to represent more than 1 initial molecule and exactly 1 initial molecule respectively.
}
\value{
A list with the class "\code{cqdat}" and the components:
\item{data}{matrix of Cq data, excluding negative partitions}
\item{h}{threshold value}
\item{n0, n1, nt0, nt1, n.lo, n.hi}{counts}
\item{nx, ny}{numbers of columns and rows}
}
\author{
Philip Wilson
}
\seealso{
\code{\link{fetch}}
}
\examples{
dat<-cbind(x=rep(1:4,each=3),
           y=rep(1:3,times=4),
           cq=c(22.859,21.638,24.152,999,   
                29.081,23.817,23.930,23.943,
                22.934,23.473,999,22.819))
cqdata<-prep(dat,failcode=999,cq.xlim=c(1,2),h=.05)
}