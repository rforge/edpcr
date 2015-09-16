\name{plotden}
\alias{plotden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density plot of Cq data
}
\description{
Produces density plot of Cq data (possibly from cqdat object). Can include line representing proprtion of negative partitions.
}
\usage{
plotden(data, detr = "xy", new = TRUE, bw = 0.05, N = 770, 
  inc0 = TRUE, col = "black", main = "Density plot of Cq data", 
  xlab = expression("C"["q"]), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
either cqdat object, vector of Cq values, or matrix/dataframe with columns named "x", "y" and "cq" cotaining x, y and Cq values.
}
  \item{detr}{
character string indicating which variables x, y, both or neither to detrend.
}
  \item{new}{
logical, if FALSE, adds plot to existing plot.
}
  \item{bw}{
the smoothing bandwidth  used by \code{density}.
}
  \item{N}{
total number of partitions, used to determine fraction of negative parititions if \code{inc0} is \code{TRUE}.
}
  \item{inc0}{
logical, if \code{TRUE} then plot includes line representing fraction of negative partitions.
}
  \item{col}{
positive integer or character string indicating line colour.
}
  \item{main}{
plot title.
}
  \item{xlab}{
x-axis label.
}
  \item{\dots}{
other plotting parameters.
}
}
\value{
A list with components:
\item{density}{an object with class "\code{density}" (see \code{density}).}
\item{p0}{proportion of negative partitions.}
%% ...
}
\author{
Philip Wilson
}
\seealso{
\code{\link{plotall}}
}
\examples{
dat<-fetch(Exp37a,1)
plotden(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}