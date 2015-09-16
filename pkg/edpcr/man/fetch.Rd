\name{fetch}
\alias{fetch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Cq data
}
\description{
Extracts Cq data from a data frame (optionally read from a file) and (optionally) creates \code{cqdat} object.
}
\usage{
fetch(data=NULL,panel=NA,dir=NULL,sep="\t",skip=11,
  header=TRUE,check.names=FALSE,cq.xlim=NA,out=c(3.5,5),
  cols=NULL,vals=NULL,Rep=NA,
  colx="Chamber ID",coly="Chamber ID",colpan="Panel ID",
  colcq="Value",colh="Threshold",colfail="Call",
  funx=function(x) substr(x,6, 7),
  funy=function(x) substr(x,10,11),
  funfail=function(x) x=="Fail",h=NA,raw=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
either data frame or filename.
}
  \item{panel}{
the panel number (alternatively use \code{cols}, \code{vals} and if necessary \code{Rep}).
}
  \item{dir}{
directory of file (if not working directory).
}
  \item{sep}{
if \code{data} is a filename then the field separator character (see \code{\link{read.table}}).  
  }
  \item{skip}{
the number of lines of the data file (excluding any header) to skip before beginning to read data.
}
  \item{header}{
logical. Whether the file or data frame contains the names of the variables as its first line.
}
  \item{check.names}{
logical. If \code{TRUE} and \code{data} is a filename then the names of the variables in the extracted data frame are checked to ensure that they are syntactically valid variable names.  
  }
  \item{cq.xlim}{
if \code{NA} then no panel columns are trimmed, if vector then \code{cq.xlim[1]} and \code{cq.xlim[2]} are the start and end columns for which Cq data is retained, with positive and negative counts only kept for remaining trimmed columns.
}
  \item{out}{
vector of 2 positive numbers defining distance from mean (after detrending) for excluding low and high outliers respectively.
}
  \item{cols}{
vector of indices and/or character strings identifying the panel to get the data for data (not used if pan specified).
}
  \item{vals}{
vector of indices and/or character strings (in order of cols) identifying the values taken by the data in the columns designated by cols for the desired panel (not used if \code{panel} specified).
}
  \item{Rep}{
if cols and vals reduce the data to more than 1 panel, then the Rep'th panel is used (not used if \code{panel} specified).
}
  \item{colx, coly, colpan, colcq, colh}{
integers and/or character strings identifying columns for x-location, y-location, panel number, Cq values, h, fail condition respectively.
}
  \item{colfail}{
integer and/or character string (or possibly vector of such) identifying columns used by funfail to identify negative partitions.
}
  \item{funx, funy}{
functions to extract x- and y-values from the columns specified by colx and coly.
}
  \item{funfail}{
function from entries of column(s) specified by colfail with logical as output to identify negative partitions.
}
  \item{h}{
threshold value (can also be extracted by specifying colh)
}
  \item{raw}{
logical. If \code{FALSE} data frame with x-, y-locations and Cq data is returned with no further processing (\code{cq.xlim} and \code{out} are ignored).
}
}
\details{
This function creates a data frame of x-location data, y-location data and Cq values and either returns just it (if \code{raw} is \code{TRUE}) or passes it to \code{prep} which returns it within a \code{cqdat} object.
 
The function extracts data for a single panel, which can either be identified by the parameter \code{panel}, or by specifying vectors \code{cols}, \code{vals} and if necessary \code{Rep}.
 
If the values of \code{cq.xlim} are not within the range of x-location values then Cq data will be excluded, but with the associated positive and negative counts retained instead.
}
\value{
Either a matrix of Cq data with locations (if \code{raw} is \code{TRUE}), or an object of class \code{cqdat} (by running \code{prep}).
}
\author{
Philip Wilson
}
\examples{
# Select panel directly. Trim data for columns 1-10 and 60-70
# (there are 70 columns). 
#Data is from Fluidigm output, so defaults are appropriate.
dat1<-fetch(Exp37a,panel=1,cq.xlim=c(11,60))

# Plot raw data
with(data.frame(fetch(Exp37a,panel=1,raw=TRUE)),plot(x,cq))

\dontrun{

# Extract data from a named tab-delimited file
dat<-fetch("dpcr.txt",sep"\t",panel=1)

# Extract data from a named comma separated values file
dat<-fetch("dpcr.csv",sep",",panel=1)


# Select panel indirectly. Select negative partitions
# to be those with Cq coded as 999.
data(dpcr2)
dat2<-fetch(dpcr,
  cols=c("Master Mix","Sample Name","Assay Name"),
  vals=c("GE","P2A","Jiang_16S"),
  Rep=1,colfail="Cq",funfail=function(x) x==999)
}
}