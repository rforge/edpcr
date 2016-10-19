fetch <-
function(data=NULL,panel=NA,dir=NULL,sep="\t",skip=11,header=TRUE,check.names=FALSE,cq.xlim=NA,out=c(3.5,5),
                cols=NULL,vals=NULL,Rep=NA,
                colx="Chamber ID",coly="Chamber ID",colpan="Panel ID",colcq="Value",
                colh="Threshold",colfail="Call",
                funx=function(x) substr(x,6, 7),funy=function(x) substr(x,10,11),
                funfail=function(x) x=="Fail",h=NA,raw=FALSE) {
  if (is.data.frame(data)) {
    if (header) colnames(data)<-as.matrix(data[skip+1,])
    data<-data[-(1:(skip+if (header) 1 else 0)),]
  }
    if (is.character(data)) {
        dir.orig<-getwd()
        setwd(dir)
        data<-read.table(data,skip=skip,sep=sep,header=header,check.names=TRUE)
        setwd(dir.orig)
    }

  if (!is.na(panel)) {
      data<-data[data[,colpan]==panel,]
    } else
    {
        for (i in 1:length(cols)) data<-data[data[,cols[i]]==vals[i],]
        if (nrow(data)==0) stop("No panels remaining after subsetting")
        if (!is.na(Rep)) {
            pans<-unique(data[,colpan])
            if (Rep>length(pans)) stop(paste(sep="","Insufficient panels after subsetting for Rep = ",Rep)) else 
                data<-data[data[,colpan]==pans[Rep],]
        }
    }
    data[funfail(data[,colfail]),colcq]<-NA

    if (is.na(h)) h<-as.numeric(as.character(data[1,colh]))

    data<-cbind(x=as.numeric(sapply(data[,colx],funx)),
                y=as.numeric(sapply(data[,coly],funy)),
                cq=as.numeric(as.character(data[,colcq])))
    data<-data[order(data[,2],data[,1]),]

    if (raw) return(data) else invisible(prep(data,NA,cq.xlim=cq.xlim,out=out,h=h))
}
