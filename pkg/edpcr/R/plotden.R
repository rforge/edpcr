plotden <-
function(data,detr="xy",new=TRUE,bw=.05,N=770,inc0=TRUE,col="black",main="Density plot of Cq data",xlab=expression("C"["q"]),...) 
{
    if (class(data)=="cqdat") data<-data[[1]]
    if (is.matrix(data)) data<-data.frame(data)
    p0<-1-nrow(data)/N # (proportion 0 molecules)
    nx<-max(data[,"x"])
    ny<-max(data[,"y"])
    if (identical(detr,"x")) {
        trend<-(lm<-robust::lmRob(cq~x,data=data))[[1]]["x"]
        data[,"cq"]<-data[,"cq"]-trend*(data[,"x"]-nx/2)
    } else 
        if (identical(detr,"y")) {
            trend<-robust::lmRob(cq~y,data=data)[[1]]["y"]
            data[,"cq"]<-data[,"cq"]-trend*(data[,"y"]-ny/2)
        } else 
            if (identical(detr,"xy")) {
                trend<-robust::lmRob(cq~x+y,data=data)[[1]]
                data[,"cq"]<-data[,"cq"]-trend["x"]*(data[,"x"]-nx/2)-trend["y"]*(data[,"y"]-ny/2)
            }
    data<-detrend(data,detr)
    den<-density(data[,"cq"],bw=bw)
    if (new) plot(den,col=col,main=main,xlab=xlab,...) else points(den,type="l",col=col,...)
    u<-par("usr")
    if (inc0) {
        lines(x=rep(u[1]+.005*diff(u[1:2]),2),y=.99*u[4]*c(0,p0),col=col,...)
        lines(x=rep(u[1]+.005*diff(u[1:2]),2),y=.99*u[4]*c(p0,1),col="grey",...)
    }
    invisible(list(density=den,p0=p0))
}
