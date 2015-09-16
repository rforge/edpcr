prep <-
function(data,failcode=NA,cq.xlim=NA,out=c(3.5,5),h=NA) {
    if (is.null(dim(data))) {
        data<-as.data.frame(cbind(x=0,y=1:length(data),data))
        data0<-data
        incxy<-FALSE
        ny<-0
        nx<-0
    } else {
        data<-data[order(data[,"x"],data[,"y"]),]
        row.names(data)<-NULL
        data0<-data<-data.frame(data)
        nx<-max(data$x)
        ny<-max(data$y)
        incxy<-TRUE
    }
    if (!is.na(failcode)) data$cq[data$cq==failcode]<-NA
    if (all(is.na(data0$cq))) {print("No finite cq");return(c())}
    if (length(cq.xlim)==2&incxy&!identical(cq.xlim,c(1,nx))&(!is.null(cq.xlim))) {
        keep<-data$x%in%seq(cq.xlim[1],cq.xlim[2])
        nt0<-sum( is.na(data[!keep,"cq"]))
        nt1<-sum(!is.na(data[!keep,"cq"]))
        data<-data[keep,]
    } else {nt0<-0;nt1<-0}
    n0<-sum(is.na(data$cq))
    n1<-sum(!is.na(data$cq))
    data<-data[!is.na(data$cq),]
    if (!identical(out,FALSE)&any(out<Inf)) {
        if (incxy) tmp<-resid(robust::lmRob(cq~x+y,data=data)) else tmp<-data$cq-mean(data$cq)
        id.hi<-which(tmp>out[1])
        id.lo<-which(tmp<(-out[2]))
        n.hi<-length(id.hi)
        n.lo<-length(id.lo)
        n1<-n1-n.hi-n.lo
        if ((n.hi+n.lo)>0) data<-data[-c(id.hi,id.lo),]
    } else {n.hi<-n.lo<-0}  
    
    tmp<-list(data=data,h=h,n0=n0,n1=n1,nt0=nt0,nt1=nt1,n.lo=n.lo,n.hi=n.hi,nx=nx,ny=ny)
    class(tmp)<-"cqdat"
    invisible(tmp)
}
