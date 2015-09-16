plotall <-
function(res,types="dimp",detpars=NULL,leg=TRUE,main=NULL,lwds=c(1,1,1,1),bw=.025,
                  c0=6,maxn0=7,cols=c("blue","orange","purple","red"),xlim=NULL,dxlim=c(-3.5,2),
                  inc0=TRUE,legx="topright",legy=NULL,inset=.025,detr=NULL,...) {
    if (!is.null(detpars)) if (detpars[1]=="m") detpars<-res[["pars.mod"]] else if (detpars[1]=="p") detpars<-res[["pars.mc"]] 
    types<-unlist(strsplit(types,""))
    h<-res[["h"]]
    new<-TRUE
    if (("trendx"%in%names(detpars)|"trendy"%in%names(detpars))&!is.null(detr))  stop("Cannot detrend through both arguments detpars and detr.")
    if ("trendx"%in%names(detpars)) res[["data"]]["cq"]<-res[["data"]]["cq"]-detpars["trendx"]*(res[["data"]]["x"]-res[["nx"]]/2)
    if ("trendy"%in%names(detpars)) res[["data"]]["cq"]<-res[["data"]]["cq"]-detpars["trendy"]*(res[["data"]]["y"]-res[["ny"]]/2)
    ids<-c()
    if ("d"%in%types) {
        ids<-c(ids,1)
        plotden(res[["data"]],col=cols[1],new=new,lwd=lwds[1],main="",bw=bw,inc0=inc0,xlab=expression("C"["q"]),ylab="Density",detr=detr,...)
        new<-FALSE
    }
    if ("i"%in%types) {
        ids<-c(ids,2)
        plotfun(pars=res[["pars0"]],h=h,c0=c0,maxn0=maxn0,new=new,col=cols[2],main="",lwd=lwds[2],inc0=inc0,xlab=expression("C"["q"]),xlim=xlim,dxlim=dxlim,...)
        new<-FALSE
    }
    if ("m"%in%types&!is.null(res[["pars.mod"]])) {
        ids<-c(ids,3)
        plotfun(pars=res[["pars.mod"]],h=h,c0=c0,new=new,col=cols[3],main="",lwd=lwds[3],inc0=inc0,xlim=xlim,dxlim=dxlim,...)
        new<-FALSE
    }
    if ("p"%in%types&!is.null(res[["pars.mc"]])) {
        ids<-c(ids,4)
        plotfun(pars=res[["pars.mc"]],h=h,c0=c0,new=new,col=cols[4],lwd=lwds[4],main="",inc0=inc0,xlim=xlim,dxlim=dxlim,...)
    }
    if (leg) {
        legend(x=legx,y=legy,legend=c("Data","Initial","Mode","Posterior")[ids],col=cols[ids],lwd=lwds[ids],inset=inset)
    }
    title(main=main)
}
