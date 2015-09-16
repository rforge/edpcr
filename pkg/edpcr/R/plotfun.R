plotfun <-
function(pars=NULL,h=.1,xlim=NULL,dxlim=c(-3,2.5),new=TRUE,npoints=10^3,c0=6,maxn0=7,inc0=TRUE,type=0,
                  xlab=expression("C"["q"]),ylab="density",col="black",...) {
    if (is.list(pars)) {
        if (type=="i") if ("pars0"%in%names(pars)) pars<-pars[["pars0"]] else stop("pars0 not available")
        if (type=="m") if ("pars.mod"%in%names(pars)) pars<-pars[["pars.mod"]] else stop("pars.mod not available")
        if (type=="p") if ("pars.mc"%in%names(pars)) pars<-pars[["pars.mc"]] else stop("pars.mc not available")
        if ("h"%in%names(pars)) h<-pars[["h"]]
    }
    pars<-pars[!names(pars)%in%c("trendx","trendy")] # plot of function excludes trends
    mpars<-conv(pars,"p","m")
    
    if (!new) {u<-par("usr"); xlim<-u[1:2]} else 
        if (is.null(xlim)) xlim<--(log(pars["A"])-log(h))/log(1+pars["E"])+if (is.null(dxlim)) c(-3,2.5) else dxlim
    xs<-data.frame(xlim[1]+diff(xlim)*1:npoints/npoints)
    names(xs)<-"cq"
    ys<-exp(c(llcq(mpars,xs,sum.=FALSE,full=FALSE,h=h,c0=c0,maxn0=maxn0,prior=FALSE))+500)
    ys<-ys/mean(ys)/diff(xlim)
    d<-cbind(xs,ys)
    if (new) {plot(d,type="l",col=col,xlab=xlab,ylab=ylab,...); u<-par("usr")} else points(d,type="l",col=col,...)
    if (inc0) {
        if ("disp"%in%names(pars)) {
            disp<-pars["disp"]
            lam<-mu2lam(pars["mu"],disp)
            p0<-1/compoisson::com.compute.z(lam,disp) } else 
                p0<-exp(-pars["mu"])
            lines(x=rep(u[1]+.01*diff(u[1:2]),2),y=.99*u[4]*c(0,p0),col=col)
            lines(x=rep(u[1]+.01*diff(u[1:2]),2),y=.99*u[4]*c(p0,1),col="grey")
    }
}
