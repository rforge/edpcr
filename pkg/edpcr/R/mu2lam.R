mu2lam <-
function(mu,disp,l0=NULL,l1=NULL) {
    if (disp==1) return(mu)
    if (is.null(l0)) l0<-(mu^disp)/2
    if (is.null(l1)) l1<-(mu^disp)
    while((f.hi<-(mu-compoisson::com.mean(l1,disp)))>0) l1<-l1*2
    while((f.lo<-(mu-compoisson::com.mean(l0,disp)))<0) l0<-l0/2
    uniroot(function(lam,mu,disp) mu-compoisson::com.mean(lam,disp),c(l0,l1),f.lower=f.lo,f.upper=f.hi,mu=mu,disp=disp)$root
}
