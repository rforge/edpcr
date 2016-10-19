rcq<-function(E=.95,E1=.95,mu=1.5,disp=1.2,N=10^4,nx=70,ny=11,h=.1,A=10^-10,c.mean=NULL,trendx=NULL,trendy=NULL,incxy=TRUE) {
  if (!is.null(nx)&!is.null(ny)) {incxy<-TRUE;N<-nx*ny} else incxy<-FALSE
  lam<-mu2lam(mu,disp)
  p0<-compoisson::dcom(0,lam,disp)
  if (!is.null(c.mean)) {
    K<-mu/(1-p0)*(1+E1)*(1+E)^(c.mean-2)
    A<-h/K
  } else {
    K<-h/A
  }
  x0s<-compoisson::rcom(N,lam,disp)
  n0<-sum(x0s==0)
  Cqs<-sapply(x0s,function(x0) {
    if (x0==0) return(NA)
    x<-c(x0+rbinom(1,x0,E1),x0)
    i<-0
    while(x[1]<K)  {
      i<-i+1
      if (x[1]>.Machine$integer.max) {
        n1<-ceiling(x[1]/.Machine$integer.max)
        x<-c(x[1]+sum(as.numeric(rbinom(n1,x[1]%/%n1,E)))+rbinom(1,x[1]%%n1,E),x[1])
      } else x<-c(x[1]+rbinom(1,x[1],E),x[1])}
    cq<-i+log(K/x[2])/log(x[1]/x[2])
  })
  if (incxy) {
    cqs<-cbind(x=rep(1:nx,each=ny),y=rep(1:ny,nx),cq=Cqs)
    if (!is.null(trendx)&!is.null(trendy)) cqs[,"cq"]<-apply(cqs,1,function(w) if (is.na(w["cq"])) NA else w["cq"]+(w["x"]-nx/2)*trendx+(w["y"]-ny/2)*trendy)
  }
  return(cqs)
}