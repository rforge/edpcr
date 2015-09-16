llcq <-
function(mpars,data=fetch(2,1),maxn0=7,c0=6,extra=c("trendx","trendy","disp","E1"),
 h=.1,n0=300,n1=465,nt0=0,nt1=0,n.lo=0,n.hi=0,nx=0,ny=0, full=TRUE,
 delta=1e-5,sum.=TRUE,report=0,probreport=0,
 prior=FALSE, mu.fun=function(mu) 1, A.fun=function(A) 1, E.fun=function(E) 1,E1.fun=function(E1) 1,
 trendx.fun=function(x) 1, trendy.fun=function(y) 1, disp.fun=function(d) 1,prior.funs=NULL,invisible=FALSE) {
  if (report>2) cat("\n",mpars)
  if (is.null(names(mpars))) names(mpars)<-c("mu","E","A",extra)
  if (!sum.) prior<-FALSE
  if (prior) {
    if (!is.null(prior.funs)) {
      mu.fun<-prior.funs[["mu.fun"]]
      A.fun<-prior.funs[["A.fun"]]
      E.fun<-prior.funs[["E.fun"]]
      E1.fun<-prior.funs[["E1.fun"]]
      trendx.fun<-prior.funs[["trendx.fun"]]
      trendy.fun<-prior.funs[["trendy.fun"]]
      disp.fun<-prior.funs[["disp.fun"]]
    }
    if (is.null(mu.fun)) mu.fun<-function(mu) 1
    if (is.null(A.fun)) A.fun<-function(A) 1
    if (is.null(E.fun)) E.fun<-function(E) 1
    if (is.null(E1.fun)) E1.fun<-function(E1) 1
    if (is.null(trendx.fun)) trendx.fun<-function(x) 1
    if (is.null(trendy.fun)) trendy.fun<-function(y) 1
    if (is.null(disp.fun)) disp.fun<-function(disp) 1
  }
  if (class(data)=="cqdat") {
    h<-data$h;nt0<-data$nt0;nt1<-data$nt1;n0<-data$n0;n1<-data$n1
    n.hi<-data$n.hi;n.lo<-data$n.lo;nx<-data$nx;ny<-data$ny
    cqdata<-data$data
  } else cqdata<-data
  if (is.null(dim(cqdata))) {cqdata<-data.frame(cqdata);names(cqdata)<-"cq";detr<-FALSE;print("not detrending")} else detr<-TRUE

  tmp<-conv(mpars,"m","c",det=TRUE)
  cpars<-tmp[[1]]
  dfac<-tmp[[2]]
  if (prior) mu<-exp(mpars["mu"]) # convert from transformed mu
  lam<-cpars["lam"]
  E<-cpars["E"]
  if ("E1"%in%names(cpars)) E1<-cpars["E1"] else E1<-E
  A<-cpars["A"]
  if (detr&"trendx"%in%names(cpars)) trendx<-cpars["trendx"] else trendx<-0
  if (detr&"trendy"%in%names(cpars)) trendy<-cpars["trendy"] else trendy<-0
  if ("disp"%in%names(cpars)) disp<-cpars["disp"] else disp<-1
  if (trendx!=0) cqdata$cq<-cqdata$cq-trendx*(cqdata$x-nx/2)
  if (trendy!=0) cqdata$cq<-cqdata$cq-trendy*(cqdata$y-ny/2)
  cs<-cqdata$cq # cq's for detects
  if (disp==1) ps<-dpois(0:maxn0,lam) else ps<-compoisson::dcom(0:maxn0,lam,disp)
  ps<-ps/sum(ps)
  p0<-ps[1]
  ps<-ps[-1]
  for (j in 1:c0) {if (j==1) Ej<-E1 else Ej<-E;ps<-rowSums(sapply(1:length(ps),function(i) 
    c(rep(0,i-1),ps[i]*dbinom(0:i,i,Ej),rep(0,2*(length(ps)-i)))))}
  m0<-(1+E)^(cs-c0)
  v0<-(1-E)*(1+E)^(cs-c0-1)*((1+E)^(cs-c0)-1)
  m1<-(1+E)^(cs+delta-c0)
  v1<-(1-E)*(1+E)^(cs+delta-c0-1)*((1+E)^(cs+delta-c0)-1)  
  tmp<-sapply(1:length(cs),function(i) {
    z0<-(h/A-(1:(maxn0*2^c0))*m0[i])/((1:(maxn0*2^c0))*v0[i])^.5
    z1<-(h/A-(1:(maxn0*2^c0))*m1[i])/((1:(maxn0*2^c0))*v1[i])^.5
    if (is.na(tmp2<-sum(ps*(pnorm(z0)-pnorm(z1))))) return(0)
    if (tmp2>0) tmp2 else if (full) {  
      tmp3<-sapply(1:(maxn0*2^c0),function(j) 
       if (z0[j]>0&z1[j]>0) (pnorm(-z1[j])-pnorm(-z0[j])) else (pnorm(z0[j])-pnorm(z1[j])))
      return(sum(ps*tmp3))
    } else return(0)
  })/delta
  if (report>0) {if (any(tmp<=0)) cat(".") else if (runif(1)<probreport) cat("'")}
  
  if (sum.) tmp2<-(nt0+n0)*log(p0)+(nt1)*log(1-p0)+n.lo*log(1-p0)+n.hi*log(ps[1])+sum(log(tmp))-log(dfac) else
  tmp2<-log(tmp)-log(dfac)
  
  if (prior) {
    lprior<-log(mu.fun(mu))+log(E.fun(E))+log(A.fun(A))
      if (detr&"trendx"%in%names(cpars)) lprior<-lprior+log(trendx.fun(trendx))
      if (detr&"trendy"%in%names(cpars)) lprior<-lprior+log(trendy.fun(trendy))
      if ("disp"%in%names(cpars)) lprior<-lprior+log(disp.fun(disp))
      if ("E1"%in%names(cpars)) lprior<-lprior+log(E1.fun(E1))
    tmp2<-tmp2+lprior
  }
  ret<-tmp2
  if (report>1) cat("   llik",ret)
  ret<-ifelse(ret==-Inf,-10^10,ret)
  
if (invisible) invisible(ret) else return(ret)
}