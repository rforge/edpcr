get.pars0 <-
function(data,extra=c("trendx","trendy","disp","E1"),h=.1,
n0=NULL,n1=NULL,nt0=0,nt1=0,n.lo=0,n.hi=0,maxn0=7,
E.init=.9,E1.init=E.init,d.init=1) {
  if (class(data)=="cqdat") {
    h<-data$h; nt0<-data$nt0; nt1<-data$nt1; n0<-data$n0; n1<-data$n1
    n.hi<-data$n.hi; n.lo<-data$n.lo
    cqdata<-data$data
  } else cqdata<-data
  pars<-rep(NA,3+length(extra))
  names(pars)<-c("mu","E","A",extra)

  p0<-(nt0+n0)/(nt0+nt1+n0+n1+n.hi+n.lo)
  if ("disp"%in%extra) pars["disp"]<-d.init
  pars["mu"]<-p2mu(p0,d.init)
  pars["E"]<-E<-E.init
  if ("E1"%in%extra) pars["E1"]<-E1<-E1.init else E1<-E.init
  ps<-dpois(1:maxn0,-log(p0))/sum(dpois(0:maxn0,-log(p0)))
  m.adj<-1/((1-p0)*log(1+E))*(sum(sapply(1:maxn0,function(i) ps[i]*log(i))))
  m1<-mean(cqdata$cq)+m.adj
  pars["A"]<-h*(1+E1)^(-1)*(1+E)^(1-m1)
  if ("trendx"%in%extra&"trendy"%in%extra) pars[c("trendx","trendy")]<-robust::lmRob(cq~x+y,data=cqdata)[[1]][c("x","y")] else
  if ("trendx"%in%extra) pars["trendx"]<-robust::lmRob(cq~x,data=cqdata)[[1]]["x"] else   
  if ("trendx"%in%extra) pars["trendy"]<-robust::lmRob(cq~y,data=cqdata)[[1]]["y"]
  return(pars)
}