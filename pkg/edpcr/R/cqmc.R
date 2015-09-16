cqmc<-
function(data,mc.rep=10^4,h=NULL,n0=NULL,n1=NULL,nt0=0,nt1=0,n.lo=0,n.hi=0,
 pars0=NULL,report=1,probreport=.1,
 extra=c("trendx","trendy","disp","E1"),c0=6,maxn0=7,mod.method="Nelder-Mead",
 mod.rep=2000,burnin=0,nreport=10,cq.xlim=NA,
 tune=1,E1.init=.8,E.init=.9,d.init=1,prior=TRUE,
 mu.fun=function(x) dgamma(x,1.5,1.5), A.fun=function(x) x^-1,E.fun=function(x) dbeta(x,60,5),
 E1.fun=E.fun,trendx.fun=function(x) 1,trendy.fun=function(x) 1,disp.fun=function(x) dgamma(x,10,10)) {
  

  if (!prior) mu.fun<-A.fun<-E.fun<-E1.fun<-trendx.fun<-trendy.fun<-disp.fun<-function(x) 1
  mc<-(mc.rep>0) 
  mod<-(mod.rep>0) 
  if (class(data)=="cqdat") {
    h<-data$h; if (is.na(h)) stop("need valid h")
    nt0<-data$nt0;nt1<-data$nt1;n0<-data$n0;n1<-data$n1
    n.hi<-data$n.hi;n.lo<-data$n.lo;nx<-data$nx;ny<-data$ny
    cqdata<-data$data 
    if (report>0) cat("counts nt0",nt0,"nt1",nt1,"n0",n0,"n1",n1,"n.hi",n.hi,"n.lo",n.lo,"\n",sep=" ")
    
  } else stop("need data with class cqdat")
  counts<-c(n0=n0,n1=n1,n.lo=n.lo,n.hi=n.hi,nt0=nt0,nt1=nt1)
  if (is.null(pars0)) pars0<-get.pars0(cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,extra=extra,maxn0=7,
                                       d.init=d.init,E.init=E.init[1],E1.init=E1.init[1])
  mpars0<-conv(pars0,"p","m")
  logval0<-unname(llcq(mpars0,cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,nx=nx,ny=ny,c0=c0,maxn0=maxn0,report=0,
             full=TRUE,extra=extra,prior=prior,
             mu.fun=mu.fun,E.fun=E.fun,E1.fun=E1.fun,A.fun=A.fun,trendx.fun=trendx.fun,trendy.fun=trendy.fun,disp.fun=disp.fun))
   if (mod|mc) {
    parscale<-mpars0*0
    parscale[1:3]<-c(.1,.1,.1)
    if ("trendx"%in%extra) parscale["trendx"]<-.001
    if ("trendy"%in%extra) parscale["trendy"]<-.001
    if ("E1"%in%extra) parscale["E1"]<-.1
    if ("disp"%in%extra) parscale["disp"]<-.01
  }
  
  if (mod) {
    cat("\nMode stage")
    modrep<-max(length(E.init),length(E1.init))
    E.init<-rep(E.init,modrep/length(E.init))
    E1.init<-rep(E1.init,modrep/length(E1.init))
    logval.mod<- -Inf
    mpars0<-mpars.mod0<-conv(pars0,"p","m")

    if (modrep>1) mod.res<-c() else mod.res<-NA
    for (repp in 1:modrep) {
      if (repp>1) if(E.init[repp]!=E.init[repp-1]|E1.init[repp]!=E1.init[repp-1]) {
        mpars.mod0<-conv(pars0<-get.pars0(cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,extra=extra,maxn0=7,
                       d.init=d.init,E.init=E.init[repp],E1.init=E1.init[repp]),"p","m")
      }
      if (report>0) cat("\n",E.init[repp],E1.init[repp])
      logval.0<-unname(llcq(mpars.mod0,cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,nx=nx,ny=ny,c0=c0,maxn0=maxn0,report=0,
                           full=TRUE,extra=extra,prior=prior,
                           mu.fun=mu.fun,E.fun=E.fun,E1.fun=E1.fun,A.fun=A.fun,trendx.fun=trendx.fun,trendy.fun=trendy.fun,disp.fun=disp.fun))    
      op<-optim(mpars.mod0,llcq,control=list(parscale=parscale,fnscale=-1,maxit=mod.rep),
       data=cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,nx=nx,ny=ny,c0=c0,method=mod.method,probreport=probreport,
       report=report,full=TRUE,extra=extra,prior=prior,
       mu.fun=mu.fun,E.fun=E.fun,E1.fun=E1.fun,A.fun=A.fun,trendx.fun=trendx.fun,trendy.fun=trendy.fun,disp.fun=disp.fun)
      if (report>1) {cat("\n");print(c(op$par,logval0=logval.0,logval=op$value))}
      if (modrep>1) mod.res<-rbind(mod.res,c(signif(conv(op$par,"m","p"),4),logval0=logval.0,logval=op$value))
      if (repp==1|(repp>1&(op$value>logval.mod))) {
        op.mod<-op
        mpars0<-mpars.mod0
        logval0<-logval.0
        logval.mod<- op$value
        mpars.mod<-op$par
        names(mpars.mod)<-names(mpars.mod0)
      }
    }
    cat("\n")
    if (modrep>1) print(mod.res<-cbind(E.init=E.init,E1.init=E1.init,mod.res)) else
      print(c(E.init=E.init,E1.init=E1.init,mpars.mod,logval.0=logval0,logval.mod=logval.mod))
    pars.mod<-conv(mpars.mod,"m","p")
  }

  if (mc) {
    cat("\nMCMC stage")
    print(system.time({
    verbose<-round((1+burnin)*mc.rep/nreport)
    if (mod) mpars0.mc<-mpars.mod else mpars0.mc<-mpars0
    if (!mod) V<-NULL else { # Following code to calculate V modified from MCMCmetrop1R to avoid need to rerun optim
      if (op.mod$convergence != 0) {
        warning("Mode and Hessian were not found with call to optim().\nSampling proceeded anyway. \n")
      }
      cat("\nEstimate variance matrix\n")
      hessian<-optimHess(mpars.mod,llcq,control=list(parscale=-parscale),
                         data=cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,nx=nx,ny=ny,c0=c0,
                         report=report,full=TRUE,extra=extra,prior=prior,
                         mu.fun=mu.fun,E.fun=E.fun,E1.fun=E1.fun,A.fun=A.fun,trendx.fun=trendx.fun,trendy.fun=trendy.fun,disp.fun=disp.fun)
      CC <- NULL
      try(CC <- chol(-1 * hessian), silent = TRUE)
      hess.new <- hessian
      hess.flag <- 0
      
      if (max(diag(hessian) == 0)) {
        for (i in 1:nrow(hess.new)) {
          if (hess.new[i, i] == 0) {
            hess.new[i, i] <- -1e-06
          }
        }
      }
      while (is.null(CC)) {
        hess.flag <- 1
        hess.new <- hess.new - diag(diag(0.01 * abs(hessian)))
        try(CC <- chol(-1 * hess.new), silent = TRUE)
      }
      
      if (hess.flag == 1) {
        warning("Hessian from call to optim() not negative definite.\nSampling proceeded after enforcing negative definiteness. \n")
      }
      V <- solve(-1 * hess.new)
      if (report>1) print(V)
    }
    mpars.chain<-MCMCpack::MCMCmetrop1R(llcq,mpars0.mc,burnin=burnin*mc.rep,mcmc=mc.rep,verbose=verbose,V=V,
     data=cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,nx=nx,ny=ny,c0=c0,report=report,prior=prior,extra=extra,
     mu.fun=mu.fun,E.fun=E.fun,E1.fun=E1.fun,A.fun=A.fun,trendx.fun=trendx.fun,trendy.fun=trendy.fun,disp.fun=disp.fun,
     force.samp=TRUE,probreport=probreport,logfun=TRUE,tune=tune,full=TRUE,optim.control=list(optim.control=list(fnscale=-parscale)))}))
    
    colnames(mpars.chain)<-names(mpars0.mc)
    pars.chain<-conv.mcmc(mpars.chain,"m","p")

    print(pars.sum<-summary(pars.chain)[[1]])
    pars.mc<-pars.sum[,1]; print(pars.mc)
    mpars.mc<-conv(pars.mc,"p","m")
    logval.mc<-unname(llcq(mpars.mc,cqdata,h=h,nt0=nt0,nt1=nt1,n0=n0,n1=n1,n.hi=n.hi,n.lo=n.lo,nx=nx,ny=ny,c0=c0,maxn0=maxn0,report=report,
                full=TRUE,extra=extra,prior=prior,
                mu.fun=mu.fun,E.fun=E.fun,E1.fun=E1.fun,A.fun=A.fun,trendx.fun=trendx.fun,trendy.fun=trendy.fun,disp.fun=disp.fun))
  if (report>0) print(mpars.mc)
  }
  if (report>0) {
    cat("\n")
    print(c(logval0=logval0))
    if (mod) print(c(logval.mod=logval.mod))
    if (mc) print(c(logval.mc=logval.mc))
  }

  prior.funs<-list(mu.fun=mu.fun,A.fun=A.fun,E.fun=E.fun,E1.fun=E1.fun,trendx.fun=trendx.fun,trendy.fun=trendy.fun,disp.fun=disp.fun)

  output<-list(data=cqdata,counts=counts,pars0=pars0,logval0=logval0,h=h,nx=nx,ny=ny)
  if (mod) output<-c(output,list(pars.mod=pars.mod,logval.mod=logval.mod,mod.res=mod.res))
  if (mc)  output<-c(output,list(pars.mc= pars.mc, logval.mc= logval.mc,pars.sum=pars.sum,pars.chain=pars.chain,
                                 mc.acc=1-coda::rejectionRate(pars.chain)[1],prior.funs=prior.funs))
  class(output)<-"cqres"
  return(output)
}