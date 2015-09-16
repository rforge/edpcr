conv <-
function(pars,from="m",to="c",det=FALSE) {
  if (from==to) if (det) return(list(new.pars=pars,dfac=1)) else return(pars)
  if (from=="c") cpars<-pars else if (from=="m") mpars<-pars
  
  if (from=="c"|(from=="p"&to=="m")) {
    if (from=="c") { # c to p
      pars<-cpars
      names(pars)[1]<-"mu"
      if ("disp"%in%names(cpars)) pars["mu"]<-compoisson::com.mean(cpars["lam"],cpars["disp"]) else pars["mu"]<-cpars["lam"]
      if (to=="p") new.pars<-pars
    }
    if (to=="m") { # p to m
      mpars<-pars
      mpars["mu"]<-log(pars["mu"])
      mpars["A"]<-log(pars["A"])#/log(1+pars["E"])
      mpars["E"]<- -log(1/(pars["E"])-1)
      if ("disp"%in%names(pars)) mpars["disp"]<-log(pars["disp"])
      if ("E1"%in%names(pars)) mpars["E1"]<- -log(1/(pars["E1"])-1)
      new.pars<-mpars
    }
  } else
  {
    if (from=="m") { # m to p
      pars<-mpars
      pars["mu"]<-exp(mpars["mu"])
      pars["E"]<- (1+exp(-mpars["E"]))^-1
      pars["A"]<- exp(mpars["A"])#(1+pars["E"])^mpars["A"]
      if ("disp"%in%names(pars)) pars["disp"]<-exp(mpars["disp"])
      if ("E1"%in%names(pars)) pars["E1"]<- (1+exp(-mpars["E1"]))^-1
      if (to=="p") new.pars<-pars
    }
    if (to=="c") { # p to c
      cpars<-pars
      names(cpars)[1]<-"lam"
      if ("disp"%in%names(cpars)) cpars["lam"]<-mu2lam(pars["mu"],pars["disp"]) else cpars["lam"]<-pars["mu"]
      new.pars<-cpars
    }
  }
  # determinant for transformation from pars to mpars
  if (det) {
    if (from!="m"|to!="c") dfac<-NA else {
      dfac<-1/pars["mu"]/(pars["E"]*(1-pars["E"]))/(pars["A"])#*log(1+pars["E"]))
      if ("disp"%in%names(pars)) dfac<-dfac/pars["disp"]
      if ("E1"%in%names(pars)) dfac<-dfac/(pars["E1"]*(1-pars["E1"]))
      dfac<-abs(prod(dfac))
    }
  }
  if (det) return(list(new.pars=new.pars,dfac=dfac)) else return(new.pars)
}