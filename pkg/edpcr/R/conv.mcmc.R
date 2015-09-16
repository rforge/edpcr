conv.mcmc <-
function(chain,from="m",to="p") {
    conv.env<-environment()
    last.orig<-c()
    last.new<-c()
    coda::as.mcmc(t(apply(chain,1,function(x) {
        last.orig<-get("last.orig",envir=conv.env)
        names(x)<-colnames(chain)
        if (identical(x,last.orig)) return(get("last.new",envir=conv.env)) else {
            tmp<-conv(x,from,to,det=FALSE)
            names(tmp)[1]<-"mu"
            assign("last.orig",x,envir=conv.env)
            assign("last.new",tmp,envir=conv.env)
            tmp
        }
    })))
}
