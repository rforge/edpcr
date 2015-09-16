p2lam <-
function(p0,disp,n=100) {
    if (disp<0) return(NA)
    if (disp==1) return(-log(p0))
    if (disp<1&-log(p0)>1) {
        l1<-1
        while(sum(exp((0:n)*log(l1)-disp*lfactorial(0:n)))-p0^-1<0) l1<-l1+p0 
    } else if (disp>1) l1<-1/p0-.99 else l1<- -log(p0)+10^-3
    lam<-uniroot(function(lam,p0,disp) {
        if (lam==0) 1 - p0^-1 else 
            sum(exp((0:n)*log(lam)-disp*lfactorial(0:n)))-p0^-1
    },
    interval=c(0,l1),p0=p0,disp=disp)$root
    return(lam)
}
