p2mu <-
function(p0,disp,n=100) {
    if (length(p0)>1&length(disp)>1) stop("Need either p0 or disp to be of length 1")
    if (length(p0)>1) sapply(p0,function(p) p2mu(p,disp,n)) else
        if (length(disp)>1) sapply(disp,function(d) p2mu(p0,d,n)) else
        {
            lps<-(0:n)*log(p2lam(p0,disp))-disp*lfactorial(0:n)
            return(sum((1:n)*exp(lps[-1]))/sum(exp(lps)))
        }
}
