detrend <-
function(data=NULL,detr="xy",nx=NULL,ny=NULL,inc.bs=FALSE) {
    if (class(data)=="cqdat") data<-data[[1]]
    if (is.null(nx)) nx<-max(data$x)
    if (is.null(ny)) ny<-max(data$y)
    if (identical(detr,"x")) {
        trend<-(lm<-robust::lmRob(cq~x,data=data))[[1]]["x"]
        b.x<-trend;b.y<-0
        data$cq<-data$cq-trend*(data$x-nx/2)
    } else 
        if (identical(detr,"y")) {
            trend<-robust::lmRob(cq~y,data=data)[[1]]["y"]
            data$cq<-data$cq-trend*(data$y-ny/2)
            b.x<-0;b.y<-trend
        } else 
            if (identical(detr,"xy")) {
                trend<-robust::lmRob(cq~x+y,data=data)[[1]]
                data$cq<-data$cq-trend["x"]*(data$x-nx/2)-trend["y"]*(data$y-ny/2)
                b.x<-trend["x"];b.y<-trend["y"]
            }
    if (!inc.bs) return(data) else return(list(data=data,b.x=b.x,b.y=b.y))
}
