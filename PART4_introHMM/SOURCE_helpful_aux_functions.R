# these are just some session specific helpful R functions (plotting, etc).
# these SHOULD NOT be used outside the tutorial context, they probably won't work
# I just made these to not distract the Workshop participants from frivolous (but pretty) functions

#posterior = post.matrix[,grep("N\\[",colnames(post.matrix))]; xaxis.lab = 1:11; col95="grey60"; col68="grey40";colmode="black";pch.mode=19;cex.mode=1.1;pch.mean=4;cex.mean=1.1;colmean="white";par.args=list(mar=c(3,3,0.5,0.5),mgp=c(1.7,0.7,0))
plot.posterior.timeseries<-function(posterior, ylim=range(posterior), main="", xlab, ylab,x.lim,xaxis.lab,col95="grey60",col68="grey40",colmode="black",pch.mode=19,cex.mode=1,pch.mean=4,cex.mean=1.1,colmean="blue",par.args=list(mar=c(3,3,3*(main!="")+0.5*(main==""),0.5),mgp=c(1.7,0.7,0)),legend.args=list(x="topleft",col=c(colmean,colmode,col68,col95), legend=c("mean","MAP","mean+/-SE","95%CI"), fill=c(NA,NA,col68,col95),bty="n",pch=c(pch.mean,pch.mode,NA,NA),cex=1.2)){
    if(!is.null(par.args)){ do.call("par",args=par.args)}
    if(sum(abs(posterior[,1]-round(posterior[,1])))<10^-7){
        mode.f <- function(x) as.integer(names(which.max(table(x))))
    } else {
        mode.f <- function(x){ d<-density(x,adjust=1.2); return(d$x[which.max(d$y)])}
    }
    s <- apply(posterior,2,function(x){ c(mean = mean(x), median = median(x), mode=mode.f(x), lcl95=quantile(x,0.025,names=FALSE), lcl68=quantile(x,pnorm(-1),names=FALSE), ucl68=quantile(x,pnorm(1),names=FALSE), ucl95=quantile(x,1-0.025,names=FALSE))})
    # blank canvas
    xrange <- ncol(posterior)
    plot(c(1,xrange),ylim,xlab=xlab,ylab=ylab,ylim=ylim,type="n", xaxt="n",main=main)
    axis(1,at = 1:xrange, labels=xaxis.lab)
                                        # make the 95% polygon
    polygon(x=c(1:xrange, xrange:1), c(s["lcl95",],rev(s["ucl95",])),border=NA,col=col95)
    # make the mean-/+SD polygon
    polygon(x=c(1:xrange, xrange:1), c(s["lcl68",],rev(s["ucl68",])),border=NA,col=col68)    
    # add the mode
    points(1:xrange, s["mode",],pch=pch.mode,cex=cex.mode,col=colmode)
    # add the mean
    lines(1:xrange, s["mean",],pch=pch.mean,cex=cex.mean,col=colmean,typ="b")
    if(!is.null(legend.args)){ do.call("legend",args=legend.args) }
}

#colmean="red";pch.mean="x";cex.mean=1.2;x.offset=0.1;col95="red"
plot.mle.points <- function(mle.mat,add=TRUE,colmean="red",pch.mean=4,cex.mean=1.5,x.offset=0.1,col95="red",legend.args=list(x="top",col="red",legend="MLE",lty=1,pch=4,bty="n")){
    xrange <- nrow(mle.mat)
    points((1:xrange)+x.offset,mle.mat[,1], cex=cex.mean,col=colmean,pch=pch.mean)
    if(any(colnames(mle.mat)=="lcl") & any(colnames(mle.mat)=="ucl")){
        for(i in 1:nrow(mle.mat)){
            lines(rep(i,2)+x.offset, mle.mat[i,c("lcl","ucl")],col=col95,lwd=2)
            lines(i+c(0,2*x.offset),mle.mat[i,c("lcl","lcl")],col=col95,lwd=2) # tick
            lines(i+c(0,2*x.offset),mle.mat[i,c("ucl","ucl")],col=col95,lwd=2) # tick 
        }
    }
    if(!is.null(legend.args)){ do.call("legend",args=legend.args) }    
}
