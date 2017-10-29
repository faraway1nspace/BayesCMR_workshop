# For making initialization of latent states
    lexiscope <-function(y.aug,T,te1=3,te2=4,re=5,de=6,out_=3,n.sex=2){
        force(y.aug);force(de);force(re);force(te1);force(te2);force(T);force(n.sex)
        retf <- function(){
            ypp <- z<- apply(y.aug,c(1,3), function(x) min(x,na.rm=TRUE))
            first <- apply(ypp,1,function(x) min(which(x!=out_)))
            first[!is.finite(first)]<-sample(2:(T-1),size=sum(!is.finite(first)),replace=TRUE)
            last <- apply(ypp,1,function(x) max(which(x!=out_)))
            last[!is.finite(last)]<-sapply(first[!is.finite(last)],function(ft) max(ft,sample(ft:T,1)))
            for(j in 1:nrow(ypp)){
                if(any(ypp[j,]!=out_)){
                    mostcommon.z <- min(which(table(ypp[j,])[c("1","2")]==max(table(ypp[j,])[c("1","2")],na.rm=TRUE)))
                } else {
                    mostcommon.z <- sample(1:2,1)
                    ypp[j,first[j]]<-mostcommon.z
                    ypp[j,last[j]]<-mostcommon.z
                }
                leastcommon.z <- c(2,1)[mostcommon.z]
                if(first[j]!=1){
                    (fakepre<-sort(sample(c(-1,mostcommon.z,leastcommon.z),size=length(1:first[j])-1,replace=TRUE,prob=c(0.25,0.5,0.25))))
                    ypp[j,1:(first[j]-1)] <- fakepre*(fakepre>0) +  rep(re,length(fakepre))*(fakepre<=0)
                }
                if(last[j]!=T){
                    ypp[j,(last[j]+1):T] <- sort(sample(c(mostcommon.z,leastcommon.z,de),size=length(last[j]:T)-1,replace=TRUE,prob=c(0.4,0.2,0.4)))
                }
                if(first[j]!=last[j]){
                    zrow <- sample(c(mostcommon.z,leastcommon.z)+c(0,0,2,0), size=length((first[j]+1):(last[j]-1)),replace=TRUE,prob=c(0.4,0.05,0.5,0.05))
                    ypp[j,(first[j]+1):(last[j]-1)] <- zrow*(ypp[j,(first[j]+1):(last[j]-1)]==out_) + ypp[j,(first[j]+1):(last[j]-1)]*((ypp[j,(first[j]+1):(last[j]-1)]!=out_))
                }
            }
            return(list(z=ypp,llambda_mu=runif(1,0.45,0.55),lgte_mu=runif(1,0.05,0.2), lgaa_mu=runif(1,0.8,0.95), lgbb_mu=runif(1,0.8,0.95), lphi_mu=runif(1,0.94,0.98),lp_mu=runif(1,0.01,0.2),psi = matrix(c(0.5,runif(T-1,0.01,0.05)),T,n.sex)))}
        return(retf)
    }
    # INITIALIZE THE LATENT STATES

# WAIC calculation (from Gelman et al 2014)
waic.jags <- function(post,n.obs,loglike.grep.txt = "lli"){
    # get the columns indices that have 'lli' in their name (for loglike-inidividual)
    ll.col <-grep(loglike.grep.txt,colnames(post[[1]]))
    # get the mcmc values for the log values
    if(length(post)>1){ lli.samp <- as.matrix(do.call("rbind",post)[,ll.col]) } else { lli.samp <- as.matrix(post[[1]][,ll.col])}
    # get the complexity penality (WAIC2; from Gelman et al 2014)
    V.loglike.i<-apply(lli.samp,2,function(loglike){ 1/(length(loglike)-1) * sum((loglike-mean(loglike))^2) })
    p_waic2 = sum(V.loglike.i)
    # get the lppd log pointwaise predictive density
    # use the log-sum-exp trick to handle underflow issues with exp(loglike) ~= 0
    logsumexp <- function(x){ max.x <- max(x); max.x - log(length(x)) +log(sum(exp(x-max.x)))} # this is equivalient to log( 1/S * sum_s[exp{ x }]))
    logsumf.i <- apply(lli.samp,2,logsumexp)
    lppd <- sum(logsumf.i)
    waic2 <- -2*(lppd-p_waic2)
    return(waic2)
}
#loglike <- c(-2,-1.2,-2.3,-0.2,-0.01) # test the logsumexp trick
#log((1/length(loglike)) *sum(exp(loglike)))
#logsumexp(loglike)

    
    
