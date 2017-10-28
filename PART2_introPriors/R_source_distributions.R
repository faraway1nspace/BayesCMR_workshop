library(MCMCpack)
library(boot)
library(truncnorm)

### other distributions
# Beta: see rbeta from R base
# Gamma: see the rgamma(n, shape, rate=1, scale=1/rate) (R base)
# Inverse-Gamma: see the rinvgamma(n, shape, rate = 1) from MCMCpack
# Dirichlet: see rdirichlet(n, alpha) from MCMCpack

# random samples from Normal, specified with tau = sd^-2 (JAGS/Bayesian parameterization)
rnorm.jags <- function(n=1,mean=0,tau=1){
    if(tau<=0 | round(n)<1 | n<1){ stop("tau must be >0; n must be a positive integer")}
    rnorm(n,mean, sd = sqrt(1/tau))
}

# random samples from Normal, with an inverse-logit transformation (maps onto the [0,1] scale), specified with tau = sd^-2 (JAGS/Bayesian parameterization) 
rlnorm.jags <- function(n=1,mean=0,tau=1){
    if(tau<=0 | round(n)<1 ){ stop("tau must be >0; n must be a positive integer")}
    inv.logit(rnorm(n,mean, sd = sqrt(1/tau)))
}

# random samples from Normal, with an PROBIT transformation (onto the [0,1] scale), specified with tau = sd^-2 (JAGS/Bayesian parameterization)logit Normal
rpnorm.jags <- function(n=1,mean=0,tau=1){
    if(tau<=0 | round(n)<1 ){ stop("tau must be >0; n must be a positive integer")}    
    pnorm(rnorm(n,mean, sd = sqrt(1/tau)))
} 

# random samples from a truncated normal, truncated between a and b. Specified with tau = sd^-2 (JAGS/Bayesian parameterization)
rtruncnorm.jags <- function(n=1,mean=0,tau=1,a=0,b=Inf){
    if(tau<=0 | round(n)<1  | a>=b ){ stop("n must be a positive integer; tau must be >0; a must be < b")}
    rtruncnorm(n, a=a,b=b,mean = mean, sd = 1/sqrt(tau)) 
}

# random samples from the Half Normal (a Normal truncated at 0 below). Specified with tau = sd^-2 (JAGS/Bayesian parameterization)
rhalfnorm.jags <- function(n=1,mean=0,tau=1){
    if(tau<=0 | round(n)<1 ){ stop("tau must be >0; n must be a positive integer")}    
    rtruncnorm(n, a=0,b=Inf,mean = mean, sd = 1/sqrt(tau)) 
}

# random samples from scaled-half-student-t ( a student-t truncated at 0 below). Specified with tau = sd^-2 (JAGS/Bayesian parameterization). Mean is 0. tau = 1/variance. nu is the degrees-of-freedom
rhalft.jags <- function(n=1,tau=1,nu=4){
   if(tau<=0 | round(n)<1 | nu <=0){ stop("n must be a positive integer; tau must be >0; nu must be >0")}        
    abs(1/sqrt(tau)*rt(n=n,df=nu))
}

summaryf <- function(x,probs = c(0.025,0.1587,0.5,0.841,0.975),plot=TRUE,breaks=20,col="grey"){
    mn = mean(x)
    s = sd(x)
    q = quantile(x,probs=probs)
    max = max(x)
    if(plot){ hist(x,breaks=breaks,col=col) }
    return(c(mean=mn,sd=s, q, max=max))}
