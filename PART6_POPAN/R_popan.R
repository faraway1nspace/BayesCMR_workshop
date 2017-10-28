# # PRACTICAL EXERCISE: POPAN in Bayesian HMM
# as part of the 2017 SMM Workshop "Bayesian Capture-Recapture" by Dr. R.W. Rankin and Mrs. K.E. Nicholson

# PART 1: Demonstration of a fully-time-varying POPAN-esque model with CMR data of Humpback dolphins (Sousa) from Hunt, T.N., Bejder, L., Allen, S.J., Rankin, R.W., Hanf, D., Parra, G.J., 2017. Demographic characteristics of Australian humpback dolphins reveal important habitat toward the south-western limit of their range. Endang Species Res 32, 71–88. doi:10.3354/esr00784. 
# PART 2: Model Selection Exercise: calculating the WAIC and Marginal Likelihoods
# THIS IS FOR EDUCATIONAL PURPOSES ONLY!!!

library(rjags)
# setwd() # set working directory to "BayesCMR_workshop/PART6_POPAN/"

# load the data
 # raw data in  "POPAN_6P37S_FINAL.txt"
# load and convert to matrix
y <- read.csv("hunt_popan_sousa.csv",header=FALSE)
T <- ncol(y) # number of primary periods (6-months intervals)
n <- nrow(y) # number of captured animals

# data augmentation (PXDA)
n.aug <- round(0.8*n) # number of ''pseudo-individuals"
m <- n.aug+n # total number of real and pseudo-individuals
y.aug <- rbind(y, matrix(0,nrow=n.aug,ncol=T))
nrow(y.aug)==m

# MODEL 1: fully-time-varying psi, phi and p
# WARNING: we must constrain p1=p2, and the last two phi's (confounding)
# priors: beta prior on p
pr.p = matrix(c(1,1,1,1,1,1,1,1,1,1),nrow=T-1,ncol=2,byrow=TRUE,dimnames=list(2:T, c("a","b"))) # rows are capture-periods, and columns are the beta shape parameters (notice there is no free p_1)
# priors: beta prior on phi
pr.phi = matrix(c(1,1,1,1,1,1,1,1),nrow=T-2,ncol=2,byrow=TRUE,dimnames=list(1:(T-2), c("a","b"))) # rows are capture-periods, and columns are the beta shape parameters (notice there is no free phi[T-1]
# priors: beta prior on psi
pr.psi = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1),nrow=T,ncol=2,byrow=TRUE,dimnames=list(1:T, c("a","b"))) # rows are capture-periods, and columns are the beta shape parameters

# jags model
jags.model.txt <- "model{
# priors on psi (time-varying)
for(t in 1:T){
  psi[t] ~ dbeta(pr.psi[t,1], pr.psi[t,2]) # prior on recruitment
} 
# priors on capture (NOTE: constraint ON p[1])
for(t in 1:(T-1)){ # (only T-1 stochastic nodes)
  # stochastic p for 2:T 
  p.rv[t] ~ dbeta(pr.p[t,1], pr.p[t,2]) # prior on p
  # set deterministic nodes to enforce constraint
  p[t+1] <- p.rv[t]
} 
p[1] <- p[2] # make p1 = p2
# priors on survival (NOTE: constraint on final phi)
for(t in 1:(T-2)){
  # stochastic phi for 1:(T-2)
  phi.rv[t] ~ dbeta(pr.phi[t,1], pr.phi[t,2]) # prior on p
  # deterministic nodes to enforce constraint on final phi
  phi[t+1] <- phi.rv[t]
} 
phi[1] <- 1 # fake node (there is no survival between t=0 and t=1)
phi[T] <- phi[T-1] # constrain final phi to equal previous phi
# TRANSITION AND EMISSION MATRICES
for(t in 1:T){
  # HMM TRANSITION MATRIX
  tr[1,1,t] <- 1-psi[t] # unborn to unborn
  tr[2,1,t] <- psi[t]   # unborn to alive
  tr[3,1,t] <- 0     # (illegal)
# FROM alive (col2) to...
  tr[1,2,t] <- 0     # (illegal)
  tr[2,2,t] <- phi[t]   # alive to alive
  tr[3,2,t] <- 1-phi[t] # alive to dead
# FROM dead (col3) to...
  tr[1,3,t] <- 0     # (illegal)
  tr[2,3,t] <- 0     # (illegal)
  tr[3,3,t] <- 1     # dead to dead
  # HMM EMISSION MATRIX
  # state 1: unborn 
  em[1,1,t]<-1 # (100% no capture)
  em[2,1,t]<-0  
# state 2: alive
  em[1,2,t]<-1-p[t]
  em[2,2,t]<-p[t]  
# state 3: dead 
  em[1,3,t]<-1 # (100% no capture)
  em[2,3,t]<-0
}
# HMM LATENT STATE PROCESS: at t=1
for(i in 1:M){
  # initialize first state in z=1
  z[i,1] ~ dcat(tr[,1,1]) 
  # CONDITIONAL LIKELIHOOD: at t=1
  y[i,1] ~ dcat(em[,z[i,1],1]) # conditional capture
  # loop though capture periods > 1
  for(t in 2:T){
     z[i,t] ~ dcat(tr[,z[i,t-1],t]) # z_t | z_t-1
     # CONDITIONAL LIKELIHOOD: for t>1
     y[i,t] ~ dcat(em[,z[i,t],t]) # conditional capture
   } # t
 }# M
# DERIVATIVES: Pop abundance and recruits
for(i in 1:M){
    # check when each individual was born/recruited
    recruit_i[i,1] <- equals(z[i,1],2)
    # check whether an individual was actually alive
    N_i[i,1] <- equals(z[i,1],2)
    for(t in 2:T){
       recruit_i[i,t] <- equals(z[i,t-1],1)*equals(z[i,t],2)
       N_i[i,t] <- equals(z[i,t],2)
    } # t
} # M
for(t in 1:T){
   # tally total population abundance per t
   N[t] <- sum(N_i[1:M,t]) # sum over all individuals
   # tally total recruits per t
   recruits[t] <- sum(recruit_i[1:M,t]) 
}
# Calculate POPAN pent (probability of entry)
cumprob[1] <- psi[1]
for(t in 2:T){  
  cumprob[t] <- psi[t]*prod(1-psi[1:(t-1)]) 
}
cumprob.norm <- sum(cumprob[1:T])
# POPAN Inclusion probabilities
for(t in 1:T){ 
   pent[t] <- cumprob[t]/cumprob.norm 
} #t
}
"
jags.file <- "JAGS_hunt_popan.JAG"
# save the JAGS model syntax to a local file
sink(file=jags.file) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection

# ASSEMBLE THE DATA: into a list
y.jags <- y.aug+1 # convert AUGMENTED data (0,1) into (1,2) for JAGS model
jags.data <- list(
    y = y.jags, # capture data
    T=T,
    M = nrow(y.jags),
    pr.psi=pr.psi, 
    pr.phi=pr.phi, 
    pr.p=pr.p 
)

# INITIALIZE ALL THE RANDOM VARIABLES (include z)
# a function to generate random values and see the MCMC chains
# MUST return a list with names of random variables
jags.inits.f <- function(){
    # initial CMR parameters
    phi.rv=rbeta(T-2,10,1)
    p.rv=rbeta(T-1,10,10)
    psi=rbeta(T,7,10)
    # initialize latent states (z)
    # NOTE: must be consistent with y.jags
    unborn_<-1; alive_<-2; dead_<-3
    z <- matrix(0,nrow(y.jags),T+2) # add two dummy columns (one at beginning and one at the end)
    for(i in 1:nrow(y.jags)){
        max.t <- if(any(y.jags[i,]==alive_)){ 2+max(which(y.jags[i,]==alive_))} else { T+2}
        min.t <- if(any(y.jags[i,]==alive_)){ min(which(y.jags[i,]==alive_))} else { 1}        
        z[i,] <- rep(alive_,ncol(z))
        z[i,1:min.t]<- unborn_
        z[i,max.t:ncol(z)]<- dead_
    }
    z <- z[,2:(ncol(z)-1)]
    return(list(phi.rv=phi.rv,
                psi=psi,
                p.rv=p.rv,
                z=z)
           )
    }

# COMPILE THE JAGS MODEL
m1 <- jags.model(file=jags.file, inits=jags.inits.f, data=jags.data,n.chains=3,n.adapt=1000)
# BURN-IN PHASE
update(m1, 5000)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m1,variable.names=c("p","phi","psi","N","recruits","cumprob.norm"), n.iter=80000,thin=200)
# collapse the 3 chains into a single matrix
post.matrix <- do.call("rbind",post) # (optional, for easily manipulation)
# check MCMC Convergence and mixing with gelman.diag, plot, acf, etc.

# CHECK THAT OUR DATA AUGMENTATION n.aug IS OKAY: need the cumumlative entry probability (called cumprob.norm' to be << 1
quantile(post.matrix[,"cumprob.norm"],c(0.5,0.9,0.95,0.975,0.99))
# 0.6902648 0.7607443 0.7795476 0.7973764 0.8289005 
# ANSWER?


######################################################################################
# PART2: Model selection with the WAIC
######################################################################################

# in this exercise, we are going to calculate the WAIC
library(rjags)
# setwd() # set working directory to "BayesCMR_workshop/PART6_POPAN/"

# load the data
 # raw data in  "POPAN_6P37S_FINAL.txt"
# load and convert to matrix
y <- read.csv("hunt_popan_sousa.csv",header=FALSE)
T <- ncol(y) # number of primary periods (6-months intervals)
n <- nrow(y) # number of captured animals

# data augmentation (PXDA)
n.aug <- round(0.8*n) # number of ''pseudo-individuals"
m <- n.aug+n # total number of real and pseudo-individuals
y.aug <- rbind(y, matrix(0,nrow=n.aug,ncol=T))
nrow(y.aug)==m

# MODEL 1: fully-time-varying psi, phi and p
# WARNING: we must constrain p1=p2, and the last two phi's (confounding)
# priors: beta prior on p
pr.p = matrix(c(1,1,1,1,1,1,1,1,1,1),nrow=T-1,ncol=2,byrow=TRUE,dimnames=list(2:T, c("a","b"))) # rows are capture-periods, and columns are the beta shape parameters (notice there is no free p_1)
# priors: beta prior on phi
pr.phi = matrix(c(1,1,1,1,1,1,1,1),nrow=T-2,ncol=2,byrow=TRUE,dimnames=list(1:(T-2), c("a","b"))) # rows are capture-periods, and columns are the beta shape parameters (notice there is no free phi[T-1]
# priors: beta prior on psi
pr.psi = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1),nrow=T,ncol=2,byrow=TRUE,dimnames=list(1:T, c("a","b"))) # rows are capture-periods, and columns are the beta shape parameters

# jags model
jags.model.txt <- "model{
# priors on psi (time-varying)
for(t in 1:T){
  psi[t] ~ dbeta(pr.psi[t,1], pr.psi[t,2]) # prior on recruitment
} 
# priors on capture (NOTE: constraint ON p[1])
for(t in 1:(T-1)){ # (only T-1 stochastic nodes)
  # stochastic p for 2:T 
  p.rv[t] ~ dbeta(pr.p[t,1], pr.p[t,2]) # prior on p
  # set deterministic nodes to enforce constraint
  p[t+1] <- p.rv[t]
} 
p[1] <- p[2] # make p1 = p2
# priors on survival (NOTE: constraint on final phi)
for(t in 1:(T-2)){
  # stochastic phi for 1:(T-2)
  phi.rv[t] ~ dbeta(pr.phi[t,1], pr.phi[t,2]) # prior on p
  # deterministic nodes to enforce constraint on final phi
  phi[t+1] <- phi.rv[t]
} 
phi[1] <- 1 # fake node (there is no survival between t=0 and t=1)
phi[T] <- phi[T-1] # constrain final phi to equal previous phi
# TRANSITION AND EMISSION MATRICES
for(t in 1:T){
  # HMM TRANSITION MATRIX
  tr[1,1,t] <- 1-psi[t] # unborn to unborn
  tr[2,1,t] <- psi[t]   # unborn to alive
  tr[3,1,t] <- 0     # (illegal)
# FROM alive (col2) to...
  tr[1,2,t] <- 0     # (illegal)
  tr[2,2,t] <- phi[t]   # alive to alive
  tr[3,2,t] <- 1-phi[t] # alive to dead
# FROM dead (col3) to...
  tr[1,3,t] <- 0     # (illegal)
  tr[2,3,t] <- 0     # (illegal)
  tr[3,3,t] <- 1     # dead to dead
  # HMM EMISSION MATRIX
  # state 1: unborn 
  em[1,1,t]<-1 # (100% no capture)
  em[2,1,t]<-0  
# state 2: alive
  em[1,2,t]<-1-p[t]
  em[2,2,t]<-p[t]  
# state 3: dead 
  em[1,3,t]<-1 # (100% no capture)
  em[2,3,t]<-0
}
# HMM LATENT STATE PROCESS: at t=1
for(i in 1:M){
  # initialize first state in z=1
  z[i,1] ~ dcat(tr[,1,1]) 
  # CONDITIONAL LIKELIHOOD: at t=1
  y[i,1] ~ dcat(em[,z[i,1],1]) # conditional capture
  # loop though capture periods > 1
  for(t in 2:T){
     z[i,t] ~ dcat(tr[,z[i,t-1],t]) # z_t | z_t-1
     # CONDITIONAL LIKELIHOOD: for t>1
     y[i,t] ~ dcat(em[,z[i,t],t]) # conditional capture
   } # t
 }# M
# DERIVATIVES: Pop abundance and recruits
for(i in 1:M){
    # check when each individual was born/recruited
    recruit_i[i,1] <- equals(z[i,1],2)
    # check whether an individual was actually alive
    N_i[i,1] <- equals(z[i,1],2)
    for(t in 2:T){
       recruit_i[i,t] <- equals(z[i,t-1],1)*equals(z[i,t],2)
       N_i[i,t] <- equals(z[i,t],2)
    } # t
} # M
for(t in 1:T){
   # tally total population abundance per t
   N[t] <- sum(N_i[1:M,t]) # sum over all individuals
   # tally total recruits per t
   recruits[t] <- sum(recruit_i[1:M,t]) 
}
# Calculate POPAN pent (probability of entry)
cumprob[1] <- psi[1]
for(t in 2:T){  
  cumprob[t] <- psi[t]*prod(1-psi[1:(t-1)]) 
}
cumprob.norm <- sum(cumprob[1:T])
# POPAN Inclusion probabilities
for(t in 1:T){ 
   pent[t] <- cumprob[t]/cumprob.norm 
} #t
# loglike: for WAIC estimation (conditional loglikelkhood, unfortunately)
for(i in 1:N.obs){
  for(t in 1:T){
     llit[i,t] <- logdensity.cat(y[i,t], em[1:E, z[i,t],t])
  } # t primary periods
  lli[i] <- sum(llit[i,1:T])
} # individuals N.obs
}
"
jags.file <- "JAGS_hunt_popan.JAG"
# save the JAGS model syntax to a local file
sink(file=jags.file) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection

# ASSEMBLE THE DATA: into a list
y.jags <- y.aug+1 # convert AUGMENTED data (0,1) into (1,2) for JAGS model
jags.data <- list(
    y = y.jags, # capture data
    T=T,
    M = nrow(y.jags),
    pr.psi=pr.psi, 
    pr.phi=pr.phi, 
    pr.p=pr.p,
    N.obs = n, # number of capture individuals
    E = 2 # number of capture events
)

# INITIALIZE ALL THE RANDOM VARIABLES (include z)
# a function to generate random values and see the MCMC chains
# MUST return a list with names of random variables
jags.inits.f <- function(){
    # initial CMR parameters
    phi.rv=rbeta(T-2,10,1)
    p.rv=rbeta(T-1,10,10)
    psi=rbeta(T,7,10)
    # initialize latent states (z)
    # NOTE: must be consistent with y.jags
    unborn_<-1; alive_<-2; dead_<-3
    z <- matrix(0,nrow(y.jags),T+2) # add two dummy columns (one at beginning and one at the end)
    for(i in 1:nrow(y.jags)){
        max.t <- if(any(y.jags[i,]==alive_)){ 2+max(which(y.jags[i,]==alive_))} else { T+2}
        min.t <- if(any(y.jags[i,]==alive_)){ min(which(y.jags[i,]==alive_))} else { 1}        
        z[i,] <- rep(alive_,ncol(z))
        z[i,1:min.t]<- unborn_
        z[i,max.t:ncol(z)]<- dead_
    }
    z <- z[,2:(ncol(z)-1)]
    return(list(phi.rv=phi.rv,
                psi=psi,
                p.rv=p.rv,
                z=z)
           )
    }

# COMPILE THE JAGS MODEL
m1 <- jags.model(file=jags.file, inits=jags.inits.f, data=jags.data,n.chains=3,n.adapt=1000)
# BURN-IN PHASE
update(m1, 5000)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m1,variable.names=c("p","phi","psi","N","recruits","cumprob.norm","lli"), n.iter=80000,thin=200)
# collapse the 3 chains into a single matrix
post.matrix <- do.call("rbind",post) # (optional, for easily manipulation)
# check MCMC Convergence and mixing with gelman.diag, plot, acf, etc.


# WAIC approximation 1 function: (WAIC1; from Gelman et al 2014)
waic1.jags <- function(post,n.obs,loglike.grep.txt = "lli"){
    logsumexp <- function(x){ max.x <- max(x); max.x - log(length(x)) +log(sum(exp(x-max.x)))} # this is equivalient to log( 1/S * sum_s[exp{ x }]))    
    # get the columns indices that have 'lli' in their name (for loglike-inidividual)
    ll.col <-grep(loglike.grep.txt,colnames(post[[1]]))
    # get the mcmc values for the log values
    if(length(post)>1){ lli.samp <- as.matrix(do.call("rbind",post)[,ll.col]) } else { lli.samp <- as.matrix(post[[1]][,ll.col])}
    # get the complexity penality (WAIC1; from Gelman et al 2014)
    logElike <- apply(lli.samp,2,logsumexp)
    Eloglike <- apply(lli.samp,2,mean)
    p_waic1 <- 2*sum(logElike-Eloglike)
    # get the lppd log pointwaise predictive density
    # use the log-sum-exp trick to handle underflow issues with exp(loglike) ~= 0
    logsumf.i <- apply(lli.samp,2,logsumexp)
    lppd <- sum(logsumf.i)
    waic1 <- -2*(lppd-p_waic1)
    return(waic1)
}

# WAIC approximation 2 function: (WAIC2; from Gelman et al 2014)
waic2.jags <- function(post,n.obs,loglike.grep.txt = "lli"){
    logsumexp <- function(x){ max.x <- max(x); max.x - log(length(x)) +log(sum(exp(x-max.x)))} # this is equivalient to log( 1/S * sum_s[exp{ x }]))    
    # get the columns indices that have 'lli' in their name (for loglike-inidividual)
    ll.col <-grep(loglike.grep.txt,colnames(post[[1]]))
    # get the mcmc values for the log values
    if(length(post)>1){ lli.samp <- as.matrix(do.call("rbind",post)[,ll.col]) } else { lli.samp <- as.matrix(post[[1]][,ll.col])}
    # get the complexity penality (WAIC2; from Gelman et al 2014)
    V.loglike.i<-apply(lli.samp,2,function(loglike){ 1/(length(loglike)-1) * sum((loglike-mean(loglike))^2) })
    p_waic2 = sum(V.loglike.i)
    # get the lppd log pointwaise predictive density
    # use the log-sum-exp trick to handle underflow issues with exp(loglike) ~= 0
    logsumf.i <- apply(lli.samp,2,logsumexp)
    lppd <- sum(logsumf.i)
    waic2 <- -2*(lppd-p_waic2)
    return(waic2)
}


# approximate the WAIC using the log.density of capture-histories outputs from JAGS
waic1 <- waic1.jags(post)
waic2 <- waic2.jags(post)
