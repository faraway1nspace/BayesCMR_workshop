# # PRACTICAL EXERCISE: PCRD in Bayesian HMM
# as part of the 2017 SMM Workshop "Bayesian Capture-Recapture" by Dr. R.W. Rankin and Mrs. K.E. Nicholson.
# Data from: Nicholson, K., Bejder, L., Allen, S.J., Krützen, M., Pollock, K.H., 2012. Abundance, survival and temporary emigration of bottlenose dolphins (Tursiops sp.) off Useless Loop in the western gulf of Shark Bay, Western Australia. Mar. Freshwater Res. 63, 1059–1068. doi:10.1071/MF12210
# Models from: Rankin, R.W., Nicholson, K.E., Allen, S.J., Krützen, M., Bejder, L., Pollock, K.H., 2016. A full-capture Hierarchical Bayesian model of Pollock’s Closed Robust Design and application to dolphins. Front Mar Sci 3. doi:10.3389/fmars.2016.00025

####################################################################################################
# PART 1: Demonstration of an individual-heterogeneity PCRD model
# DON'T RUN: takes too long:set the working direct to 
library(boot)
library(rjags)
setwd(??)

# load some helpful files
source("R_PCRD_JAGS_SOURCE.R")

MARK.file.name <- "Nicholson_pcrd_ch.inp"
T2 <- c(5, 5, 10, 5,3) # number of secondary periods per primary period
T <- length(T2) # number of primary periods
capture.histories = import.mark.inp(MARK.file.name,T2,header=FALSE) # import MARK inp
Y.tt <- capture.histories[["Y.tt"]] # get the 3D array N x T2 x T
n = nrow(Y.tt) # number of (observed) individuals 

# jags model text
jags.txt <-"model{
  # HYPERPRIORS
  sigma.eps ~ dt(0,pr.sigma.eps[1],pr.sigma.eps[2]) T(0,) # controls dispersion of individual heterogeneity in capture
  sigma.p.t ~ dt(0,pr.sigma.p.t[1],pr.sigma.p.t[2]) T(0,) # controls dispersion of primary periods capture
  sigma.p.tt ~ dt(0,pr.sigma.p.tt[1],pr.sigma.p.tt[2]) T(0,) # controls dispersion of 2ndary periods capture
  # hierarchical capture process
  lp.mu ~ dnorm(pr.lp.mu[1],pr.lp.mu[2]) # mean capture (logit)
  for(t in 1:T){ # loop through primary periods
    lp.t[t] ~ dnorm(0, pow(sigma.p.t,-2))
    for(s in 1:T2[t]){ # loop through secondary periods
      lp.tt[s,t] ~ dnorm(0, pow(sigma.p.tt,-2))
    }
  }
  # PRIORS: on process parameters
  lphi.mu ~ dnorm(pr.lphi.mu[1],pr.lphi.mu[2]) # survival mean (logit scale)
  lg1.mu ~ dnorm(pr.g1.mu[1],pr.g1.mu[2]) # gamma-prime, mean (logit scale)
  lg2.mu ~ dnorm(pr.g2.mu[1],pr.g2.mu[2]) # gamma-prime-prime, mean (logit scale)
  # time-varying psi
  for(t in 1:T){
    psi[t] ~ dunif(pr.psi[t,1], pr.psi[t,2])
  }
  # shunt time-constant parameters to time-varying
  for(t in 1:(T-1)){ #loop through primary periods
    g1[t] <- ilogit(lg1.mu) #prob migrate out|out (probability)
    g2[t] <- ilogit(lg2.mu) #prob migrate in|out (probability)
    phi[t] <- ilogit(lphi.mu) # survival (probability)
  }
  # HIDDEN MARKOV MODEL: TRANSITION MATRIX
  # first state transition (pure nusance; strictly from outside-pop to part of marked-pop)
  tr0[1] <- (1-psi[1]) #remains not-yet-in-pop
  tr0[2] <- psi[1] # recruits
  tr0[3] <- 0 # cannot go to TE
  tr0[4] <- 0 # cannot go to dead
  # transition matrix from (2:T)
  for(t in 1:(T-1)){
     # from state 1 (un-recruited)
     tr[1,1,t] <- 1-psi[t+1] # does not recruit
     tr[2,1,t] <- psi[t+1] # recruits
     tr[3,1,t] <- 0 # illegal
     tr[4,1,t] <- 0 # illegal
     # from state 2 (onsite)
     tr[1,2,t] <- 0
     tr[2,2,t] <- phi[t]*(1-g2[t]) # stays onsite
     tr[3,2,t] <- phi[t]*g2[t] # leaves to TE
     tr[4,2,t] <- 1-phi[t] # dies
     # from state 3 (TE)
     tr[1,3,t] <- 0
     tr[2,3,t] <- phi[t]*(1-g1[t]) # returns onside
     tr[3,3,t] <- phi[t]*g1[t] # stays in TE
     tr[4,3,t] <- 1-phi[t] # dies
     # death
     tr[1,4,t] <- 0
     tr[2,4,t] <- 0
     tr[3,4,t] <- 0
     tr[4,4,t] <- 1 # can only stay dead
  } #t
  #loop through M individuals
  for(i in 1:M){
    # BUILD EMISSION MATRIX (for educational purposes, there are better ways of going this)
    # Notice that because we have individual heterogeneity, we have individual-specific emission matrices
    lp.i[i] ~ dnorm(0,pow(sigma.eps,-2))
    for(t in 1:T){ # loop through primary periods
      for(s in 1:T2[t]){ # loop secondary periods
        # unrecruited: state 1
        em[1,1,t,s,i] <- 1 # no capture
        em[2,1,t,s,i] <- 0 # capture illegal
        # onsite: state 2
        em[2,2,t,s,i] <- 1/(1+exp(-1*(lp.mu+lp.t[t]+lp.tt[s,t]+lp.i[i]))) # capture
        em[1,2,t,s,i] <- 1-em[2,2,t,s,i] # no capture
        # TE: state 3
        em[1,3,t,s,i] <- 1 # 100% no capture
        em[2,3,t,s,i] <- 0 # 
        # dead: state 4
        em[1,4,t,s,i] <- 1 # no capture
        em[2,4,t,s,i] <- 0 # 
      } # s 
    } # t
    # LATENT STATE FOR INDIVIDUAL i at t=1
    z[i,1]~ dcat(tr0[1:4]) 
    # CONDITIONAL LIKELIHOOD for i at t=2
    for(s in 1:T2[1]){ #loop through secondary periods
       y[i,s,1] ~ dcat(em[,z[i,1],1,s,i])
    }
    # HMM PROCESS for t>1
    for(t in 2:T){ 
      # LATENT STATE 
      z[i,t] ~ dcat(tr[1:4, z[i,t-1], t-1])
      # EMISSION
      for(s in 1:T2[t]){ # loop through secondary periods
        y[i,s,t] ~ dcat(em[,z[i,t],t,s,i])
      }
    } #t
  } # i
} #end model
"
# save jags text/code in object jags.txt.3 to local disk as a file to import into external program JAGS
modname <- "JAGS_hierarchical_pcrd.JAG"
sink(file=modname) # this creates a new file called 'JAGS_HierBayes.JAG'
cat(jags.txt,fill=TRUE) # this places the text from object jags.txt.3 into file JAGS_HierBayes.JAG
sink() # this cloes the connection to file 'JAGS_HierBayes.JAG'

# HYPERPRIORS: half-student-t (tau (precision, nu (degrees-of-freedom))
pr.sigma.eps<-c(0.05^(-2), 3) 
pr.sigma.p.t<-c(0.05^(-2), 3)
pr.sigma.p.tt<-c(0.05^(-2), 3)
# PRIORS: logit-normal
pr.lp.mu<-c(0,1.5^(-2))
pr.lphi.mu<-c(logit(0.92),1.5^(-2))
pr.g1.mu<-c(0,1^(-2))
pr.g2.mu<-c(0,1.5^(-2))
# priors on recruitment: dunif
pr.psi<-matrix(c(0,1),nrow=T,ncol=2,byrow=TRUE,dimnames=list(1:T,c("a","b")))

# DATA AUGMENTATION
n.aug <- n # a N
m <- n.aug+n
Y.aug <- array(NA,c(m, max(T2),T))
mm <- nrow(Y.aug)
dimnames(Y.aug)[[1]] <- c(row.names(Y.tt),paste0("aug",1:n.aug))
Y.aug[1:n,,]<-Y.tt # insert observed capture-histories
Y.aug[(n+1):mm,,]<-0*Y.tt[1:n.aug,,] # all-zero capture histories
Y.jags <- Y.aug+1 # convert so 1 = nocapture and 2 = capture

# ASSEMBLE THE JAGS DATA
jags.data <- list(
    y = Y.jags,
    T=T,
    T2=T2,
    M=m,
    pr.sigma.eps=pr.sigma.eps,
    pr.sigma.p.t=pr.sigma.p.t,
    pr.sigma.p.tt=pr.sigma.p.tt,
    pr.lp.mu=pr.lp.mu,
    pr.lphi.mu=pr.lphi.mu,
    pr.g1.mu=pr.g1.mu,
    pr.g2.mu=pr.g2.mu,
    pr.psi=pr.psi)

# INITIALIZATION FUNCTION
init.func = function(){
    RET=list(
     lphi.mu=runif(1,0.87,0.96), 
     lg1.mu=logit(runif(1,0.2,0.8)),
     lg2.mu=logit(runif(1,0.2,0.8)),
     lp.mu =logit(runif(1,0.1,0.5)),
     sigma.p.t=runif(1,0.01,0.025)^0.5,
     sigma.p.tt=runif(1,0.005,0.015)^0.5,
     sigma.eps = runif(1,0.0005,0.0023)^0.5
    );
    lp.i=rnorm(m,0,sd=RET$sigma.eps);
    lp.t = rnorm(T,0,RET$sigma.p.t);
    lp.s = matrix(NA,max(T2),T);   
    for(t_ in 1:T){
        lp.s[1:T2[t_],t_]<-rnorm(T2[t_],0,RET$sigma.p.tt)
    };
    gen.z =z = generate.z.psi(y=Y.aug,T2=T2,first.capture=FALSE,z.priors = list(phi.beta=c(shape1=30,shape2=5),g1.beta=c(shape1=20,shape2=20),g2.beta=c(shape1=20,shape2=20), pd.beta=c(shape1=12,shape2=65)), exclude_=1, in_=2, out_=3, dead_=4)
    RET=c(RET, list(
                   lp.tt=lp.s,
                   lp.t = lp.t,
                   lp.i=lp.i,
                   psi=gen.z$psi,
                   z = gen.z$z
                   ));
    return(RET)}

if(nchains == 1){
    jags.inits = init.func()
} else {
    jags.inits = lapply(1:nchains, function(x) init.func())
}

# MCMC parameters
  nchains <- 3 # number of MCMC chains
  nadapt <-40000 # adaption phase
  nburn <- 40000 # discard these draws
  niter <- 300000 # length of MCMC chains
  nsamp <- 1000 # number of samples to take from each MCMC chains
  thin_ <- round(niter/nsamp) # only


# compile model
  m <- jags.model(file=modname,data=jags.data,inits=jags.inits,n.chains=nchains,n.adapt=nadapt)
# burn-in phase
  update(m,nburn) # discard nburn iterations (ensure that chains reach equilibium states)
# which variables to summarize
  variable.names <- c("sigma.phi","sigma.g1","sigma.g2","sigma.p.t","sigma.p.tt","sigma.eps","phi","g1","g2","lp.mu","lp.t","lp.tt")
  post <- coda.samples(m,variable.names=variable.names,n.iter=niter,thin=thin_) 

#######################################################################################
# PART 2: Demonstration of Random effects model for PCRD model (random effects for p by secondary periods and primary periods)
library(boot)
library(rjags)
setwd(??)

# load some helpful files
source("R_PCRD_JAGS_SOURCE.R")

MARK.file.name <- "Nicholson_pcrd_ch.inp"
T2 <- c(5, 5, 10, 5,3) # number of secondary periods per primary period
T <- length(T2) # number of primary periods
capture.histories = import.mark.inp(MARK.file.name,T2,header=FALSE) # import MARK inp
Y.tt <- capture.histories[["Y.tt"]] # get the 3D array N x T2 x T
n = nrow(Y.tt) # number of (observed) individuals 

# jags model text
jags.txt <-"model{
  # HYPERPRIORS
  sigma.p.t ~ dt(0,pr.sigma.p.t[1],pr.sigma.p.t[2]) T(0,) # controls dispersion of primary periods capture
  sigma.p.tt ~ dt(0,pr.sigma.p.tt[1],pr.sigma.p.tt[2]) T(0,) # controls dispersion of 2ndary periods capture
  # hierarchical capture process
  lp.mu ~ dnorm(pr.lp.mu[1],pr.lp.mu[2]) # mean capture (logit)
  for(t in 1:T){ # loop through primary periods
    lp.t[t] ~ dnorm(0, pow(sigma.p.t,-2))
    for(s in 1:T2[t]){ # loop through secondary periods
      lp.tt[s,t] ~ dnorm(0, pow(sigma.p.tt,-2))
    }
  }
  # PRIORS: on process parameters
  lphi.mu ~ dnorm(pr.lphi.mu[1],pr.lphi.mu[2]) # survival mean (logit scale)
  lg1.mu ~ dnorm(pr.g1.mu[1],pr.g1.mu[2]) # gamma-prime, mean (logit scale)
  lg2.mu ~ dnorm(pr.g2.mu[1],pr.g2.mu[2]) # gamma-prime-prime, mean (logit scale)
  # time-varying psi
  for(t in 1:T){
    psi[t] ~ dunif(pr.psi[t,1], pr.psi[t,2])
  }
  # shunt time-constant parameters to time-varying
  for(t in 1:(T-1)){ #loop through primary periods
    g1[t] <- ilogit(lg1.mu) #prob migrate out|out (probability)
    g2[t] <- ilogit(lg2.mu) #prob migrate in|out (probability)
    phi[t] <- ilogit(lphi.mu) # survival (probability)
  }
  # HIDDEN MARKOV MODEL: TRANSITION MATRIX
  # first state transition (pure nusance; strictly from outside-pop to part of marked-pop)
  tr0[1] <- (1-psi[1]) #remains not-yet-in-pop
  tr0[2] <- psi[1] # recruits
  tr0[3] <- 0 # cannot go to TE
  tr0[4] <- 0 # cannot go to dead
  # transition matrix from (2:T)
  for(t in 1:(T-1)){
     # from state 1 (un-recruited)
     tr[1,1,t] <- 1-psi[t+1] # does not recruit
     tr[2,1,t] <- psi[t+1] # recruits
     tr[3,1,t] <- 0 # illegal
     tr[4,1,t] <- 0 # illegal
     # from state 2 (onsite)
     tr[1,2,t] <- 0
     tr[2,2,t] <- phi[t]*(1-g2[t]) # stays onsite
     tr[3,2,t] <- phi[t]*g2[t] # leaves to TE
     tr[4,2,t] <- 1-phi[t] # dies
     # from state 3 (TE)
     tr[1,3,t] <- 0
     tr[2,3,t] <- phi[t]*(1-g1[t]) # returns onside
     tr[3,3,t] <- phi[t]*g1[t] # stays in TE
     tr[4,3,t] <- 1-phi[t] # dies
     # death
     tr[1,4,t] <- 0
     tr[2,4,t] <- 0
     tr[3,4,t] <- 0
     tr[4,4,t] <- 1 # can only stay dead
  } #t
  # BUILD EMISSION MATRIX (for educational purposes, there are better ways of going this)
  for(t in 1:T){ # loop through primary periods
      for(s in 1:T2[t]){ # loop secondary periods
        # unrecruited: state 1
        em[1,1,t,s] <- 1 # no capture
        em[2,1,t,s] <- 0 # capture illegal
        # onsite: state 2
        em[2,2,t,s] <- 1/(1+exp(-1*(lp.mu+lp.t[t]+lp.tt[s,t]))) # capture
        em[1,2,t,s] <- 1-em[2,2,t,s] # no capture
        # TE: state 3
        em[1,3,t,s] <- 1 # 100% no capture
        em[2,3,t,s] <- 0 # 
        # dead: state 4
        em[1,4,t,s] <- 1 # no capture
        em[2,4,t,s] <- 0 # 
      } # s 
  } # t
  for(i in 1:M){
    # LATENT STATE FOR INDIVIDUAL i at t=1
    z[i,1]~ dcat(tr0[1:4]) 
    # CONDITIONAL LIKELIHOOD for i at t=2
    for(s in 1:T2[1]){ #loop through secondary periods
       y[i,s,1] ~ dcat(em[,z[i,1],1,s])
    }
    # HMM PROCESS for t>1
    for(t in 2:T){ 
      # LATENT STATE 
      z[i,t] ~ dcat(tr[1:4, z[i,t-1], t-1])
      # EMISSION
      for(s in 1:T2[t]){ # loop through secondary periods
        y[i,s,t] ~ dcat(em[,z[i,t],t,s])
      }
    } #t
  } # i
} #end model
"
# save jags text/code in object jags.txt.3 to local disk as a file to import into external program JAGS
modname <- "JAGS_hierarchical_pcrd.JAG"
sink(file=modname) # this creates a new file called '"JAGS_hierarchical_pcrd.JAG"
cat(jags.txt,fill=TRUE) # this places the text from object jags.txt
sink() # this cloes the connection to file 

# HYPERPRIORS: half-student-t (tau (precision, nu (degrees-of-freedom))
pr.sigma.p.t<-c(0.05^(-2), 3)
pr.sigma.p.tt<-c(0.05^(-2), 3)
# PRIORS: logit-normal
pr.lp.mu<-c(0,1.5^(-2))
pr.lphi.mu<-c(logit(0.92),1.5^(-2))
pr.g1.mu<-c(0,1^(-2))
pr.g2.mu<-c(0,1.5^(-2))
# priors on recruitment: dunif
pr.psi<-matrix(c(0,1),nrow=T,ncol=2,byrow=TRUE,dimnames=list(1:T,c("a","b")))

# DATA AUGMENTATION
n.aug <- n # a N
m <- n.aug+n
Y.aug <- array(NA,c(m, max(T2),T))
mm <- nrow(Y.aug)
dimnames(Y.aug)[[1]] <- c(row.names(Y.tt),paste0("aug",1:n.aug))
Y.aug[1:n,,]<-Y.tt # insert observed capture-histories
Y.aug[(n+1):mm,,]<-0*Y.tt[1:n.aug,,] # all-zero capture histories
Y.jags <- Y.aug+1 # convert so 1 = nocapture and 2 = capture

# ASSEMBLE THE JAGS DATA
jags.data <- list(
    y = Y.jags,
    T=T,
    T2=T2,
    M=m,
    pr.sigma.p.t=pr.sigma.p.t,
    pr.sigma.p.tt=pr.sigma.p.tt,
    pr.lp.mu=pr.lp.mu,
    pr.lphi.mu=pr.lphi.mu,
    pr.g1.mu=pr.g1.mu,
    pr.g2.mu=pr.g2.mu,
    pr.psi=pr.psi)

# INITIALIZATION FUNCTION
init.func = function(){
    RET=list(
     lphi.mu=runif(1,0.87,0.96), 
     lg1.mu=logit(runif(1,0.2,0.8)),
     lg2.mu=logit(runif(1,0.2,0.8)),
     lp.mu =logit(runif(1,0.1,0.5)),
     sigma.p.t=runif(1,0.01,0.025)^0.5,
     sigma.p.tt=runif(1,0.005,0.015)^0.5
    );
    lp.t = rnorm(T,0,RET$sigma.p.t);
    lp.s = matrix(NA,max(T2),T);   
    for(t_ in 1:T){
        lp.s[1:T2[t_],t_]<-rnorm(T2[t_],0,RET$sigma.p.tt)
    };
    gen.z =z = generate.z.psi(y=Y.aug,T2=T2,first.capture=FALSE,z.priors = list(phi.beta=c(shape1=30,shape2=5),g1.beta=c(shape1=20,shape2=20),g2.beta=c(shape1=20,shape2=20), pd.beta=c(shape1=12,shape2=65)), exclude_=1, in_=2, out_=3, dead_=4)
    RET=c(RET, list(
                   lp.tt=lp.s,
                   lp.t = lp.t,
                   psi=gen.z$psi,
                   z = gen.z$z
                   ));
    return(RET)}

# MCMC parameters
  nchains <- 2 # number of MCMC chains
  nadapt <-5000 # adaption phase
  nburn <- 40000 # discard these draws
  niter <- 300000 # length of MCMC chains
  nsamp <- 1000 # number of samples to take from each MCMC chains
  thin_ <- round(niter/nsamp) # only

if(nchains == 1){
    jags.inits = init.func()
} else {
    jags.inits = lapply(1:nchains, function(x) init.func())
}


# compile model
  m <- jags.model(file=modname,data=jags.data,inits=jags.inits,n.chains=nchains,n.adapt=nadapt)
# burn-in phase
  update(m,nburn) # discard nburn iterations (ensure that chains reach equilibium states)
# which variables to summarize
  variable.names <- c("sigma.phi","sigma.p.tt","sigma.p.t","phi","g1","g2")
  post <- coda.samples(m,variable.names=variable.names,n.iter=niter,thin=thin_) 
