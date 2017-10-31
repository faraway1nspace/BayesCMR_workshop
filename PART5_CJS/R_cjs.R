# PRACTICAL EXERCISE PART 5: Cormack-Jolly-Seber in JAGS
# as part of the 2017 SMM Workshop "Bayesian Capture-Recapture" by Dr. R.W. Rankin and Mrs. K.E. Nicholson
# There are TWO DEMONSTRATIONS below and TWO EXERCISES
# DEMONSTRATION 1: phi(sex)p(sex) with beta priors
# DEMONSTRATION 2: phi(sex)p(sex) with logit-normal priors
# EXERCISE 1: phi(t,sex)p(t,sex) with beta priors
# EXERCISE 2: phi(t,sex)p(t,sex) with logit-normal priors

# Data from Nicholson et al 2012. Abundance, survival and temporary emigration of bottlenose dolphins (Tursiops sp.) off Useless Loop in the western gulf of Shark Bay, Western Australia. Mar. Freshwater Res. 63, 1059–1068. doi:10.1071/MF12210 . In the original study was for a PCRD system (multiple secondary periods per primary period). For this exercise, everything has been flattened to a single capture period
# THIS IS FOR EDUCATIONAL PURPOSES ONLY!!!

##############################################################
# CJS DEMONSTRATION 1: sex-varying CJS (Beta parameterization)
# phi(sex)p(sex) 
# in this demonstration, we will run a sex-varying CJS with BETA priors on phi and p
# (in the next demonstration, we will run a sex-varying CJS with LOGIT-NORMAL priors on phi and p)

library(rjags)
setwd("") # set your working directory to BayesCMR_workshop/PART5_CJS
# READ IN THE CAPTURE HISTORIES (Nicholson, for CJS)
d <- read.csv(file = "data_nicholson_cjs.csv", header=TRUE,stringsAsFactors=FALSE) #
y <- d[,c('p1','p2','p3','p4','p5')]
T <- 5 # five capture periods
# sex vector
sex <- d[,"sex"] # females f, males m, unknowns u
sex <- ifelse(d[,"sex"]=="f",1,ifelse(d[,"sex"]=="m",2,3)) # females f, males m, unknowns u
n.sex <- length(unique(sex)) # number of sexes

# for CJS, we also need the period at FIRST CAPTURE
first <- apply(y, 1, function(x) min(which(x !=0)))

# filter: for CJS, we DO NOT model individuals whose first capture was at T
cjs.keep.ix <- which(first < T) # find those rows where first is NOT at T
y.cjs <- y[cjs.keep.ix,]
y.jags <- y.cjs +1 # fors jags: make 1=no capture and 2=capture
sex.jags <- sex[cjs.keep.ix]
first.jags <- first[cjs.keep.ix]

# priors on capture
pr.p<-matrix(c(1,1,1,1,1,1),nrow=n.sex,ncol=2,dimnames=list(c("f","m","u"),c("a","b"))) # rows are sexes, columns are the shape parameters
pr.phi<-matrix(c(1,1,1,1,1,1),nrow=n.sex,ncol=2,dimnames=list(c("f","m","u"),c("a","b")))

# model with sex differences
jags.model.txt <- "model{
for(x in 1:n.sex){ # loop through sexes
  phi[x] ~ dbeta(pr.phi[x,1],pr.phi[x,2])
  p[x] ~ dbeta(pr.p[x,1],pr.p[x,2])
  # FROM alive to...
  tr[1,1,x] <- phi[x]   # alive
  tr[2,1,x] <- 1-phi[x] # dead
  # FROM dead to ...
  tr[1,2,x] <- 0 # alive (illegal)
  tr[2,2,x] <- 1 # dead to dead
  # HMM EMISSION MATRIX
  # state 1: alive
  em[1,1,x]<-1-p[x] # miss
  em[2,1,x]<- p[x]  # capture
  # state 2: dead 
  em[1,2,x]<- 1  # miss
  em[2,2,x]<- 0  # capture
} # sex
# latent proces
for(i in 1:n){
  # HMM LATENT STATE PROCESS: at t=1
  z[i,first[i]+1] ~ dcat(tr[,1,sex.v[i]]) # initialize state 1 (alive)
  # loop though capture periods after first capture
  for(t in (first[i]+1):(T-1)){
     z[i,t+1] ~ dcat(tr[,z[i,t],sex.v[i]]) # z_t | z_t-1
  } # t
}# i 
# conditional likelihood
for(i in 1:n){
  # loop though capture periods after first capture
  for(t in (first[i]+1):T){
     y[i,t] ~ dcat(em[,z[i,t],sex.v[i]])
  } # t
} # i
# Abundance estimate: horvitz-type estimator: n_obs 
for(t in 2:T){ # loop through time
  for(x in 1:n.sex){ # loop through sex
    for(i in 1:n){ # loop through individuals
      N_i[i,t-1,x] <- equals(sex.v[i],x)*equals(y[i,t],2)/p[sex.v[i]]
    } # i
    # tally estimate for sex x at time t (offset by -1)
    N[t-1,x] <- sum(N_i[,t-1,x]) # sum over all individuals of sex x
   } # x
} # t
}"
jags.file <- "cjs_phi-dot_p-dot_beta.jag"
# save the JAGS model syntax to a local file
sink(file=jags.file) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection

# ASSEMBLE THE DATA: into a list
jags.data <- list(
    n = nrow(y.jags), # number of captures
    n.sex = n.sex, # number of sexes (3: f,m,unkn)
    y = y.jags, # capture data
    T=T, # number of primary periods
    first = first.jags, # vector of primary period at first capture for each i
    sex.v = sex.jags, #  vectof of sex for each i
    pr.phi=pr.phi, # priors on phi (per sex)
    pr.p=pr.p # priors on p (per sex)
)

# INITIALIZE ALL THE RANDOM VARIABLES (include z)
# a function to generate random values and see the MCMC chains
# MUST return a list with names of random variables
jags.inits.f <- function(){
    # initial CMR parameters
    phi=rbeta(n.sex,10,7)
    p=rbeta(n.sex,10,10)
    # initialize latent states (z)
    # NOTE: must be consistent with y.jags
    alive_<-1; dead_<-2; capture_ <- 2
    z <- matrix(0,nrow(y.jags),T+2) # add two dummy columns (one at beginning and one at the end)
    for(i in 1:nrow(y.jags)){
        max.t <- if(any(y.jags[i,]==capture_)){ 2+max(which(y.jags[i,]==capture_))} else { T+2}
        min.t <- if(any(y.jags[i,]==capture_)){ min(which(y.jags[i,]==capture_))} else { 1}        
        z[i,] <- rep(alive_,ncol(z))
        z[i,1:(min.t+1)]<- NA
        z[i,max.t:ncol(z)]<- dead_
    }
    z <- z[,2:(ncol(z)-1)]
    return(list(phi=phi,
                p=p,
                z=z )
           )
    }

# COMPILE THE JAGS MODEL
m1 <- jags.model(file=jags.file, inits=jags.inits.f, data=jags.data,n.chains=3,n.adapt=1000)
# BURN-IN PHASE
update(m1, 5000)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m1,variable.names=c("p","phi","N"), n.iter=40000,thin=400)
# collapse the 3 chains into a single matrix
post.matrix <- do.call("rbind",post) # (optional, for easier manipulation)

library(RMark)
y.mark <- data.frame(ch=apply(y.cjs,1,function(x) paste0(x,collapse="")),sex=factor(sex.jags),stringsAsFactors=FALSE)
mmark=mark(y.mark,model="CJS",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~sex)),threads=7,groups="sex")
ptable = mmark$results$real[4:6,]
# calculate the nb of recaptured individiduals / occasion
obs = list(f=colSums(y.cjs[which(sex.jags==1),-1]),
           m=colSums(y.cjs[which(sex.jags==2),-1]),
           u=colSums(y.cjs[which(sex.jags==3),-1]))
Nhat = list(f=obs[[1]]/ptable[1,1],
                    m=obs[[2]]/ptable[2,1],
                    u=obs[[3]]/ptable[3,1]) 

######################################################################
# CJS DEMONSTRATION 2: sex-varying CJS (logit-normal parameterization)
# phi(logit(sex))p(logit(sex,effort)) 
# in this demonstration, we will run a sex-varying CJS with logit-normal specification of both phi and p
# the logit normal allows complex LINEAR MODEL decompositions
# e.g., we will make logit(p) = mu + effort + sex 
# ... this is a typical logistic regression. It is different from the above Beta parameterization, in which sexes are COMPLETELY independent (no sharing of information)

library(rjags)
library(boot)
setwd("~/Documents/science/presentations/BayesCMR_workshop/PART5_CJS/")
# READ IN THE CAPTURE HISTORIES (Nicholson, for CJS)
d <- read.csv(file = "data_nicholson_cjs.csv", header=TRUE,stringsAsFactors=FALSE) #
y <- d[,c('p1','p2','p3','p4','p5')]
T <- 5 # five capture periods
# convert sex to 1,2,3 integers (for jags)
sex <- ifelse(d[,"sex"]=="f",1,ifelse(d[,"sex"]=="m",2,3)) # females f, males m, unknowns u
n.sex <- length(unique(sex))
effort <- c(1.609438,2.302585,1.609438,1.098612) # effort covariate for t2,t3,t4,t5 (t1 not necessary)

# for CJS, we also need the period at FIRST CAPTURE
first <- apply(y, 1, function(x) min(which(x !=0)))

# filer: for CJS, we DO model individuals whose first capture is at T
cjs.keep.ix <- which(first < T) # find those rows where first is NOT at T
y.cjs <- y[cjs.keep.ix,]
y.jags <- y.cjs +1 # fors jags 
sex.jags <- sex[cjs.keep.ix]
first.jags <- first[cjs.keep.ix]

# priors on capture
pr.logit.p.mu <- c(0, 1.55^(-2))  # logit-normal prior: flat on probability scale
pr.logit.p.effort <-c(0,1.55^(-2))
pr.logit.p.sex<-matrix(c(0,1.55^(-2)),2,2,byrow=TRUE,dimnames=list(c("m","u"),c("mean","precision"))) # logit-normal prior: flat on probability scale
pr.logit.phi<-matrix(c(0,1.55^(-2)),3,2,byrow=TRUE,dimnames=list(c("f","m","u"),c("mean","precision"))) # logit-normal prior: flat on probability scale

# model with sex differences
jags.model.txt <- "model{
# mean capture probability on logit scale (for females)
logit.p.mu ~dnorm(pr.logit.p.mu[1],pr.logit.p.mu[2])
# effect of effort
logit.p.effort ~dnorm(pr.logit.p.effort[1],pr.logit.p.effort[2])
# sex effect (for males and unknowns)
for(x in 1:(n.sex-1)){
  logit.p.sex[x] ~ dnorm(pr.logit.p.sex[x,1],pr.logit.p.sex[x,2])
}
# combine mean, sex, and effort effects into a single p (time-varying)
for(t in 1:(T-1)){
   # females 
   p[t,1] <- 1/(1+exp(-1*(logit.p.mu + EFFORT[t]*logit.p.effort))) # inverse logit transformation
   for(x in 1:(n.sex-1)){
     # for males and unknowns
     p[t,x+1] <- 1/(1+exp(-1*(logit.p.mu + logit.p.sex[x] + EFFORT[t]*logit.p.effort))) # ilogit
   }
}
for(x in 1:n.sex){ # loop through sexes
  # logit normal survival
  logit.phi[x] ~ dnorm(pr.logit.phi[x,1],pr.logit.phi[x,2])
  phi[x] <- 1/(1+exp(-logit.phi[x])) # logit transformation
  # FROM alive to...
  tr[1,1,x] <- phi[x]   # alive
  tr[2,1,x] <- 1-phi[x] # dead
  # FROM dead to ...
  tr[1,2,x] <- 0 # alive (illegal)
  tr[2,2,x] <- 1 # dead to dead
}
# HMM EMISSION MATRIX
for(x in 1:n.sex){
  for(t in 1:(T-1)){
    # state 1: alive
    em[1,1,t,x]<-1-p[t,x] # miss
    em[2,1,t,x]<- p[t,x]  # capture
    # state 2: dead 
    em[1,2,t,x]<- 1  # miss
    em[2,2,t,x]<- 0  # capture
  } # t
} # sex
# latent proces
for(i in 1:n){
  # HMM LATENT STATE PROCESS: at t=1
  z[i,first[i]+1] ~ dcat(tr[,1,sex.v[i]]) # initialize state 1 (alive)
  # loop though capture periods after first capture
  for(t in (first[i]+1):(T-1)){
     z[i,t+1] ~ dcat(tr[,z[i,t],sex.v[i]]) # z_t | z_t-1
  } # t
}# i 
# conditional likelihood
for(i in 1:n){
  # loop though capture periods after first capture
  for(t in (first[i]+1):T){
     y[i,t] ~ dcat(em[,z[i,t],t-1,sex.v[i]])
  } # t
} # i
# Abundance estimate: horvitz-type estimator: n_obs 
for(t in 2:T){ # loop through time
  for(x in 1:n.sex){ # loop through sex
    for(i in 1:n){ # loop through individuals
      N_i[i,t-1,x] <- equals(sex.v[i],x)*equals(y[i,t],2)/p[t-1,sex.v[i]]
    } # i
    # tally estimate for sex x at time t (offset by -1)
    N[t-1,x] <- sum(N_i[,t-1,x]) # sum over all individuals of sex x
   } # x
} # t
}"
jags.file <- "cjs_phi-dot_p-dot_logitNormal.jag"
# save the JAGS model syntax to a local file
sink(file=jags.file) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection

# ASSEMBLE THE DATA: into a list
jags.data <- list(
    n = nrow(y.jags), # number of captures
    n.sex = n.sex, # number of sexes (3: f,m,unkn)
    y = y.jags, # capture data
    T=T, # number of primary periods
    first = first.jags, # vector of primary period at first capture for each i
    sex.v = sex.jags, #  vectof of sex for each i
    EFFORT = effort,
    pr.logit.phi=pr.logit.phi, # priors on phi (per sex)
    pr.logit.p.mu =pr.logit.p.mu ,
    pr.logit.p.effort=pr.logit.p.effort,
    pr.logit.p.sex=pr.logit.p.sex
)

# INITIALIZE ALL THE RANDOM VARIABLES (include z)
# a function to generate random values and see the MCMC chains
# MUST return a list with names of random variables
jags.inits.f <- function(){
    # initial CMR parameters
    logit.phi=rnorm(n.sex,1.4,0.1)
    logit.p.mu = rnorm(1,0,0.1) # female logit capture probability
    logit.p.effort = rnorm(1,0,0.1) # female logit capture probability    
    logit.p.sex = rnorm(2,0,0.1) # males and unknowns logit capture probability
    # initialize latent states (z)
    # NOTE: must be consistent with y.jags
    alive_<-1; dead_<-2; capture_ <- 2
    z <- matrix(0,nrow(y.jags),T+2) # add two dummy columns (one at beginning and one at the end)
    for(i in 1:nrow(y.jags)){
        max.t <- if(any(y.jags[i,]==capture_)){ 2+max(which(y.jags[i,]==capture_))} else { T+2}
        min.t <- if(any(y.jags[i,]==capture_)){ min(which(y.jags[i,]==capture_))} else { 1}        
        z[i,] <- rep(alive_,ncol(z))
        z[i,1:(min.t+1)]<- NA
        z[i,max.t:ncol(z)]<- dead_
    }
    z <- z[,2:(ncol(z)-1)]
    return(list(logit.phi=logit.phi,
                logit.p.mu=logit.p.mu,
                logit.p.effort=logit.p.effort,
                logit.p.sex=logit.p.sex,
                z=z)
           )
    }

# COMPILE THE JAGS MODEL
m1 <- jags.model(file=jags.file, inits=jags.inits.f, data=jags.data,n.chains=3,n.adapt=1000)
# BURN-IN PHASE
update(m1, 5000)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m1,variable.names=c("logit.phi","logit.p.mu","logit.p.sex","logit.p.effort","p","phi","N"), n.iter=40000,thin=400)
# collapse the 3 chains into a single matrix
post.matrix <- do.call("rbind",post) # (optional, for easier manipulation)



###################################################################################################
# CJS Exercise 1

# In this exercise, you will modify the phi(sex)p(sex) model and make it a fuly-time-varying model. I.e., the script belong uses constant phi and constant p (only variation in sex). Try to edit the JAGS code to have time-varying survival (phi) and time-varying capture p.
# *NOTE: remember that: phi[T-1] and p[T] cannot be seperated
# *NOTE: remember that: there is no ``p_1'' conceptually (because the CJS conditions on first capture), BUT in the JAGS syntax we cannot have a vector that skips the first value, i.e., there must be a p[1]. So, just keep in mind that all the JAGS p's are offset by one entry.
# *HINTS: - the transition matrix WAS 3 dimensional [state,state,sex] but will now be 4-dimensional [state,state,sex,time]. Similarly for the emission matrix
# *HINTS: - previously, phi was just a vector with three values (female,male, unknown); but now it should be a MATRIX with rows representing TIME and columns representing SEX. A similar idea applies to p
# *HINTS: - previously, the priors on phi was encoded as a matrix, with rows as SEX and columns as the Beta shape parameters (a and b). Now, the priors of phi should be a 3D ARRAY with rows represent TIME, columns represent SEX and slices representing SHAPE-PARAMETERS. A similar idea applies to p.
# *HINTS: - you must enforce a constraint on the last value of phi's. Enforce phi[T-2]=phi[T-2]. Unfortunately, you CANNOT mix random (~) and deterministic (<-) elements in a single vector in JAGS. In other words, you can't do...
# phi[1] ~ dbeta(a1,b1)
# phi[2] ~ dbeta(a2,b2)
# phi[3] ~ dbeta(a3,b3)
# phi[4] <- phi[3]
# ... the best way is to make two different objects: one for random variables only, and then shunt those values into another "deterministic" object. E.g.,
# phi.rv[1] ~ dbeta(a1,b1)
# phi.rv[2] ~ dbeta(a2,b2)
# phi.rv[3] ~ dbeta(a3,b3)
# phi[1] <- phi.rv[1]
# phi[2] <- phi.rv[2]
# phi[3] <- phi.rv[3]
# phi[4] <- phi.rv[3] # NOTICE THE ENFORCEMENT OF CONSTRAINT HERE

[: COPY AND PASTE THE CODE FROM DEMONSTRATION 1 :]
[: THEN MODIFY!:]



###################################################################################################
# CJS Exercise 2
# In this exercise, you will modify the sex-varying logit-normal model phi(logit(sex))p(logit(sex)) make it a fuly-time-varying model. You will KEEP the effort:sex decomposition for p: it is already time-varying because of the use of EFFORT covariate (which is the time-varying thing). So, this exercise only involved modification of phi to make it sex and time varying.
# *NOTE: there are TWO WAYS to do a sex and time variation: i) the "main effects" model includes a main effort for time and a main effect of sex (being a male), in which case the sex effect is the same magnitude each year. ii) the "interaction" model has independent sex-effects each year, so the sex effect of being a male in period one is different from the sex effect in period 2. You must decide which you would like to do. The main effects model is more CONSTRAINED and easier to do, the INTERACTION model is likely more realistic but may be overly-complex
# *HINT: previously, the priors for phi were stored in a matrix with rows representing sex and columns representing the Normal parameters MEAN and PRECISION. Your prior structure will be different now that phi will be time-varying. In the main effects model: there is a PRIOR for the intercept (females in year 2), a prior for the main effect of being male, and a prior for the main effect of being 'unsexed/unknown', then there will be a prior for the main effect of year 3, a prior for the year effect of year 4, and a prior for the year effect of year 5. Alternatively, in the interaction model, there will be prior for EVERY COMBINATION: f-year2, f-year3, f-year4, f-year5, m-year2, m-year3, m-year4, m-year5, u-year2, u-year3, u-year4, u-year5.

[: COPY AND PASTE THE CODE FROM DEMONSTRATION 2 :]
[: THEN MODIFY!:]
