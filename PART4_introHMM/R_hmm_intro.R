# PRACTICAL EXERCISE: INTRODUCTION TO HIDDEN MARKOV MODELS (FOR CMR)
# as part of the 2017 SMM Workshop "Bayesian Capture-Recapture" by Dr. R.W. Rankin and Mrs. K.E. Nicholson

# In this exercise, we will model a single capture history as a three-state HMM.
# It is reminiscient of the POPAN model, with:
# - a ``recruitment" probability psi (time-constant) to model state1->state2
# - survival probability phi (time-constant) to model state2->state3 
# - a capture probability p (time-constant) 
# NOTE: * Real POPAN models would always have time-varying recruitment probabilities
# THIS IS FOR EDUCATIONAL PURPOSES ONLY!!!
library(rjags)
library(boot)
# the data
y<-c(0,0,0,0,1,0,1,1,0,0,0)
T<-11 # number of capture periods

# PRIORS (on psi,phi,p)
pr.p=c(1,1)  # flat Beta prior on capture probability
pr.phi=c(1,1)# flat Beta prior on survival
pr.psi=c(1,1) # flat Beta prior on recruitment

# JAGS SYNTAX
jags.model.txt <- "model{
# PRIORS
p ~ dbeta(pr.p[1], pr.p[2]) # prior on capture history
phi ~ dbeta(pr.phi[1], pr.phi[2]) # prior on survival
psi ~ dbeta(pr.psi[1], pr.psi[2]) # prior on recruitment

# HMM TRANSITION MATRIX
# FROM unborn (col1) to...
tr[1,1] <- 1-psi # unborn to unborn
tr[2,1] <- psi   # unborn to alive
tr[3,1] <- 0     # (illegal)
# FROM alive (col2) to...
tr[1,2] <- 0     # (illegal)
tr[2,2] <- phi   # alive to alive
tr[3,2] <- 1-phi # alive to dead
# FROM dead (col3) to...
tr[1,3] <- 0     # (illegal)
tr[2,3] <- 0     # (illegal)
tr[3,3] <- 1     # dead to dead

# HMM EMISSION MATRIX
# state 1: unborn (100% no capture)
em[1,1]<-1  
em[2,1]<-0  
# state 2: alive
em[1,2]<-1-p
em[2,2]<-p  
# state 3: dead (100% no capture)
em[1,3]<-1  
em[2,3]<-0

# HMM LATENT STATE PROCESS: at t=1
z[1] ~ dcat(tr[,1]) # initialize state 1 
# CONDITIONAL LIKELIHOOD: at t=1
y[1] ~ dcat(em[,z[1]]) # conditional capture

# loop though capture periods > 1
for(t in 2:T){
   z[t] ~ dcat(tr[,z[t-1]]) # z_t | z_t-1
   # CONDITIONAL LIKELIHOOD: for t>1
   y[t] ~ dcat(em[,z[t]]) # conditional capture
} # t
}"
jags.file <- "JAGS_hmm_intro.JAG"
# save the JAGS model syntax to a local file
sink(file=jags.file) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection

# ASSEMBLE THE DATA: into a list
y.jags <- y+1 # convert (0,1) into (1,2) for JAGS model
jags.data <- list(
    y = y.jags, # capture data
    T=T,
    pr.psi=pr.psi, 
    pr.phi=pr.phi, 
    pr.p=pr.p 
)

# INITIALIZE ALL THE RANDOM VARIABLES (include z)
# a function to generate random values and see the MCMC chains
# MUST return a list with names of random variables
jags.inits.f <- function(){
    # initial CMR parameters
    phi=rbeta(1,10,7)
    psi=rbeta(1,7,10)
    p=rbeta(1,10,10)
    # initialize latent states (z)
    # NOTE: must be consistent with y
    y<-c(0,0,0,0,1,0,1,1,0,0,0)+1
    unborn_<-1; alive_<-2; dead_<-3
    max.t <- max(which(y==2))+1 # first capture -1
    min.t <- min(which(y==2))-1 # last capture +1
    z <- rep(alive_,length(y))
    z[1:min.t]<- unborn_
    z[max.t:length(y)]<- dead_
    return(list(phi=phi,
                psi=psi,
                p=p,
                z=z)
           )
    }

# COMPILE THE JAGS MODEL
m1 <- jags.model(file=jags.file, inits=jags.inits.f, data=jags.data,n.chains=3,n.adapt=1000)
# BURN-IN PHASE
update(m1, 5000)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m1,variable.names=c("p","phi","psi","z"), n.iter=80000,thin=400)
# collapse the 3 chains into a single matrix
post.matrix <- do.call("rbind",post) # (optional, for easily manipulation)

# PLOT MCMC CHAIN: look for convergence, good mixing of p, psi, phi
plot(post,ask=TRUE)
# get summaries of posteriors
post.sum <- summary(post)
# print mean and sd of the parameters
sprintf("for phi: the posterior mean and SE are %0.2f and %0.2f",post.sum[[1]]["phi","Mean"],post.sum[[1]]["phi","SD"])
sprintf("for p: the posterior mean and SE are %0.2f and %0.2f",post.sum[[1]]["p","Mean"],post.sum[[1]]["p","SD"])
sprintf("for psi: the posterior mean and SE are %0.2f and %0.2f",post.sum[[1]]["psi","Mean"],post.sum[[1]]["psi","SD"])

# EDUCATIONAL PURPOSES ONLY: Let's look at the latent states, just for fun (normally, you wouldn't care about their values: just a nuissance process)
# get just the posteriors of the latent states
post.z <- post.matrix[,c('z[1]','z[2]','z[3]','z[4]','z[5]','z[6]','z[7]','z[8]','z[9]','z[10]','z[11]')] # smarter: post.matrix[,grep("z\\[",colnames(post.matrix))]
# get the percentage of time spent in each state (1,2,3) per primary periods
z.percentages<-apply(post.z,2,FUN=function(zcol){
    prop.table(tabulate(zcol, nbins=3))
})
# make a composition time-series plot
barplot(z.percentages,legend.text = c("Unborn","Alive","Dead"),col=c("pink","green","grey30"),ylab="State Marginal Probabiilities", xlab="Primary Periods")

# Question: Why is there a 100% marginal probability of being alive in primary periods 5,6,7,8?



###################################################################################################
# PART 2: PARTICIPANT EXERCISE: Mastering JAGS HMMs and 'for loops'
###################################################################################################
# In the next task, you must try to update the previous JAGS Syntax to:
# - 1) handle multiple individuals (many capture histories).
# - 2) time-varying psi (recruitment parameter varies per capture period)
# - 3) because psi is now time-varying, the transition matrix ('tr') will be different for each capture period 1 to T (because psi[t-1] != psi[t]). Therefore, there is not ONE tranition matrix, but T different transition matrices. Rather than literally make T different matrices, you should try to make a 3D array that concatenates all T matrices. I.e., tr[j,k,t]. Think of it as a 'book' where each slice tr[,,t] is a transition matrix for period t. Therefore, the dimensions of this array will be 3 dimensional. 

# Inputs:
# DATA: see the capture histories in capture history in the .RDS file:
# PRIORS: you don't have any prior information, so just use flat priors 

# TIPS:
# - previously, there was just a single capture history. Now, there will be M individuals. Use a 'for loop' to loop over all M individuals and evaluate each individuals' conditional likelihood and HMM latent process. In other words, generalize the previous code to handle any individual i
# - for specifying time-varying psi parameters, try to make a for loop as well. You COULD (but shouldn't) do something like...
# psi[1] ~ dbeta(pr.beta[1,1], pr.beta[1,2])
# psi[2] ~ dbeta(pr.beta[2,1], pr.beta[2,2])
# psi[3] ~ dbeta(pr.beta[3,1], pr.beta[3,2])
# ....
# ... but there is an easier way using a 'for loop'.
# - the 3D array 'tr' can also be constructed via a 'for loop'. Once again, you must merely generalize the code from the previous SECTION, and loop for all capture periods t. 


library(rjags)

# PLEASE SET THE WORKING DIRECTORY to the "???/BayesCMR_workshop/PART4_introHMM/"
setwd(???)

# the data
y<-read.table(file = "fake_popan_data.txt",sep=",")
# windows people: y<-read.table(file = "fake_popan_data_windows.txt",sep=",")
T<-ncol(y) # number of capture periods

# PRIORS (on psi,phi,p)
pr.p=c(1,1)  # flat Beta prior on capture probability 9
pr.phi=c(1,1)# flat Beta prior on survival
pr.psi=matrix(1,T,2) # flat Beta prior on recruitment
# notice that pr.psi is a MATRIX
print(pr.psi)

jags.model.txt <- "model{
# PRIORS
p ~ dbeta(pr.p[1], pr.p[2]) # prior on capture history
phi ~ dbeta(pr.phi[1], pr.phi[2]) # prior on survival
for(t in 1:T){
  psi[t] ~ dbeta(pr.psi[t,1], pr.psi[t,2]) # prior on recruitment
} # t

# HMM TRANSITION MATRIX
# FROM unborn (col1) to...
for(t in 1:T){
  tr[1,1,t] <- 1-psi[t] # unborn to unborn
  tr[2,1,t] <- psi[t]   # unborn to alive
  tr[3,1,t] <- 0     # (illegal)
   # FROM alive (col2) to...
  tr[1,2,t] <- 0     # (illegal)
  tr[2,2,t] <- phi   # alive to alive
  tr[3,2,t] <- 1-phi # alive to dead
# FROM dead (col3) to...
  tr[1,3,t] <- 0     # (illegal)
  tr[2,3,t] <- 0     # (illegal)
  tr[3,3,t] <- 1     # dead to dead
}
# HMM EMISSION MATRIX
# state 1: unborn (100% no capture)
em[1,1]<-1  
em[2,1]<-0  
# state 2: alive
em[1,2]<-1-p
em[2,2]<-p  
# state 3: dead (100% no capture)
em[1,3]<-1  
em[2,3]<-0

# HMM LATENT STATE PROCESS: at t=1
for(i in 1:M){
  z[i,1] ~ dcat(tr[,1,1]) # initialize state 1 
  # CONDITIONAL LIKELIHOOD: at t=1
  y[i,1] ~ dcat(em[, z[i,1] ]) # conditional capture
  # loop though capture periods > 1
  for(t in 2:T){
     z[i,t] ~ dcat(tr[,z[i,t-1],t]) # z_t | z_t-1
     # CONDITIONAL LIKELIHOOD: for t>1
     y[i,t] ~ dcat(em[,z[i,t]]) # conditional capture
  } # t
} # M
}"
jags.file <- "JAGS_hmm_intro.JAG"
# save the JAGS model syntax to a local file
sink(file=jags.file) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection

# ASSEMBLE THE DATA: into a list
y.jags <- y+1 # convert (0,1) into (1,2) for JAGS model
jags.data <- list(
    y = y.jags, # capture data
    T=T,
    pr.psi=pr.psi, 
    pr.phi=pr.phi, 
    pr.p=pr.p 
)


# ASSEMBLE THE DATA: into a list
y.jags <- y+1 # convert (0,1) into (1,2) for JAGS model
jags.data <- list(
    y = y.jags, # capture data
    T=T,
    M = nrow(y),
    pr.psi=pr.psi, 
    pr.phi=pr.phi, 
    pr.p=pr.p 
)

# INITIALIZE ALL THE RANDOM VARIABLES (include z)
# a function to generate random values and see the MCMC chains
# MUST return a list with names of random variables
jags.inits.f <- function(){
    # initial CMR parameters
    phi=rbeta(1,10,7)
    p=rbeta(1,10,10)
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
    return(list(phi=phi,
                psi=psi,
                p=p,
                z=z)
           )
    }

# COMPILE THE JAGS MODEL
m1 <- jags.model(file=jags.file, inits=jags.inits.f, data=jags.data,n.chains=3,n.adapt=1000)
# BURN-IN PHASE
update(m1, 5000)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m1,variable.names=c("p","phi","psi"), n.iter=80000,thin=200)
# collapse the 3 chains into a single matrix
post.matrix <- do.call("rbind",post) # (optional, for easily manipulation)

# check MCMC Convergence and mixing
gelman.diag(post) # convergence rule-of-thumb: values <1.1
par(mfrow=c(3,1));
acf(post.matrix[,"p"]); # auto-correlation (need independent samples)
acf(post.matrix[,"phi"]); # auto-correlation
acf(post.matrix[,"psi[1]"]) # auto-correlation

# TODO:
# summarize the posterior distributions of p and phi
# qu: what is a mean and 95%CI for phi?
# qu: what is the probability that phi <0.9?
# qu: make plots to compare to the true values:
# - phi.true <-  0.901
# - p.true <- 0.503

###########################################
# PART3: ADD extra POPAN-LIKE DERIVATIVES
# -- population abundance, births, probability of entries 

# Add the following lines of code to your previous JAGS file
# (take care to insert the code BEFORE the final bracket }



'# DERIVATIVES: Pop abundance and recruits
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
for(t in 2:T){  cumprob[t] <- psi[t]*prod(1-psi[1:(t-1)]) }
cumprob.norm <- sum(cumprob[1:T])
# POPAN Inclusion probabilities
for(t in 1:T){ pent[t] <- cumprob[t]/cumprob.norm } #t
}
'

# remember: to add
# compile model: 
m1 <- jags.model(file=jags.file, inits=jags.inits.f, data=jags.data,n.chains=3,n.adapt=1000)
# BURN-IN PHASE
update(m1, 5000)
# SAMPLE FROM POSTERIORS
# remember to add new Derivates to the 'variable.names' argument 
# variable.names=c("recruits","pent","N")
post <- coda.samples(m1,variable.names=c("p","phi","psi","recruits","pent","N"), n.iter=80000,thin=200)
post.matrix <- do.call("rbind", post)

# plot mcmc and do diagnostics
plot(post,ask=TRUE)
gelman.diag(post[,c("phi","p","pent[1]","pent[2]","pent[3]","pent[4]","pent[5]")])

# summary()

# Let's compare the outputs from our Bayesian model vs. PROGRAM MARK
library(RMark)
y.mark <- data.frame(ch=apply(y,1,function(row_) { paste0(row_,collapse="")}),stringsAsFactors=FALSE)
# remove ALL ZERO captures
y.mark<-y.mark[which(y.mark$ch!="00000000000"),,drop=FALSE] # mark mark data
mark.processed<-process.data(y.mark,model="POPAN") # process mark dadta
mark.ddl<-make.design.data(mark.processed) # make design data
# let's fix a few parameters
mark.ddl$pent
mle<-mark(mark.processed,mark.ddl, # run mark model
        model.parameters=list(Phi=list(formula=~1),# time-invariant survival
                              p=list(formula=~1,share=TRUE), # time-invariant capture
                              pent=list(formula=~time)), # time-varying pent
        retry=20,invisible=TRUE,adjust=FALSE,brief=TRUE,delete=TRUE)
# get mark estimates of survival and p
mle.phi <- mle$results$real["Phi g1 a0 t1",]
mle.p <- mle$results$real["p g1 a0 t1",]
# get mark estimates of pent (entry probabilities)
mle.pent <- mle$results$real[grep("pent",rownames(mle$results$real)),]
# get mark's estimates of births (notice T1 isn't reported (why?))
mle.recruits <- mle$results$derived[["B Net Births"]]
# get Mark's estimates of pent
mle.N <- mle$results$derived[["N Population Size"]]

# PLOTS: comparing MLE's to posterior point estimates
source("SOURCE_helpful_aux_functions.R") # add some plotting functions

# PLOT ABUNDANCE
par(mfrow=c(2,1))
plot.posterior.timeseries(posterior = post.matrix[,grep("N\\[",colnames(post.matrix))], main="N",xlab="Primary Periods", ylab="Counts",xaxis.lab = 1:11) # plot the posterior
plot.mle.points(mle.N) # add mark estimates

# PLOT RECRUITS
plot.posterior.timeseries(posterior = post.matrix[,grep("recruits\\[",colnames(post.matrix))[-1]],main="Apparent Births", xlab="Primary Periods", ylab="Counts",xaxis.lab = 2:11) # plot the posterior
plot.mle.points(mle.recruits) # add mark estimates

# PLOT THE ENTRY PROBABILITIES
plot.posterior.timeseries(posterior = post.matrix[,grep("pent\\[",colnames(post.matrix))[-1]], main="Entry Probabilities pent",xlab="Primary Periods", ylab="probability",xaxis.lab = 2:11) # plot the posterior
# add mark
plot.mle.points(mle.pent)
# Qu: which pent and birth values have BOUNDARY-VALUE MLEs? (i.e., close to their natural boundaries, like 0), and what is the consequence of this on the MLEs versus the posteriors

# plot the history of recruitment at t=3 (demonstrate difference between mode,mean, MLE)
#jpeg(filename="/tmp/POPAN_demo_recruit2.jpg",quality=50,pointsize=12,width=600,height=480)
h<-hist(post.matrix[,"recruits[3]"],breaks=15,col="grey80",main="Apparent Births at t=3",xlab="counts of births at t=3", xlim=c(0,15))
abline(v = mean(post.matrix[,"recruits[3]"]),col="blue",lwd=2); text(x=mean(post.matrix[,"recruits[3]"]),y=100, labels="Posterior mean",srt=90,col="blue",pos=4,font=2,cex=1.3)
h.r3.mode.x <-mean(h$breaks[which.max(h$counts)+c(0,1)]) # for MAP (histogram x position)
h.r3.mode.y <-h$counts[which.max(h$counts)] # for MAP (histogram y position)
text(x=h.r3.mode.x,y=h.r3.mode.y-10, labels="MAP",col="blue",pos=3,font=2,cex=1.3) # plot MAP
polygon(x=mle.recruits[2,c("lcl","lcl","ucl","ucl")],y=c(0,h.r3.mode.y,h.r3.mode.y,0),col="pink",lwd=2) # 95%\CI
text(x=mle.recruits[2,c("ucl")],y=100, labels="95% Confidence Interval",col="red",pos=4,font=2,srt=90,cex=1.3) # plot MAP
abline(v = mle.recruits[2,"estimate"],col="red",lwd=2); text(x=mle.recruits[2,"estimate"],y=100, labels="MLE",srt=90,col="red",pos=1,font=2,lwd=2,lty=2,cex=1.3)
#dev.off()
####################################################
# PART4: DATA GENERATION:
# THIS IS HOW I GENERATED THE DATA FOR PART 3

library(rjags)

# PLEASE SET THE WORKING DIRECTORY to the "???/BayesCMR_workshop/PART4_introHMM/"

# the data
T<-11 # number of capture periods

# PRIORS (on psi,phi,p)
pr.p=c(2000,2000)  # expected p=0.5
pr.phi=c(9000,1000)# expected phi=0.9
pr.psi=matrix(0,T,2) # Beta prior for time-varying psi
pr.psi[,1]<-c(400,(1:10)*5)
pr.psi[,2]<-c(600,(10:1)*20)

# JAGS SYNTAX
jags.model.txt <- "model{
# PRIORS
p ~ dbeta(pr.p[1], pr.p[2]) # prior on capture history
phi ~ dbeta(pr.phi[1], pr.phi[2]) # prior on survival
for(t in 1:T){
  psi[t] ~ dbeta(pr.psi[t,1], pr.psi[t,2]) # prior on recruitment
}
# HMM TRANSITION MATRIX
# FROM unborn (col1) to...
for(t in 1:T){
tr[1,1,t] <- 1-psi[t] # unborn to unborn
tr[2,1,t] <- psi[t]   # unborn to alive
tr[3,1,t] <- 0     # (illegal)
# FROM alive (col2) to...
tr[1,2,t] <- 0     # (illegal)
tr[2,2,t] <- phi   # alive to alive
tr[3,2,t] <- 1-phi # alive to dead
# FROM dead (col3) to...
tr[1,3,t] <- 0     # (illegal)
tr[2,3,t] <- 0     # (illegal)
tr[3,3,t] <- 1     # dead to dead
}
# HMM EMISSION MATRIX
# state 1: unborn (100% no capture)
em[1,1]<-1  
em[2,1]<-0  
# state 2: alive
em[1,2]<-1-p
em[2,2]<-p  
# state 3: dead (100% no capture)
em[1,3]<-1  
em[2,3]<-0
# HMM LATENT STATE PROCESS: at t=1
for(i in 1:M){
  z[i,1] ~ dcat(tr[,1,1]) # initialize state 1 
  # CONDITIONAL LIKELIHOOD: at t=1
  y[i,1] ~ dcat(em[,z[i,1]]) # conditional capture
  # loop though capture periods > 1
  for(t in 2:T){
     z[i,t] ~ dcat(tr[,z[i,t-1],t]) # z_t | z_t-1
     # CONDITIONAL LIKELIHOOD: for t>1
     y[i,t] ~ dcat(em[,z[i,t]]) # conditional capture
   } # t
 }# M
}"
jags.file <- "JAGS_hmm_exercise2.JAG"
# save the JAGS model syntax to a local file
sink(file=jags.file) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection

# ASSEMBLE THE DATA: into a list
jags.data <- list(
    y = matrix(NA,350,T), # capture data
    T=T,
    M=350,
    pr.psi=pr.psi, 
    pr.phi=pr.phi, 
    pr.p=pr.p 
)

# INITIALIZE ALL THE RANDOM VARIABLES (include z)
# a function to generate random values and see the MCMC chains
# MUST return a list with names of random variables

# COMPILE THE JAGS MODEL
m1 <- jags.model(file=jags.file, inits=NULL, data=jags.data,n.chains=1,n.adapt=1000)
# BURN-IN PHASE
update(m1, 1000)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m1,variable.names=c("p","phi","psi","y"), n.iter=1000,thin=500)
# collapse the 3 chains into a single matrix
y <- matrix(0,350,T)
for(k in colnames(post[[1]])[grep("y\\[",colnames(post[[1]]))]){
    eval(parse(text = paste0(k,"<-",post[[1]][2,k])))
}
first <- apply(y, 1, function(x) min(which(x==2)))
y <- y[order(first),]-1    

write.table(y,file="fake_popan_data.txt",row.names=FALSE,col.names=FALSE,sep=",",eol="\r\n")

# true values
# psi.tru<-c(0.374,0.033,0.04,0.059,0.106,0.2,0.264,0.279,0.354,0.597,0.793)
# phi.tru <- c(0.901)
# p.tru <- c(0.513)
# pent.true <- do.call(function(psi){ c(psi[1],sapply(2:11,function(k) { prod(1-psi[1:(k-1)])*psi[k]})) },args=list(psi=psi.tru))
# c(0.374,0.033,0.04,0.059,0.106,0.2,0.264,0.279,0.354,0.597,0.793)
