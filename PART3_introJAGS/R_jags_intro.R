# # PRACTICAL EXERCISE: Introduction to JAGS
# as part of the 2017 SMM Workshop "Bayesian Capture-Recapture" by Dr. R.W. Rankin and Mrs. K.E. Nicholson
# There are three exercises
# - Exercise 1: A simple coin-flip survival model with Beta priors
# - Exercise 2: A simple coin-flip survival model with logit-Normal priors
# - Exercise 3: A logistic regression model
####################################################################################
# EXERCISE 1: Mean survival 
# This model will be a simple ''coin-flip'' estimation model. We have 30 observations of animals dying or surviving. The goal is to build a simple Bernoulli model and get the posterior survival rate, assumin our Prior probability for phi (survival) was Beta(9,1). Because we are assuming a Bernoulli distribution for each outcome, this implies that the likelihood (data density) for each observation will be y[i] ~ dbern(phi). 

library(rjags)

# DATA: counts of animals surviving (1) or dying (0) in 1 year
y <- c(1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1)

# PRIOR: Beta distribution of phi, with shape1 and shape2
pr.phi <- c(shape1 = 9, shape2=1) # prior expectation of 0.9

# jags model syntax
jags.txt <- "model{
# priors
phi ~ dbeta(pr.phi[1],pr.phi[2])
# likelihood
for(i in 1:length(y)){ # loop through each observation
  y[i] ~ dbern(phi)
}
}"
# file to save the JAGS model
jags.model.file <- "JAGS_bernoulli.JAG"
sink(file = jags.model.file) # open connection
cat(jags.txt,fill=TRUE) # send text to file
sink() # close the connection

# collect the jags data: a named list
jags.data <- list(
    y = y,
    pr.phi = pr.phi
    )
# jags expects the list to have names corresponding to all the data-items in the JAGS syntax

# initial values for MCMC chain:
jags.inits.f <- function(){
    phi = rbeta(1, pr.phi[1],pr.phi[2])
    # return a named list with the same names as the random variables in the JAGS model
    ret <- list(phi = phi)
    return(ret)
}

# MCMC parameters
n.chains <- 3 # number of MCMC chains
n.adapt <- 1000 # adaption phase for MCMC sampler
n.burn <- 1000 # burn-in phase to ensure convergence of chains
n.iter <- 10000 # number of MCMC iterations per chain
thin_ <- 10 # thin the chain to reduce auto-correlation

# COMPILE THE JAGS MODEL
m <- jags.model(file = jags.model.file,data=jags.data, inits = jags.inits.f, n.chains=n.chains, n.adapt=n.adapt)
# BURN-IN PHASE
update(m,n.burn)
# GET POSTERIOR SAMPLES
post <- coda.samples(model=m, variable.names = c("phi"), n.iter=n.iter, thin=thin_)
# (optional) collapse the three chains into a chain matrix
post.matrix <- do.call("rbind",post)

# inspect chains
plot(post)
# summarize the posterior distributions
summary(post)
# posterior expected value
phi.mean <- colMeans(post.matrix)
# posterior 95%CI
phi.95CI <- quantile(post.matrix[,"phi"], c(0.025, 0.975))
# standard error
phi.se <- sqrt(var(post.matrix[,"phi"]))
# posterior mode (most likely value) aka MAP estimate
phi.density<-density(post.matrix[,"phi"])
phi.MAP <- phi.density$x[which.max(phi.density$y)]

# plot the results
plot(phi.density,main=expression(phi))
abline(v=phi.MAP,col="blue") # MAP
abline(v=phi.mean,col="red") # expected value
abline(v = phi.95CI,col="grey")


########################################################
# EXERCISE TWO: Repeat the above, but instead with a LOGIT-NORMAL prior for phi
library(rjags)

# DATA: counts of animals surviving (1) or dying (0) in 1 year
y <- c(1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1)

# PRIOR: Beta distribution of phi, with shape1 and shape2
pr.phi <- ??????????

# jags model syntax
jags.txt <- "model{
# priors
phi ~ ????????????????
# likelihood
for(i in 1:length(y)){ # loop through each observation
  y[i] ~ dbern(phi)
}
}"
# file to save the JAGS model
jags.model.file <- "JAGS_logitNormal.JAG"
sink(file = jags.model.file) # open connection
cat(jags.txt,fill=TRUE) # send text to file
sink() # close the connection

# collect the jags data: a named list
jags.data <- list(
    y = y,
    pr.phi = pr.phi
    )
# jags expects the list to have names corresponding to all the data-items in the JAGS syntax

# initial values for MCMC chain:
jags.inits.f <- function(){
    phi = rbeta(1, pr.phi[1],pr.phi[2])
    # return a named list with the same names as the random variables in the JAGS model
    ret <- list(phi = phi)
    return(ret)
}

# MCMC parameters
n.chains <- 3 # number of MCMC chains
n.adapt <- 1000 # adaption phase for MCMC sampler
n.burn <- 1000 # burn-in phase to ensure convergence of chains
n.iter <- 10000 # number of MCMC iterations per chain
thin_ <- 10 # thin the chain to reduce auto-correlation

# COMPILE THE JAGS MODEL
m <- jags.model(file = jags.model.file,data=jags.data, inits = jags.inits.f, n.chains=n.chains, n.adapt=n.adapt)
# BURN-IN PHASE
update(m,n.burn)
# GET POSTERIOR SAMPLES
post <- coda.samples(model=m, variable.names = c("phi"), n.iter=n.iter, thin=thin_)
# (optional) collapse the three chains into a chain matrix
post.matrix <- do.call("rbind",post)

#############################################################################
# EXERCISE THREE: A Logistic Regression
# Now, instead of just a mean survival parameter, we will make a linear-regression equation based on two independent covariates (x1 and x2). Each individual has a value for x1 and x2. The regression coefficients (beta1 and beta2) will have normal priors. The task for you is to specify the priors for the regression coefficients beta1 and beta2.
# TIP: if beta1 = 0, then x1 has NO MARGINAL EFFECT ON SURVIVAL.
library(rjags)

# DATA: counts of animals surviving (1) or dying (0) in 1 year
y <- c(1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1)
# independent variable x1
x1 <- c(1.29,1.61,0.79,0.96,0.52,-1.14,-1.1,-0.65,-1.69,-0.62,1.09,-0.47,-1.74,0.1,-1.28,0.7,0.19,1.66,-0.51,-0.71,0.73,0.68,-1.13,-0.36,-1.14,0.06,1.6,-0.2,0.33,0.44)
# independent variable x2
x2 <- c(1.9,1.63,0.43,1.57,-0.98,1.29,1.14,-0.58,1.01,1.68,0.09,1.19,1.23,1.18,-0.74,-0.15,1.98,1.94,-0.77,-0.15,-0.19,-0.71,0.25,1.34,1.22,1.62,-0.02,0.27,0.36,-0.45)
# PRIORS: need three priors: for beta1, beta2 and the intercept beta0
????

# jags model syntax
jags.txt <- "model{
# priors
beta0 ~ ??????
beta1 ~ ??????
beta2 ~ ??????
# likelihood
for(i in 1:length(y)){ # loop through each observation
    # linear regression on logit scale
   logit.mu[i] <- beta0 + beta1*x1[i] + beta2*x2[i]
   # convert to probability scale
   phi[i] <- 1/(1+exp(-logit.mu[i]))
   # likelihood
   y[i] ~ dbern(phi[i])
}
}"
# file to save the JAGS model
jags.model.file <- "JAGS_logistic_regression.JAG"
sink(file = jags.model.file) # open connection
cat(jags.txt,fill=TRUE) # send text to file
sink() # close the connection

# collect the jags data: a named list
jags.data <- list(
    y = y,
    x1=x1,
    x2=x2,
    pr.beta0 = ???,
    pr.beta1 = ???,
    pr.beta2 = ???
    )
# jags expects the list to have names corresponding to all the data-items in the JAGS syntax

# initial values for MCMC chain:
jags.inits.f <- function(){
    beta0=rnorm(1,2,0.2)
    beta1=rnorm(1,0,0.2)
    beta2=rnorm(1,0,0.2)    
    # return a named list with the same names as the random variables in the JAGS model
    ret <- list(beta0=beta0,
                beta1=beta1,
                beta2=beta2)
    return(ret)
}

# MCMC parameters
n.chains <- 3 # number of MCMC chains
n.adapt <- 1000 # adaption phase for MCMC sampler
n.burn <- 1000 # burn-in phase to ensure convergence of chains
n.iter <- 10000 # number of MCMC iterations per chain
thin_ <- 10 # thin the chain to reduce auto-correlation

# COMPILE THE JAGS MODEL
m <- jags.model(file = jags.model.file,data=jags.data, inits = jags.inits.f, n.chains=n.chains, n.adapt=n.adapt)
# BURN-IN PHASE
update(m,n.burn)
# GET POSTERIOR SAMPLES
post <- coda.samples(model=m, variable.names = c("phi"), n.iter=n.iter, thin=thin_)
# (optional) collapse the three chains into a chain matrix
post.matrix <- do.call("rbind",post)

