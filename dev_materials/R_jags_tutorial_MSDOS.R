# MURUG Bayesian Workshop, Nov 19th 2015
# Robert W. Rankin
# robertw.rankin@gmail.com

### EXAMPLE 1: model of the mean
# simple model: heights are distributed according to: y ~ Normal(mu,tau)

# data
y <- c(183.46, 182.32, 178.31, 181.36, 165.12, 185.68, 170.47, 178.11, 174.86, 182.03, 180.09, 172.88, 177.94, 177.26, 182.58, 171, 173.74, 177.78, 180.02, 163.05)

library(rjags) # load jags

model.txt <- '
model{
# priors
mu ~ dnorm(100,0.000001) # vague
tau ~ dgamma(0.001, 0.001)
# likelihodd
for(i in 1:20){
  y[i]~dnorm(mu, tau)
}
}'
sink("murug_example1.JAGS")
cat(model.txt,fill=TRUE)
sink() # save the file locally (for JAGS)

# set up JAGS input arguments
jags.data <- list(y=y) # named list of data
jags.inits <- list(mu = rnorm(1,175,3), tau = 1/(runif(1,2,5))^2) # list of initializations
nchains <- 3 # number of MCMC chains
nadapt <- 1000 # number of iterations of adapting the MCMC algorithm
nburn <- 1000 # burnin phas
niter <- 10000 # number of iterations of MCMC
thin <- 5 # take sample every 5th iteration
# compile and run model (adaption phas#)
m <- jags.model(file="murug_example1.JAGS", data=jags.data, inits=jags.inits, n.chains=nchains, n.adapt=nadapt)
update(m,nburn) # burn in phase
# take samples from the posterior
post.samp <- coda.samples(m, variable.names=c("mu","tau"),n.iter=niter,thin=thin)

# concatenate the 3 chains into one chain
# post.samp.cat is a MATRIX of posterior samples: a column for mu, and a colum for sigma
post.samp.cat <- do.call("rbind",post.samp) # convert to a single matrix of samples
head(post.samp.cat)

# check convergence and good mixing
plot(post.samp) # look for NO trend line, no clustering
acf(post.samp.cat) # look for low autocorrelation (want independent draws)
gelman.diag(post.samp) # should be ~1, and < 1.1

# How do you want to use the posterior?
# get some point estimates
# posterior expectation mean, S.E. and 95%Credibility Intervals / Highest Posterior Density

post.point = c("Post.mean" = mean(post.samp.cat[,"mu"]), # posterior expectation
               "Post.SE" = sd(post.samp.cat[,"mu"]), # S.E. of estimate
               "HPD95lo" = HPDinterval(as.mcmc(post.samp.cat[,"mu"]),prob=0.95)[1], # 95% interval low
               "HPD95hi" = HPDinterval(as.mcmc(post.samp.cat[,"mu"]),prob=0.95)[2]) # 95% interval hi

# compare to the MLE
mle.model=glm(y~1); 
mle.hat=summary(mle.model)$coefficients
mle=c(mu =mle.hat[1,"Estimate"],se=mle.hat[1,"Std. Error"], lo95CI=qnorm(0.025,mle.hat[1,"Estimate"],mle.hat[1,"Std. Error"]), hi95CI=qnorm(0.975,mle.hat[1,"Estimate"],mle.hat[1,"Std. Error"]))

print(post.point)
print(mle)



###########################################################
# EXAMPLE TWO: regression
# same data, but now we have regression information
# data: y heights, x1 and x2 are two external covariates
# we want to know the MARGINAL effect of x1,x2 on the y

# data
library(rjags)
y <- c(183.46, 182.32, 178.31, 181.36, 165.12, 185.68, 170.47, 178.11, 174.86, 182.03, 180.09, 172.88, 177.94, 177.26, 182.58, 171, 173.74, 177.78, 180.02, 163.05)
x1 <- c(0.2, 1.22, 1.37, 1.05, -0.8, 0.13, 0.86, 1.73, -1.6, 1.01, -0.46, -0.82, -1.07, -1.43, 2.31, -1.02, -0.33, -0.07, 0.58, -1.83)
x2 <- c(49.77, 48.39, 49.31, 48.77, 50.01, 51.93, 51.2, 48.5, 49.97, 51.26, 50.55, 49.97, 49.91, 51.08, 51.22, 49.04, 51.1, 50.31, 51.35, 50.95)


model2.txt <- '
model{
# priors
intercept ~ dnorm(0,0.000001)
beta1 ~ dnorm(0, 0.0000001)
beta2 ~ dnorm(0, 0.0000001)
tau ~ dgamma(0.001,0.001)
# likelihood
for(i in 1:20){
  mu[i] <- intercept + beta1*x1[i] + beta2*x2[i]
  y[i] ~ dnorm(mu[i], tau)
}
}'
sink("murug_example2.JAGS")
cat(model2.txt,fill=TRUE)
sink() # save the file locally (for JAGS)

# set up JAGS input arguments
jags.data <- list(y=y,x1=x1,x2=x2) # named list of data
jags.inits <- list(intercept=0, beta1=0,beta2=0, tau = 1/(runif(1,2,5))^2) # list of initializations
nchains <- 3 # number of MCMC chains
nadapt <- 1000 # number of iterations of adapting the MCMC algorithm
nburn <- 1000 # burnin phas
niter <- 10000 # number of iterations of MCMC
thin <- 5 # take sample every 5th iteration
# compile and run model (adaption phas#)
m <- jags.model(file="murug_example2.JAGS", data=jags.data, inits=jags.inits, n.chains=nchains, n.adapt=nadapt)
update(m,nburn) # burn in phase
# take samples from the posterior
post.samp <- coda.samples(m, variable.names=c("intercept","beta1","beta2","sigma"),n.iter=niter,thin=thin)
post.samp.cat <- do.call("rbind",post.samp) # convert to a single matrix of samples

# check convergence and good mixing
plot(post.samp) # look for NO trend line, no clustering
acf(post.samp.cat) # look for low autocorrelation (want independent draws)
gelman.diag(post.samp) # should be ~1, and < 1.1

# LOOKS TERRIBLE! Poor Mixing, Did NOT converge
# NEED LONGER CHAINS

nchains <- 2 # number of MCMC chains
nadapt <- 10000 # number of iterations of adapting the MCMC algorithm
nburn <- 10000 # burnin phas
niter <- 1000000 # number of iterations of MCMC
thin <- 100 # take sample every 5th iteration
m <- jags.model(file="murug_example2.JAGS", data=jags.data, inits=jags.inits, n.chains=nchains, n.adapt=nadapt)
update(m,nburn) # burn in phase
# take samples from the posterior
post.samp <- coda.samples(m, variable.names=c("intercept","beta1","beta2","sigma"),n.iter=niter,thin=thin)
plot(post.samp) # look for NO trend line, no clustering

# compare to GLM
summary(glm(y~x1+x2))
summary(post.samp)[[1]]
