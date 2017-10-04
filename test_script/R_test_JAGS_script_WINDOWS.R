# TEST SCRIPT: for the SMM2017 Workshop on Bayesian Capture-Recapture
# By Robert W. Rankin robertw.rankin@gmail.com

# STEP1: please ensure you have already installed R and JAGS
# STEP2: then install the 'rjags' package to call JAGS from R: install.packages("rjags")
# STEP3: run all of the following lines of code, in R.
# If it generates a plot with title "Successful! You installed JAGS", then you were successful at installing JAGS and rjags.

library(rjags) # load the rjags library
# some fake survival data
y <- c(1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1)
# make the JAGS script
jags.script <- 'model{
# priors
phi ~ dbeta(pr.a, pr.b)
# likelihood
for(i in 1:length(y)){
   y[i] ~ dbern(phi)
} # i 
}'
# NAME of the jags script
jags.file <- "test_script.JAG"
# save the jags text to the jags file: the following lines just 'sink' the text to the file
sink(file=jags.file)
cat(jags.script,fill=TRUE)
sink() # 
# how you should have a file called 'test_script.JAG' in your working directory
# check that the file exists
file.exists.yn <- file.exists(jags.file)
if(file.exists.yn){ print("Congratulations. You made a JAGS file")
} else { print("Failed. You didn't make the JAGS file. The model won't run") }
# assemble the data into an R list for JAGS
jags.data <- list(y = y, # observed data
                  pr.a = 6, # prior for phi (beta shape parameter),
                  pr.b = 1 # prior for phi (beta shape parameter),
                  )
# compile the jags model
m <- jags.model(file = jags.file, data = jags.data,n.chains=1,n.adapt=1000)
# burn-in phase of MCMC sampling
update(m, 2000)
# sample from the posterior
posteriors <- coda.samples(m,variable.names="phi",n.iter=10000,thin=5)
plot(posteriors,main=c("Successful!","You installed JAGS"))
# ... if it generates a plot with title "Successful! You installed JAGS", then you were successful in installing JAGS and rjags.
