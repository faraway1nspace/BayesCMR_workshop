# This is an example multi-state closed robust design (MSCRD) model with Humpback Dolphins from Queensland, Australia, as studied by Dr. Daniele Cagnazzi and Dr. Robert Rankin Unpublished
# IMPORTANT FEATURES OF MODEL:
# 1) two sexes,
# 2) two strata (north coast and south coast) and two temporary-emigration (unobservable) strata.
# 3) a covariate (flooding) which influences the capture probability (makes the water very murky).
# 4) unequal sampling area, and so we adjust the capture probabilities by the computed sampling area per secondary period. Called PROP.AREA (proportion of study area covered)
# 5) hierarchical model of time-varying secondary period capture probabilities (in addition to covariate flooding and the PROP.AREA
# NOTE: the data has been censored for the last five years, so this is just half the actual data that will make it into a publication.

library(rjags)
library(boot)
source("SOURCE_helper_functions.R")

T<-6 # 6 primary periods
T2<-c(5,6,3,5,5,5) # secondary periods per primary period
TTT.tot <- sum(T2) # total number of seconary periods
tt.ind <- matrix(c(1,2,3,4,5,NA,6,7,8,9,10,11,12,13,14,NA,NA,NA,15,16,17,18,19,NA,20,21,22,23,24,NA,25,26,27,28,29,NA),nrow=6,ncol=T) # cumulative secondary periods, arrange with columns equal to primary periods

# read in the data
data_<-readRDS("data_cagnazzi_sousa.RDS")
# capture data
y <- data_$y
N.obs <- nrow(y) #
# augment the data for JAGS

# sex data
n.sex <- 2 # two sexes
sex.vector <- data_$sex.ind

# study area: what proportion of the study area was sampled per secondary periods?
PROP.AREA <- data_$PROP.AREA

# covariates: flooding by strata (columns) and before/after
cov.p.flooding <- data_$cov.p.flood

# PRIORS
priors <- list(
    pr.llambda_mu=c(0,0.5^-2),
                                        #gte, temp emigration rate
    pr.lgte_mu= c(logit(0.2),0.5^-2),
                                        # gaa, movement to strata A->A
    pr.lgaa_mu= c(logit(0.8),0.66^-2),
                                        # gbb, movement to strata B->B
    pr.lgbb_mu= c(logit(0.8),0.66^-2),
                                        # phi, survival
    pr.lphi_mu= c(logit(0.9129279),0.6832291^(-2)),
    pr.dlphi_sex=c(0,0.4^(-2)),
                                        # p, capture
    pr.lp_mu=c(logit(0.2),0.66^(-2)),
    pr.dlp_sex=c(0,0.66^(-2)),
    pr.dlp_flo_dur=c(0,0.66^(-2)),
    pr.dlp_flo_aft=c(0,0.66^(-2)),
    pr.dlp_strata = c(0,0.66^(-2)),
    pr.dlp_stratXflo = c(0,0.5^(-2)),
    # hyper-prior on student-t distribution for p-random-effects (by time)
    pr.sig_p_t = c((0.05^(-2)),3),
    # beta prior on psi (recruitment, or 'removal-entry')
    pr.psi =array(c(3.5,0.4,0.444,0.5,0.571,0.667,3.5,3.6,3.556,3.5,3.429,3.333,3.5,0.4,0.444,0.5,0.571,0.667,3.5,3.6,3.556,3.5,3.429,3.333),c(T,2,2)) # slices are shape parameter
) # priors

jags.txt <- "
model{
  # recruitment process (technically, removal entry process)
  for(x in 1:n.sex){
    for(ti in 1:T){
      psi[ti,x] ~ dbeta(pr.psi[ti,x,1], pr.psi[ti,x,2])
    } # pp
  } # sex
  # all processes: mean, on logit
  llambda_mu ~ dnorm(pr.llambda_mu[1],pr.llambda_mu[2])
  lgte_mu ~ dnorm(pr.lgte_mu[1],pr.lgte_mu[2])
  lgaa_mu ~dnorm(pr.lgaa_mu[1],pr.lgaa_mu[2])
  lgbb_mu ~dnorm(pr.lgbb_mu[1],pr.lgbb_mu[2])
  lphi_mu ~dnorm(pr.lphi_mu[1],pr.lphi_mu[2])
  lp_mu ~dnorm(pr.lp_mu[1],pr.lp_mu[2])
  # all processes: hyper-prior: variance
  sig_p_all ~ dt(0,pr.sig_p_t[1], pr.sig_p_t[2]) T(0,1.5)
# priors on fixed effects
  # all processes: sex effects
  dlp_sex~dnorm(pr.dlp_sex[1],pr.dlp_sex[2])
  dlphi_sex~dnorm(pr.dlphi_sex[1],pr.dlphi_sex[2])
  # all processe: flooding effects [before]
  dlp_flo_dur ~ dnorm(pr.dlp_flo_dur[1],pr.dlp_flo_dur[2])
  # all processe: flooding effects [after]
  dlp_flo_aft ~ dnorm(pr.dlp_flo_aft[1],pr.dlp_flo_aft[2])
  # other strata effects on capture probability
  dlp_strata ~ dnorm(pr.dlp_strata[1],pr.dlp_strata[2])
  dlp_stratXflo ~ dnorm(pr.dlp_stratXflo[1],pr.dlp_stratXflo[2])
  # time effects on demographic processes 
  for(ti in 1:(T-1)){
      # reassemble on to the probability scale
      for(x in 1:n.sex){
        # model of assortment
        lambda[ti,x] <- 1/(1+exp(-llambda_mu ))
        # model of A->B
        gaa[ti,x] <- 1/(1+exp(-lgaa_mu )) # 
        # model of B->A
        gbb[ti,x] <- 1/(1+exp(-lgbb_mu )) # 
        # model of temporary emigration A->TE
        gate[ti,x] <- 1/(1+exp(-lgte_mu )) # -
        # model of temporary emigration B->TE
        gbte[ti,x] <- 1/(1+exp(-lgte_mu )) 
        # model on survival
        phi[ti,x] <- 1/(1+exp(-lphi_mu - dlphi_sex))
     } # time ti
  } # sex
  for(x in 1:n.sex){
     lambda[T,x] <- 1/(1+exp(-llambda_mu ))
  }
  # capture process
  for(k in 1:(E-1)){
     for(tti in 1:TTT.tot){
        dlp_t[tti,k] ~dnorm(0,pow(sig_p_all,-2))
        for(x in 1:n.sex){
          p[tti,k,x] <- PROP.AREA[tti,k]/(1+exp(-lp_mu - dlp_sex*equals(x,2) - dlp_strata*equals(k,2) - dlp_flo_dur*cov.p.flo[tti,k,1] - dlp_stratXflo*equals(k,2)*cov.p.flo[tti,k,1] - dlp_flo_aft*cov.p.flo[tti,k,2] - dlp_t[tti,k])) # hierarchical specification for all capture probabilities
        } # sex
     } # tii
  } # k
for(x in 1:n.sex){
  # first transition process: t0 to t1
  # states 1:A; 2:B; 3:TE; 4:non-recruitment; 5=dead
  tr[1,1,1,x] <- 1
  tr[2,1,1,x] <- 0
  tr[3,1,1,x] <- 0
  tr[4,1,1,x] <- 0
  tr[5,1,1,x] <- 0
  tr[6,1,1,x] <- 0
  tr[1,2,1,x] <- 0
  tr[2,2,1,x] <- 1
  tr[3,2,1,x] <- 0
  tr[4,2,1,x] <- 0
  tr[5,2,1,x] <- 0
  tr[6,2,1,x] <- 0
  tr[1,3,1,x] <- 0
  tr[2,3,1,x] <- 0
  tr[3,3,1,x] <- 1
  tr[4,3,1,x] <- 0
  tr[5,3,1,x] <- 0
  tr[6,3,1,x] <- 0
  tr[1,4,1,x] <- 0
  tr[2,4,1,x] <- 0
  tr[3,4,1,x] <- 1
  tr[4,4,1,x] <- 0
  tr[5,4,1,x] <- 0
  tr[6,4,1,x] <- 0
  tr[1,5,1,x] <- lambda[1,x]*psi[1,x] # recruit to A
  tr[2,5,1,x] <- (1-lambda[1,x])*psi[1,x] # recruit to B
  tr[3,5,1,x] <- 0 # TE latent states
  tr[4,5,1,x] <- 0 # TE latent states
  tr[5,5,1,x] <- 1-psi[1,x] # do not recruit
  tr[6,5,1,x] <- 0
  tr[1,6,1,x] <- 0
  tr[2,6,1,x] <- 0
  tr[3,6,1,x] <- 0
  tr[4,6,1,x] <- 0
  tr[5,6,1,x] <- 0
  tr[6,6,1,x] <- 0
  for(ti in 2:T){
    # the transition process: t to t+1
      tr[1,1,ti,x] <- (gaa[ti-1,x])*phi[ti-1,x]*(1-gate[ti-1,x]) # stays in A strata
      tr[2,1,ti,x] <- (1-gaa[ti-1,x])*phi[ti-1,x]*(1-gate[ti-1,x]) # goes to B strata
      tr[3,1,ti,x] <- gaa[ti-1,x]*gate[ti-1,x]*phi[ti-1,x] # goes to TE
      tr[4,1,ti,x] <- (1-gaa[ti-1,x])*gate[ti-1,x]*phi[ti-1,x] # goes to TE
      tr[5,1,ti,x] <- 0
      tr[6,1,ti,x] <- 1-phi[ti-1,x] # dies
      # from B strata to...
      tr[1,2,ti,x] <- (1-gbb[ti-1,x])*phi[ti-1,x]*(1-gbte[ti-1,x]) # goes to B strata
      tr[2,2,ti,x] <- (gbb[ti-1,x])*phi[ti-1,x]*(1-gbte[ti-1,x]) # stays in A strata
      tr[3,2,ti,x] <- (1-gbb[ti-1,x])*gbte[ti-1,x]*phi[ti-1,x] # A TE
      tr[4,2,ti,x] <- (gbb[ti-1,x])*gbte[ti-1,x]*phi[ti-1,x] # B TE
      tr[5,2,ti,x] <- 0
      tr[6,2,ti,x] <- 1-phi[ti-1,x] # dies
      # from TE-A to ...
      tr[1,3,ti,x] <- gaa[ti-1,x]*phi[ti-1,x]*(gate[ti-1,x]) # returns to A
      tr[2,3,ti,x] <- (1-gaa[ti-1,x])*phi[ti-1,x]*(gate[ti-1,x]) # return to B
      tr[3,3,ti,x] <- (1-gate[ti-1,x]) *phi[ti-1,x] # stays in TE
      tr[4,3,ti,x] <- 0
      tr[5,3,ti,x] <- 0
      tr[6,3,ti,x] <- 1-phi[ti-1,x] # dies
      # from TE-B to...
      tr[1,4,ti,x] <- (1-gbb[ti-1,x])*phi[ti-1,x]*(gbte[ti-1,x]) # back to A
      tr[2,4,ti,x] <- gbb[ti-1,x]*phi[ti-1,x]*(gbte[ti-1,x]) # back to B
      tr[3,4,ti,x] <- 0
      tr[4,4,ti,x] <- (1-gbte[ti-1,x])*phi[ti-1,x] # stays in TE-B
      tr[5,4,ti,x] <- 0
      tr[6,4,ti,x] <- 1-phi[ti-1,x] # dies
      # from recruitment to...
      tr[1,5,ti,x] <- lambda[ti,x]*psi[ti,x] # recruit to A
      tr[2,5,ti,x] <- (1-lambda[ti,x])*psi[ti,x] # recruit to B
      tr[3,5,ti,x] <- 0
      tr[4,5,ti,x] <- 0
      tr[5,5,ti,x] <- 1-psi[ti,x] # do not recruit
      tr[6,5,ti,x] <- 0
      tr[1,6,ti,x] <- 0
      tr[2,6,ti,x] <- 0
      tr[3,6,ti,x] <- 0
      tr[4,6,ti,x] <- 0 # 
      tr[5,6,ti,x] <- 0 # 
      tr[6,6,ti,x] <- 1 # stays dead
  } # done transition process
  # the emission process (capture process)
  for(tti in 1:TTT.tot){
      em[1,1,tti,x] <- p[tti,1,x]
      em[2,1,tti,x] <- 0
      em[3,1,tti,x] <- 1-p[tti,1,x]
      em[1,2,tti,x] <- 0
      em[2,2,tti,x] <- p[tti,2,x]
      em[3,2,tti,x] <- 1-p[tti,2,x]
      em[1,3,tti,x] <- 0
      em[2,3,tti,x] <- 0
      em[3,3,tti,x] <- 1 # te, unrecruited
      em[1,4,tti,x] <- 0
      em[2,4,tti,x] <- 0
      em[3,4,tti,x] <- 1 # te, unobservable
      em[1,5,tti,x] <- 0
      em[2,5,tti,x] <- 0
      em[3,5,tti,x] <- 1 # unrecruited, unobservable
      em[1,6,tti,x] <- 0
      em[2,6,tti,x] <- 0
      em[3,6,tti,x] <- 1 # dead, unobservable
  }
} # sex 
for(i in 1:M){
  # first primary period
  z[i,1] ~ dcat(tr[1:K,5,1,x.ind[i]]) # first primary period
  for(s in 1:T2[1]){ 
     y[i,s,1] ~ dcat(em[1:E, z[i,1], tt.ind[s,1], x.ind[i]])
  }
  # 2-T captures occassion
  for(ti in 2:T){
     z[i,ti] ~ dcat(tr[1:K,z[i,ti-1],ti,x.ind[i]]) # markov
     for(s in 1:T2[ti]){ 
        y[i,s,ti] ~ dcat(em[1:E, z[i,ti], tt.ind[s,ti], x.ind[i]])
     } # s
  } # ti primary periods
} # individuals M
# loglike: for WAIC estimation (conditional loglikelkhood, unfortunately)
for(i in 1:N.obs){
  for(ti in 1:T){
     for(s in 1:T2[ti]){ 
        llits[i,s,ti] <- logdensity.cat(y[i,s,ti], em[1:E, z[i,ti], tt.ind[s,ti], x.ind[i]])
     } # s
     llit[i,ti] <- sum(llits[i,1:T2[ti],ti])
  } # ti primary periods
  lli[i] <- sum(llit[i,1:T])
} # individuals N.obs
# POPULATION ABUNDANCE, RECRUITS, SUPER-POP
for(i in 1:M){
  for(x_ in 1:n.sex){
    alive_i[i,1,x_] <- step(4-z[i,1])*equals(x.ind[i],x_) # is animal alive or not  
    recruit_i[i,1,x_] <- alive_i[i,1,x_]*equals(x.ind[i],x_) # new recruit or
    for(k in 1:n.strata){
      N.onsite_i[i,1,k,x_]<-equals(z[i,1],k)*equals(x.ind[i],x_) # is on-site or not
    }
    for(ti in 2:T){
      alive_i[i,ti,x_] <- step(4-z[i,ti])*equals(x.ind[i],x_) # is animal alive or not
      recruit_i[i,ti,x_] <- alive_i[i,ti,x_]*equals(z[i,ti-1],5)*equals(x.ind[i],x_) # new recruit? 
      for(k in 1:n.strata){
        N.onsite_i[i,ti,k,x_]<-equals(z[i,ti],k)*equals(x.ind[i],x_) # is on site or not
      }
    }
  }
}
for(t_ in 1:T){
  for(x_ in 1:n.sex){ # loop through each sex
     for(k in 1:n.strata){
        N.onsite[t_,x_,k] <- sum(N.onsite_i[1:M,t_,k,x_])
     } # k
     alive[t_,x_] <- sum(alive_i[1:M,t_,x_])
     recruits[t_,x_] <- sum(recruit_i[1:M,t_,x_])
  } # x_
} # t_ 
}"
modname <- "JAGS_hierarchical_mcrd_cagnazzi.JAG"
sink(file=modname) # 
cat(jags.txt,fill=TRUE) 
sink() # this cloes the connection to file 

# DATA AUGMENTATION
n.aug.x <- round(0.85*table(sex.vector)) # females and males
n.aug <- sum(n.aug.x)
m <- n.aug+N.obs
sex.aug.v <- c(sex.vector,unlist(mapply(1:n.sex,n.aug.x, FUN=rep))) # basically double the amount of pseudo females and pseudo males
#table(sex.aug.v)
y.aug <- array(NA, c(m,dim(y)[2],dim(y)[3]))
y.aug[1:N.obs,,]<-y # 
y.aug[(N.obs+1):m,,]<-y[1:n.aug,,]*0+3 # all zero capture histories (3 is no capture)
x.ind <- sex.aug.v

# assemble all the data together
jags.data <- c(list(
    y=y.aug,
    M=nrow(y.aug),
    N.obs=N.obs,
    n.sex=n.sex,
    x.ind = x.ind, # sex of individuals
    T=T, # primary periods
    tt.ind=tt.ind, # map secondary periods to primary periods
    T2=T2,
    n.strata =2,
    TTT.tot=TTT.tot, 
    K=6, # number of latent states 
    E=3, # number of 'event's (seen in strata1, strat2 and no-capture)
    PROP.AREA=PROP.AREA, # proprotion of area surveyed              
    cov.p.flo=cov.p.flooding), # flooding covariate
    priors) 


# mcmc PARAMETERS
n.chains = 3
n.adapt <- 4000
n.burn <- 20000
n.samp <- 4000
n.iter <- 200000
thin_ = round(n.iter/n.samp)

# initialize the random variables (make a function for JAGS to call)
init.f <- lexiscope(y.aug,T=T,te1=3,te2=4,re=5,de=6,out_=3,n.sex=n.sex)
inits <- init.f()

# SINGLE CHAIN    
m<-jags.model(file = modname, data=jags.data, init=inits, n.chains=1,n.adapt=n.adapt)
# BURN-IN
update(m,n.burn)
post <- coda.samples(m,variable.names=c("lli","phi","psi", "gate","gbte","p","lambda","gaa","gbb","dlphi_sex", "dlp_sex", "dlp_flo_dur","dlp_flo_aft","dlp_strata", "dlp_stratXflo", "sig_p_all","N.onsite","recruits","alive"),n.iter=n.iter,thin=thin_)

# PARALLEL CHAINS
b2 <- jags.parallel(data = jags.data, inits=init.f, parameters.to.save=c("lli","phi","psi", "gate","gbte","p","lambda","gaa","gbb","dlphi_sex", "dlp_sex", "dlp_flo_dur","dlp_flo_aft","dlp_strata", "dlp_stratXflo", "sig_p_all","N.onsite","recruits","alive"),model.file=modname, n.chains=n.chains,n.iter=n.burn+round(n.iter/n.chains),n.burnin=n.burn,n.thin=thin_,export_obj_names = c("thin_","n.iter","n.chains","jags.data","n.burn","init.f","y.aug"))
post <- as.mcmc(b2)
