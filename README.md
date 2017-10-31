# Bayesian Capture-Recapture in JAGS: Workshop for SMM 2017 Biennial Conference
Welcome to the SMM2017 workshop on Bayesian Capture-recapture. This page contains materials from the full-day workshop, including
- 8 lectures and tutorials
- R and JAGS code

The focus of the workshop was to introduce [JAGS](http://www.google.com/url?q=http://mcmc-jags.sourceforge.net/&sa=U&ved=0ahUKEwjf38Gb6tXWAhUIVLwKHejnA2EQFggdMAQ&usg=AOvVaw3VPi0Ffru14OG--3erpJZh) for Bayesian Capture-Mark-Recapture (including models CJS, POPAN, PCRD, MSCRD as a more general type of Hidden Markov Model). The underlying idea is based on the paper of [Rankin and Nicholson et al. 2016](http://journal.frontiersin.org/article/10.3389/fmars.2016.00025) and [Rankin's 2017 Ph.D. thesis](http://researchrepository.murdoch.edu.au/id/eprint/38257/). Researchers may find the lectures and R/JAGS code useful to make their own CJS, POPAN, PCRD, and MSCRD models in JAGS.

Feel free to submit feature requests and report issues in the Github issue tracker page (see link above)

OUTLINE
-------

![Outline](/img/outline.png)

PART 1: History, Philosophy and Properties of Bayesian Inference
----------------------------------------------------------------

A brief introduction to Bayesianism. What is *posterior inference*? What divides *subjective* vs objective Bayesians? How do sample-size and *prior information* influence the frequency properties of a Bayes estimator?  What is sample-based approximation? [PDF presentation here](https://github.com/faraway1nspace/BayesCMR_workshop/blob/master/PART1_introBayes/bayesian_intro.pdf). No R/JAGS code 


PART 2: Introduction to Priors and Probabilities
------------------------------------------------

This was a simple exercise for participants to study some common probability distributions available in [JAGS](http://www.google.com/url?q=http://mcmc-jags.sourceforge.net/&sa=U&ved=0ahUKEwjf38Gb6tXWAhUIVLwKHejnA2EQFggdMAQ&usg=AOvVaw3VPi0Ffru14OG--3erpJZh). A PDF lecture is here. A few helpful R functions are available here. 


PART 3: Introduction to JAGS
----------------------------
This part included 3 exercises to familiarize participants with the JAGS syntax and R workflow. The emphasis in on i) encoding prior beliefs and ii) encoding the likelihood (joint probability distribution of the data), together which serve as the basic skeleton for all subsequent JAGS models. There are three JAGS exercises as part of this section:
- a simple Bernoulli model for average annual survival of 30 dolphins (with Beta Priors)
- a Bernoulli model for average annual survival of 30 dolphins (with logit-Normal priors)
- a logistic-regression of average annual survival of 30 dolphins (with logit-Normal priors)

See the notes at the end of the PDF from Part 1. 

PART 4: Intro to Hidden Markov Models, a Unifying Framework for CMR
-------------------------------------------------------------------
This tutorial shows how the CJS, POPAN, PCRD, and MSCRD are just variants of a more general Hidden-Markov Model. Participants learn about how to specif
- 

how to run a simple HMM (like a POPAN model) in JAGS for a single observation, then work to generalize the 



INSTALLATION
------------
Participants should first:
- install [JAGS](http://www.google.com/url?q=http://mcmc-jags.sourceforge.net/&sa=U&ved=0ahUKEwjf38Gb6tXWAhUIVLwKHejnA2EQFggdMAQ&usg=AOvVaw3VPi0Ffru14OG--3erpJZh)
- then the R package rjags > install.packages("rjags").
- check that your installation was successful. Please navigate to the *test_script* directory and run the file *R_test_JAGS_script.R* (NOTE: windows users please use the WINDOWS version). If successful, the script should make a plot with title "Success!"


The workshop is based heavily on the publication: <b>[Rankin RW, Nicholson K, Allen SJ, Krützen M, Bejder L, Pollock KH. 2016. A full-capture Hierarchical Bayesian model of Pollock's Closed Robust Design and application to dolphins. Frontiers in Marine Science 3:25, doi: 10.3389/fmars.2016.00025](http://journal.frontiersin.org/article/10.3389/fmars.2016.00025)</b>.

Sincerely, \\
Robert W. Rankin, Dr. \\
Krista E Nicholson, Ph.D Candidate
