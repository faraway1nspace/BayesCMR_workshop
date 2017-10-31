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

A brief introduction to Bayesianism, in general. PDF presentation [here](http://foo.bar). No R/JAGS code 


PART 2: Introduction to Priors and Probabilities
------------------------------------------------

This was a simple exercise for participants to study some common probability distributions available in [JAGS](http://www.google.com/url?q=http://mcmc-jags.sourceforge.net/&sa=U&ved=0ahUKEwjf38Gb6tXWAhUIVLwKHejnA2EQFggdMAQ&usg=AOvVaw3VPi0Ffru14OG--3erpJZh). A PDF lecture is here. A few helpful R functions are available here. 


PART 3: Introduction to JAGS
----------------------------
This part included 3 exercises to familiarize themselves with the JAGS syntax and R workflow:
- a simple Bernoulli model for average annual survival of 30 dolphins (with Beta Priors)
- a Bernoulli model for average annual survival of 30 dolphins (with logit-Normal priors)
- a logistic-regression of average annual survival of 30 dolphins (with logit-Normal priors)

See the exercise notes at the end of the PDF from Part 1

PART 4: Intro to Hidden Markov Models, a Unifying Framework for CMR
-------------------------------------------------------------------




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
