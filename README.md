# Bayesian Capture-Recapture in JAGS: Workshop for SMM 2017 Biennial Conference
Welcome to the SMM2017 workshop on Bayesian Capture-recapture. This page contains materials from the full-day workshop, including
- 8 lectures and tutorials
- R and JAGS code

The focus of the workshop was to introduce [JAGS](http://mcmc-jags.sourceforge.net/) for Bayesian Capture-Mark-Recapture (including models CJS, POPAN, PCRD, MSCRD as a more general type of Hidden Markov Model). The underlying idea is based on the paper of [Rankin and Nicholson et al. 2016](http://journal.frontiersin.org/article/10.3389/fmars.2016.00025) and [Rankin's 2017 Ph.D. thesis](http://researchrepository.murdoch.edu.au/id/eprint/38257/). Researchers may find the lectures and R/JAGS code useful to make their own CJS, POPAN, PCRD, and MSCRD models in JAGS.

Feel free to submit feature requests and report issues in the Github issue tracker page (see link above)

**Data** A special thanks to those researchers who let us their odontocetes capture-recapture data-set for some of the practical exercises, including Dr. Tim Hunt (*Humpback dolphins*, POPAN, Part 6), Krista Nicholson (*Bottlenose dolphins*, PCRD, Part 7), and Dr. Daniele Cagnazzi (*Humpback dolphins*, MSCRD, Part 8, unpublished). These data-sets are detailed in their respective publications: [Hunt et al](http://www.int-res.com/abstracts/esr/v32/p71-88/), [Nicholson et al](http://dx.doi.org/10.1071/MF12210), and unpublished data from Cagnazzi (but see his [Doctoral Thesis](http://epubs.scu.edu.au/theses/344/)). 

OUTLINE
-------

![Outline](/img/outline.png)

PART 1: Bayesian Inference: History, Philosophy, and Properties
--------------------------------------------------------------

A brief introduction to Bayesianism. What is *posterior inference*? What divides *subjective* vs objective Bayesians? How do sample-size and *prior information* influence the frequency properties of a Bayes estimator?  What is sample-based approximation? [PDF presentation here](https://github.com/faraway1nspace/BayesCMR_workshop/blob/master/PART1_introBayes/bayesian_intro.pdf). No R/JAGS code 


PART 2: Priors and Probabilities
------------------------------------------------

This was a simple exercise for participants to study some common probability distributions available in [JAGS](http://mcmc-jags.sourceforge.net/). A PDF lecture is here. A few helpful R functions are available here. 


PART 3: JAGS
----------------------------
This part included 3 exercises to familiarize participants with the JAGS syntax and R workflow. The emphasis in on i) encoding prior beliefs and ii) encoding the likelihood (joint probability distribution of the data), together which serve as the basic skeleton for all subsequent JAGS models. There are three JAGS exercises in the [R file](./PART3_introJAGS/) as part of this section:
- a simple Bernoulli model for average annual survival of 30 dolphins (with *Beta* Priors)
- a Bernoulli model for average annual survival of 30 dolphins (with *logit-Normal* priors)
- a logistic-regression of average annual survival of 30 dolphins (with logit-Normal priors)

See the notes at the end of the PDF from Part 1. 

PART 4: Hidden Markov Models, a Unifying Framework for CMR
-------------------------------------------------------------------

This tutorial shows how the CJS, POPAN, PCRD, and MSCRD are just variants of a more general Hidden-Markov Model. Participants learn about how to specify transmission matrices and emission matrices in the JAGS syntax as well as the rationale behind the latent state markov process and the complete data likelihood to recast a variety of CMR models as HMMs. 
1. exercise 1 is to run a simple HMM (like a POPAN model) for a single capture-history, then
2. exercise 2 generalizes the HMM for multiple capture histories and time-varying parameters, then
3. exercise 3 add's POPAN-like derivatives to the JAGS script (Super-population, population abundance births, and probability of entry)

See the [lecture PDF](https://github.com/faraway1nspace/BayesCMR_workshop/blob/master/PART4_introHMM/hmm_intro.pdf) and the [PDF describing the tutorials](https://github.com/faraway1nspace/BayesCMR_workshop/blob/master/PART4_introHMM/hmm_practical.pdf). The R code is in [PART 4 directory](./PART4_introHMM/).

**Key concepts and techniques**:
- HMM transmission and emission matrices
- latent states and the latent Markov process
- complete data likelihood

PART 5: Cormack-Jolly-Seber as a HMM
------------------------------------

This tutorial introduces the idea of *conditioning on the first capture* and re-casts the CJS model as a very simple HMM. Capture recapture data is from [Nicholson et al](http://dx.doi.org/10.1071/MF12210). Also discussed are the Horvitz-Thompson estimator of abundance and linear-model-type specifications of capture-histories (e.g., for individual heterogeneity).

See the [lecture PDF](./PART5_CJS/cjs.pdf) and [R/JAGS tutorial files](./PART5_CJS/).

**Key concepts and techniques**:
- first-capture vs full-capture models
- Horvitz-Thompson abundance estimation
- linear models for CMR parameters


PART 7: POPAN and Model-Selection
---------------------------------
foobar 

PART 8: Pollock's Closed Robust Design, and Hierarchical Bayes
--------------------------------------------------------------
This tutorial introduces the Kendall model for temporary emigration and the Pollock's Closed Robust Design. See the [lecture PDF](./PART7_PCRD/pcrd.pdf) and [R/JAGS tutorial files](./PART7_PCRD/) called `R_pcrd.R`.

**Key concepts and techniques**:
- secondary periods and population-closure
- nested for loops
- Hierarchical Bayes and random-effects
- hyperpriors and shrinkage
- 


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
