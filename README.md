# Case studies for "How to build your latent Markov model - the role of time and space"

This repository contains R code for the case studies presented in the paper "How to build your latent Markov model - the role of time and space". The paper provides a unified view on hidden Markov models (HMMs), state-space models (SSMs) and Markov modulated (marked) Poisson processes (MM(M)PPs), discussing model choices regarding the operationalization of time and the state space. The code provided here illustrates the concepts discussed in the paper by fitting the models to real-world data.

The code is organized as follows: Each case study presented in the paper has its own `.R` file. The files are:

* `HMM_Application.R`: Illustrating discrete-time HMMs by analysing the movement track of an elephant
from Etosha National Park.
* `ctHMM_Application.R`: Illustrating continuous-time HMMS by analysing data collected on forced expiratory volumes (FEV)
measured for lung transplant recipients.
* `SSM_Application.R`: Illustrating SSMs by fitting a stochastic volatility model to daily Bitcoin returns.
* `ctSSM_Application.R`: Illustrating continuous-time SSMs by analysing seven-metre throws in Handball to find evidence for a hot hand effect.
* `MMMPP_Application.R`: Illustrating MMPPs by analysing the surfacing times of minke whales.

Additionally, the file `likelihood_functions.R` contains all negative log-likelihood functions used to fit the above models via direct numerical maximum likelihood estimation as explained in the paper. For each fitted model, we provide a base R version as well as a high performance version of the likelihood by using the R package `LaMa` developed for fast statistical inference in all model classes discussed.

The figures presented in the paper are saved in the folder `figs` and all data necessary will either be downloaded automatically by running the case study code, or is included in the folder `data` and will be loaded automatically.
