# This case study illustrates MMPPs by analyzing the surfacing times of minke whales.

# Loading packages --------------------------------------------------------

# install.packages("LaMa")
library(LaMa)
# install.packages("expm")
library(expm)


# Loading the minke whale data --------------------------------------------

surftimes = readRDS("./data/minkewhales.rds")


# Some EDA ----------------------------------------------------------------

pdf("./figs/mmpp_divetimes.pdf", width = 7, height = 4)
par(mfrow = c(4, 1), mar = c(4.1, 2, 0, 2) + 0.1)
for(i in 1:4){
  plot(surftimes[[i]], rep(1, length(surftimes[[i]])), ylim = c(0, 1.5), type = "h", 
       yaxt = "n", xaxt = "n", xlab = "Dive time in minutes", ylab = "", 
       bty = "n", xlim = c(0, 3500))
  axis(1, at = seq(0, 3500, 500), labels = seq(0, 3500, 500))
}
dev.off()


# Loading the negative log-likelihood function ----------------------------

source("likelihood_functions.R")


# Model fitting -----------------------------------------------------------

lambda0 = 1 / c(quantile(unlist(surftimes), 0.2), quantile(unlist(surftimes), 0.8))
qs0 = rep(mean(lambda0) / 5, 2)

# initial values
theta.star = log(c(lambda0, qs0))

mods = list() # separate models for all animals
for(i in 1:length(surftimes)){
  mods[[i]] = nlm(mllk_mmpp_fast, theta.star, timediff = diff(surftimes[[i]]), N = 2,
          print.level = 1, iterlim = 1000, stepmax = 10, hessian = TRUE)
}

Lambdas = matrix(NA, nrow = length(mods), ncol = 2)
for(i in 1:length(mods)){
  Lambdas[i, ] = exp(mods[[i]]$estimate[1:2])
}
Qs = matrix(NA, nrow = length(mods), ncol = 2)
for(i in 1:length(mods)){
  Qs[i, ] = exp(mods[[i]]$estimate[4:3])
}

# Rates for each animal
Lambdas
Qs

# average waiting times for each animal
round(1 / Lambdas, 2)
# average sojourn time for each animal
round(1 / Qs, 2)

