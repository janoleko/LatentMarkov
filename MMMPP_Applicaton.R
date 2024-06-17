
# Loading packages --------------------------------------------------------

# install.packages("LaMa")
library(LaMa)
# install.packages("expm")
library(expm)


# Loading the minke whale data --------------------------------------------

mink = scan("./data/divepatterns.txt")

splitInd = which(mink == 1) # split into list by splitInd
dive = split(mink, cumsum(mink == 1)) # list of 92 individual tracks

# filter for more than 60 observations
surftimes = dive[which(unlist(lapply(dive, length)) >= 60)]
# 7 track have more than 60 observations

# we only consider whales 2, 3, 4 and 5 for the application, as the other whales
# lead to the special case of a homogeneous Poisson process -> absorbing state
surftimes = surftimes[2:5]


# Some EDA ----------------------------------------------------------------

pdf("./figs/mmpp_divetimes.pdf", width = 7, height = 4)
par(mfrow = c(4,1), mar = c(4.1,2,0,2) + 0.1)
for(i in 1:4){
  plot(surftimes[[i]], rep(1, length(surftimes[[i]])), ylim = c(0,1.5), type = "h", 
       yaxt = "n", xaxt = "n", xlab = "Dive time in minutes", ylab = "", bty = "n", xlim = c(0, 3500))
  axis(1, at = seq(0, 3500, 500), labels = seq(0, 3500, 500))
}
dev.off()


# Loading the negative log-likelihood function ----------------------------

source("likelihood_functions.R")

### longer version in base R
# mllk_mmpp_slow = function(theta.star, timediff, N=2){
#   lambda = exp(theta.star[1:N]) # state specific rates
#   Q = diag(N) # generator matrix
#   Q[!Q] = exp(theta.star[N+1:(N*(N-1))])
#   diag(Q) = 0
#   diag(Q) = -rowSums(Q)
#   delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
#   # we split the Omega matrix into (Q-Lambda)*dt and Lambda
#   Qube = array(dim = c(N,N,length(timediff)))
#   for(t in 1:length(timediff)){
#     Qube[,,t] = expm::expm((Q-diag(lambda))*timediff[t])
#   }
#   allprobs = matrix(lambda, nrow = length(timediff)+1, ncol = N, byrow = T) # Lambda
#   allprobs[1,] = 1 # no Lambda matrix at the start of the likelihood
#   # forward algorithm in R
#   foo = delta%*%diag(allprobs[1,])
#   l = log(sum(foo))
#   phi = foo / sum(foo)
#   for(t in 1:length(timediff)){
#     foo = phi %*% Qube[,,t] %*% diag(allprobs[t+1,])
#     l = l + log(sum(foo))
#     phi = foo / sum(foo)
#   }
#   -l
# }

### short version using the package LaMa
# mllk_mmpp_fast = function(theta.star, timediff, N=2){
#   lambda = exp(theta.star[1:N]) # state specific rates
#   Q = diag(N) # generator matrix
#   Q[!Q] = exp(theta.star[N+1:(N*(N-1))])
#   diag(Q) = 0
#   diag(Q) = -rowSums(Q)
#   delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary initial distribution
#   # we split the Omega matrix into (Q-Lambda)*dt and Lambda for efficiency
#   Qube = LaMa::tpm_cont(Q-diag(lambda), timediff) # (Q-Lambda)*dt
#   allprobs = matrix(lambda, nrow = length(timediff)+1, ncol = N, byrow = T) # Lambda
#   allprobs[1,] = 1 # no Lambda matrix at the start of the likelihood
#   # forward algorithm in C++ using LaMa
#   -LaMa::forward_g(delta, Qube, allprobs)
# }


# Model fitting -----------------------------------------------------------

lambda0 = 1 / c(quantile(unlist(surftimes), 0.4), quantile(unlist(surftimes), 0.6))
qs0 = rep(mean(lambda0)/5, 2)

# initial values
theta.star = log(c(lambda0, qs0))

mods = list() # separate models for all animals
for(i in 1:length(surftimes)){
  mods[[i]] = nlm(mllk_mmpp_fast, theta.star, timediff = diff(surftimes[[i]]), N = 2,
          print.level = 2, iterlim = 1000, stepmax = 10, hessian = TRUE)
}

Lambdas = matrix(NA, nrow = length(mods), ncol = 2)
for(i in 1:length(mods)){
  Lambdas[i,] = exp(mods[[i]]$estimate[1:2])
}
Qs = matrix(NA, nrow = length(mods), ncol = 2)
for(i in 1:length(mods)){
  Qs[i,] = exp(mods[[i]]$estimate[4:3])
}

# Rates for each animal
Lambdas
Qs

# average waiting times for each animal
round(1 / Lambdas, 2)
# average sojoun time for each animal
round(1 / Qs, 2)

