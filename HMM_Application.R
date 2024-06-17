# This case study illustrates discrete-time HMMs by analysing the movement track of an elephant from Etosha National Park.

# Loading packages --------------------------------------------------------

# install.packages("LaMa")
library(LaMa)
# install.packages("CircStats")
library(CircStats)
# install.packages("moveHMM")


# Loading the data --------------------------------------------------------

elephant_data = read.csv("http://www.rolandlangrock.com/elephant_data.csv")

# Originally, 15 elephants from the Etosha National Park where fitted with GPS 
# collars to collect hourly movement data (@tsalyuk2019temporal). 
# For simplicity, the raw data have already been preprocessed and we consider 
# the track of only one individual.


# Some preprocessing ------------------------------------------------------

## Defining time of day variable (could also be done with lubridate)
elephant_data$timestamp = strptime(elephant_data$timestamp, 
                                   "%Y-%m-%d %H:%M:%S", tz = "GMT")
hours = as.numeric(format(elephant_data$timestamp, "%H"))
minutes = as.numeric(format(elephant_data$timestamp, "%M"))
elephant_data$timeOfDay = hours + (minutes/60)

## calculating step lengths and turning angles
data = moveHMM::prepData(trackData = elephant_data, 
                         coordNames = c("location.long", "location.lat"), 
                         type = "LL")
data$timeOfDay2 = round(data$timeOfDay)

# we only have 2 step lengths that are exactly zero. 
sum(na.omit(data$step) == 0)
# This does not really justify using a zero inflated step-length distribution
# and we just replace these zeros with very small values
data$step[which(data$step == 0)] = 1e-5


# EDA ---------------------------------------------------------------------

nrow(data)
summary(data$step)
data$timestamp[1]; data$timestamp[nrow(data)]


# Loading the negative log-likelihood function ----------------------------

source("likelihood_functions.R")


# Fitting a homogeneous HMM -----------------------------------------------

# initial parameters for numerical optimisation
stepMean0 = c(0.1, 0.4, 1)
stepSd0 = c(0.1, 0.4, 1)
angleMean0 = rep(0,3)
angleConcentration0 = c(0.1, 0.9, 2)

# working scale initial parameter vector
theta.star0 = c(log(stepMean0), log(stepSd0), # steppars
                angleMean0, log(angleConcentration0), # anglepars
                rep(-2, 6)) # tpmpars

# fitting the model by numerically optimising the log-likelihood with nlm
mod_elephant1 = nlm(mllk_fast, theta.star0, X = data, N = 3, 
                    print.level = 2, hessian = TRUE)

# obtaining the estimated parameters
## parameters for state-dependent distributions
N = 3
(mu = exp(mod_elephant1$estimate[1:N]))
(sigma = exp(mod_elephant1$estimate[N+1:N]))
(mu.turn = mod_elephant1$estimate[2*N+1:N])
(kappa = exp(mod_elephant1$estimate[3*N+1:N]))
## transition probability matrix
Gamma = tpm(mod_elephant1$estimate[4*N+1:(N*(N-1))])
round(Gamma, 3)
## stationary distribution (proportion of time spent in each state)
(delta = LaMa::stationary(Gamma)) # there is also a stationary() function in moveHMM, but here we need LaMa

# visualising the fitted state-dependent and marginal distributions
color = c("orange", "deepskyblue", "seagreen2") # color vector

# histogram and the marginal distribution
par(mfrow = c(1,2))
## step length
hist(data$step, prob = TRUE, border = "white", xlim = c(0,4), breaks = 50,
     main = "", xlab = "step length", ylab = "density", cex = 0.9)
shape = mu^2/sigma^2; scale = sigma^2/mu
# component distributions
for(j in 1:N){
  curve(delta[j]*dgamma(x, shape=shape[j], scale=scale[j]), 
        add = T, col = color[j], lwd = 2, n = 300)
}
# marginal as sum
curve(delta[1]*dgamma(x, shape=shape[1], scale=scale[1])+
        delta[2]*dgamma(x, shape=shape[2], scale=scale[2])+
        delta[3]*dgamma(x, shape=shape[3], scale=scale[3]), 
      add = T, lty = 2, lwd = 2, n = 300)
legend("topright", col = color, lwd = 2, legend = paste("state", 1:3), bty = "n")

## turning angle
hist(data$angle, prob = TRUE, border = "white", 
     main = "", xlab = "turning angle", ylab = "density", cex = 0.9)
# component distributions
for(j in 1:N){
  curve(delta[j]*CircStats::dvm(x, mu.turn[j], kappa[j]), 
        add = T, col = color[j], lwd = 2, n = 300)
}
# marginal as sum
curve(delta[1]*CircStats::dvm(x, mu.turn[1], kappa[1])+
        delta[2]*CircStats::dvm(x, mu.turn[2], kappa[2])+
        delta[3]*CircStats::dvm(x, mu.turn[3], kappa[3]), 
      add = T, lty = 2, lwd = 2, n = 300)
legend("topright", col = color, lwd = 2, legend = paste("state", 1:3), bty = "n")


# state decoding
# for general state decoding we can also use LaMa. We must however first calculate the allprobs matrix
allprobs = matrix(1, nrow(data), N)
ind = which(!is.na(data$step) & !is.na(data$angle))
for(j in 1:N){
  allprobs[ind,j] = dgamma(data$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
    CircStats::dvm(data$angle[ind], mu.turn[j], kappa[j])
}
mod_elephant1$states = LaMa::viterbi(delta, Gamma, allprobs) # in moveHMM there is also a viterbi() function. Avoid confusion

# install.packages("scales")
library(scales) # for color transparency
par(mfrow = c(1,1))
plot(data$x, data$y, type = "l", xlab = "x", ylab = "y")
points(data$x, data$y, pch = 20, col = alpha(color[mod_elephant1$states], 0.4))

par(mfrow = c(1,1))
plot(data$timestamp, data$step, col = color[mod_elephant1$states], bty = "n", 
     ylab = "step length", xlab = "time")
plot(data$timestamp, data$angle, col = color[mod_elephant1$states], bty = "n", 
     ylab = "turning angle", xlab = "time")


# Covariate effect: Time of day -------------------------------------------

# We suspect the animal's behaviour varies with the time of day. 
# Thus we include the time of day as a covariate in our model.


# initial parameters for numerical optimisation
stepMean0 = c(0.1, 0.4, 1)
stepSd0 = c(0.1, 0.4, 1)
angleMean0 = rep(0,3)
angleConcentration0 = c(0.1, 0.9, 2)

# working scale initial parameter vector
set.seed(123)
theta.star0 = c(log(stepMean0), log(stepSd0), # steppars
                angleMean0, log(angleConcentration0), # anglepars
                rep(-2, 6), runif(6*2, -0.1, 0.1)) # tpmpars

# fitting the model by numerically optimising the log-likelihood with nlm
mod_elephant2 = nlm(mllk_HMM_fast, theta.star0, X = data, N = 3, 
                    print.level = 2, iterlim = 1000, hessian = TRUE)

# getting state process parameters
N = 3
beta = matrix(mod_elephant2$estimate[4*N+1:(3*N*(N-1))], ncol = 3)

# Visualising periodically stationary distribution
len = 300
todseq = seq(0, 24, length=len)
Delta = matrix(NA, len, N)
# here we interpolate to get a smooth version of the periodically stationary distribution
for(t in 1:len){
  Gamma = tpm_p(tod = (todseq[t]+1:24) %% 24, beta = beta)
  Delta[t,] = stationary_p(Gamma, t = 1)
}

# Obtaining confidence bands via Monte Carlo
# drawing from the approximate normal distribution of the MLE
thetas = mvtnorm::rmvnorm(1000, mean=mod_elephant2$estimate, 
                          sigma=solve(mod_elephant2$hessian))
Deltas = array(dim = c(len, N, 1000))
for(i in 1:1000){
  # cat("\n", i)
  for(t in 1:len){
    beta = matrix(thetas[i,4*N+1:(3*N*(N-1))], ncol = 3)
    Gamma = tpm_p(tod = (todseq[t]+1:24) %% 24, beta = beta)
    Deltas[t,,i] = stationary_p(Gamma, t = 1)
  }
}
DeltaCI = array(dim = c(len, N, 2))
for(j in 1:N){
  for(t in 1:len){
    DeltaCI[t,j,] = quantile(Deltas[t,j,], probs = c(0.025, 0.975))
  }
}


# Full Figure -------------------------------------------------------------
# as shown in the paper

pdf("./figs/hmm_elephant.pdf", width = 8, height = 3)

m = matrix(c(1,1,1,2,3,4), nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(0.13, 1, 1))
par(mar = c(0,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", inset = 0, col = color, 
       legend = c("resting", "foraging/ searching", "travelling"), lwd = 2,
       horiz = T, bty = "n", text.width = c(0, 0.15, 0))

par(mar = c(5.5,4,2,1), xpd = F)

# histogram and the marginal distribution
## step length
hist(data$step, prob = TRUE, border = "white", xlim = c(0,4), breaks = 50,
     main = "", xlab = "step length", ylab = "density", cex = 0.9)
shape = mu^2/sigma^2; scale = sigma^2/mu
# component distributions
for(j in 1:N){
  curve(delta[j]*dgamma(x, shape=shape[j], scale=scale[j]), 
        add = T, col = color[j], lwd = 2, n = 300)
}
# marginal as sum
curve(delta[1]*dgamma(x, shape=shape[1], scale=scale[1])+
        delta[2]*dgamma(x, shape=shape[2], scale=scale[2])+
        delta[3]*dgamma(x, shape=shape[3], scale=scale[3]), 
      add = T, lty = 2, lwd = 2, n = 300)

## turning angle
hist(data$angle, prob = TRUE, border = "white", 
     main = "", xlab = "turning angle", ylab = "density", cex = 0.9)
# component distributions
for(j in 1:N){
  curve(delta[j]*CircStats::dvm(x, mu.turn[j], kappa[j]), 
        add = T, col = color[j], lwd = 2, n = 300)
}
# marginal as sum
curve(delta[1]*CircStats::dvm(x, mu.turn[1], kappa[1])+
        delta[2]*CircStats::dvm(x, mu.turn[2], kappa[2])+
        delta[3]*CircStats::dvm(x, mu.turn[3], kappa[3]), 
      add = T, lty = 2, lwd = 2, n = 300)

# periodically stationary distribution
plot(NA, xlim = c(0,24), ylim = c(0,1), bty = "n", xaxt = "n",
     xlab = "time of day", ylab = "state occupancy probabilities")
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
for(j in 1:N){ 
  polygon(c(todseq, rev(todseq)), c(DeltaCI[,j,1], rev(DeltaCI[,j,2])), col = alpha(color[j], 0.2), border = F)
  lines(todseq, Delta[,j], lwd = 2, col = color[j]) 
}
dev.off()


