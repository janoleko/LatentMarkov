
# Loading packages --------------------------------------------------------

# install.packages("fHMM")
library(fHMM) # makes downloading financial data really easy
# install.packages("LaMa")
library(LaMa)
# install.packages("expm")
library(expm) # to compute matrix powers


# Downloading data --------------------------------------------------------

data = download_data("BTC-USD", to = "2024-05-31")
data$return = c(NA, diff(log(data$Close)))
data$Date[1]; data$Date[nrow(data)]


# Some EDA ----------------------------------------------------------------

plot(as.POSIXct(data$Date), data$Close, type = "l", 
     bty = "n", xlab = "date", ylab ="close")

plot(as.POSIXct(data$Date), data$return, type = "l", 
     bty = "n", xlab = "date", ylab = "returns")

hist(data$return, xlim = c(-0.1,0.1), prob = TRUE, border = "white", 
     breaks = 100, main = "Histogram of returns", xlab = "returns")
# heavier tails than a normal distribution -> stochastic volatility


# Loading the negative log-likelihood function ----------------------------

source("likelihood_functions.R")

### long version in base R
# mllk_ssm_slow = function(theta.star, y, bm, m){
#   # transforming parameters to natural
#   phi = plogis(theta.star[1])
#   sigma = exp(theta.star[2])
#   beta = exp(theta.star[3])
#   # defining intervals and gridpoints for numerical integration
#   b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
#   h = b[2]-b[1] # interval width
#   bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
#   # approximating tpm resulting from midpoint quadrature
#   Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
#   delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
#   # approximating state-dependent density based on midpoints
#   allprobs = matrix(1, length(y), m)
#   ind = which(!is.na(y))
#   allprobs[ind] = t(sapply(y[ind], dnorm, mean = 0, sd = beta * exp(bstar/2)))
#   # forward algorithm to calculate the approximate log-likelihood recursively
#   foo = delta%*%diag(allprobs[1,])
#   l = log(sum(foo))
#   phi = foo / sum(foo)
#   for(t in 2:nrow(X)){
#     foo = phi%*%Gamma%*%diag(allprobs[t,])
#     l = l + log(sum(foo))
#     phi = foo / sum(foo)
#   }
#   return(-l)
# }

### short version using the package LaMa
# mllk_ssm_fast = function(theta.star, y, bm, m){
#   # transforming parameters to natural
#   phi = plogis(theta.star[1])
#   sigma = exp(theta.star[2])
#   beta = exp(theta.star[3])
#   mu = theta.star[4]
#   # defining intervals and gridpoints for numerical integration
#   b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
#   h = b[2]-b[1] # interval width
#   bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
#   # approximating tpm resulting from midpoint quadrature
#   Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
#   # Gamma = Gamma / rowSums(Gamma)
#   delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
#   # approximating state-dependent density based on midpoints
#   allprobs = matrix(1, length(y), m)
#   ind = which(!is.na(y))
#   allprobs[ind,] = t(sapply(y[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))
#   # forward algorithm to calculate the approximate log-likelihood recursively
#   -forward(delta, Gamma, allprobs)
# }


# Model fitting -----------------------------------------------------------

# eyeballing initial parameters for numerical optimisation
par(mfrow = c(1,1))
hist(data$return, xlim = c(-0.1,0.1), prob = TRUE, border = "white", breaks = 100)
curve(dnorm(x, 0, sd(data$return[-1])), add = TRUE)
curve(dnorm(x, 0, 0.005), add = TRUE)
curve(dnorm(x, 0, 0.05), add = TRUE)

phi0 = 0.9
sigma0 = 0.65
beta0 = 0.02
mu0 = 0.001
qnorm(0.025, 0, sigma0/sqrt(1-phi0^2)) # 0.025 quantiale of stationary distribution -> to find a good essential range

# working scale initial parameter vector
theta.star0 = c(qlogis(phi0), log(sigma0), log(beta0), mu0)

mod_BC1 = nlm(mllk_ssm_fast, theta.star0, y = data$return, bm = 3.5, m = 200,
                    print.level = 2, hessian = TRUE, stepmax = 10)

# obtaining the estimated parameters
(phi = plogis(mod_BC1$estimate[1]))
(sigma = exp(mod_BC1$estimate[2]))
(beta = exp(mod_BC1$estimate[3]))
(mu = mod_BC1$estimate[4]) # expected return slightly positive
qnorm(0.025, 0, sigma/sqrt(1-phi^2))

I = solve(mod_BC1$hessian)
sd = sqrt(diag(I))

# confidence intervals

phiCI = plogis(mod_BC1$estimate[1] + 1.96 * c(-1,1) * sd[1])
sigmaCI = exp(mod_BC1$estimate[2] + 1.96 * c(-1,1) * sd[2])
betaCI = exp(mod_BC1$estimate[3] + 1.96 * c(-1,1) * sd[3])
muCI = mod_BC1$estimate[4] + 1.96 * c(-1,1) * sd[4]

# estimated parameters and confidence intervals

round(beta, 3)
round(betaCI, 3)

round(phi, 3)
round(phiCI, 3)

round(sigma, 3)
round(sigmaCI, 3)

round(mu, 4)
round(muCI, 4)


# Visualizing the results -------------------------------------------------

## state decoding
bm = 3.5; m = 200
# defining intervals and gridpoints for numerical integration
b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
h = b[2]-b[1] # interval width
bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
# approximating tpm resulting from midpoint quadrature
Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
# approximating state-dependent density based on midpoints
allprobs = matrix(1, length(data$return), m)
ind = which(!is.na(data$return))
allprobs[ind,] = t(sapply(data$return[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))

mod_BC1$rawstates = LaMa::viterbi(delta, Gamma, allprobs)
mod_BC1$states = bstar[mod_BC1$rawstates]


alpha = (1-(pnorm(1)-pnorm(-1)))/2 # one standard deviation
stateprobs = stateprobs(delta, Gamma, allprobs)
cumstateprobs = t(apply(stateprobs, 1, cumsum))
lower = apply(cumstateprobs, 1, function(x){
  ind = which(x <= alpha)
  ind[length(ind)]
})
med = apply(cumstateprobs, 1, function(x) ind = which(x >= 0.5)[1] )
Mean = rowSums(stateprobs * matrix(bstar, nrow = nrow(stateprobs), ncol = m, byrow = TRUE))
upper = apply(cumstateprobs, 1, function(x) ind = which(x >= 1-alpha)[1] )
quant = cbind(bstar[lower], bstar[upper])

pdf("./figs/ssm_decoded.pdf", width = 6.5, height = 5)
par(mfrow = c(1,1), mar = c(4.5,4,2,4)+0.1)
date = as.POSIXct(data$Date)
ind = (length(date)-1000):length(date)
plot(date[ind], data$return[ind], yaxt = "n", 
     type = "l", bty = "n", ylab = "", ylim = c(-0.4, 0.2), xlab = "date")

lines(date[ind], (2*beta*exp(Mean/2)-0.4)[ind], type = "l", lwd = 1, col = "orange") # mean
Interval = 2*beta*exp(quant/2)-0.4
polygon(c(date[ind], rev(date[ind])), c(Interval[ind,1], rev(Interval[ind,2])), 
        col = scales::alpha("orange", 0.5), border = F)
axis(2, at = seq(-0.2, 0.2, by = 0.1), labels = seq(-0.2, 0.2, by = 0.1))
axis(4, at = seq(-0.4, -0.1, by = 0.1), labels = 0.5*(seq(-0.4, -0.1, by = 0.1)+0.4))
mtext("volatility", side=4, line=3, at = -0.3)
mtext("return", side=2, line=3, at = 0)
dev.off()


# Calculating VaR for model validation ------------------------------------

# defining training and test data
data$date = as.POSIXct(data$Date)
trainInd = which(data$date < "2021-01-01")
testInd = which(data$date >= "2021-01-01")

mod_BC2 = nlm(mllk_ssm_fast, theta.star0, y = data$return[trainInd], bm = 3.5, m = 200,
              print.level = 2, hessian = TRUE, stepmax = 10)

# obtaining the estimated parameters
phi = plogis(mod_BC2$estimate[1])
sigma = exp(mod_BC2$estimate[2])
beta = exp(mod_BC2$estimate[3])
mu = mod_BC2$estimate[4]

# calculating allprobs and Gamma for estimated parameters
bm = 3.5; m = 200
# defining intervals and gridpoints for numerical integration
b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
h = b[2]-b[1] # interval width
bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
# approximating tpm resulting from midpoint quadrature
Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
# approximating state-dependent density based on midpoints
allprobs = matrix(1, length(data$return), m)
ind = which(!is.na(data$return))
allprobs[ind,] = t(sapply(data$return[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))

# calculating the scaled forward variables
# first unscaled, but on log scale
lalpha = LaMa:::logalpha_cpp(allprobs, delta, array(Gamma, dim = c(m,m,nrow(allprobs)-1)))

# defining a function that evaluates the density of the forecast distribution
dforecast = function(x, delta, Gamma, lalpha, day, step = 1, mu, beta, bm){
  # calculating the scaled forward variable at time point of interest
  lalphaT = lalpha[day,]
  ma = max(lalphaT)
  phiT = exp(lalphaT-ma)/sum(exp(lalphaT-ma))
  # details for numerical integration
  m = dim(Gamma)[1]
  b = seq(-bm, bm, length = m+1)
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # calculating state-dependent probs for all values in x
  allprobs = t(sapply(x, dnorm, mean = mu, sd = beta * exp(bstar/2)))
  
  rowSums(matrix(phiT%*%(Gamma%^%step), nrow = nrow(allprobs), ncol = m, byrow = TRUE) * allprobs)
}

xseq = seq(-0.5, 0.5, length = 500) # range of the returns
mVaR = rep(NA, length(testInd)) # minus value at risk vector for each day in the test data

# defining the confidence level
alpha = 0.01

# calculating forecast distribution for each test day and comparing to actual observed return
for(i in 1:length(testInd)){
  pred = dforecast(xseq, delta, Gamma, lalpha, day = testInd[i]-1, mu=mu, beta=beta, bm=bm)
  pred = pred / sum(pred)
  smaller = which(cumsum(pred) < alpha) # find the values that have smaller cum prob
  mVaR[i] = xseq[max(smaller)]
}

round(sum(data$return[testInd] < mVaR) / length(testInd), 4)
