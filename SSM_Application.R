
# Downloading data --------------------------------------------------------

library(fHMM) # makes downloading financial data really easy

# data = download_data("^GSPC")
data = download_data("BTC-USD")
data$return = c(NA, diff(log(data$Close)))
data$Date[1]; data$Date[nrow(data)]


## EDA
plot(as.POSIXct(data$Date), data$Close, type = "l", 
     bty = "n", xlab = "date", ylab ="close")

plot(as.POSIXct(data$Date), data$return, type = "l", 
     bty = "n", xlab = "date", ylab = "returns")

hist(data$return, xlim = c(-0.1,0.1), prob = TRUE, border = "white", 
     breaks = 100, main = "Histogram of returns", xlab = "returns")
# heavier tails than a normal distribution

# Likelihood functions ----------------------------------------------------

### long version
mllk_ssm_long = function(theta.star, y, bm, m){
  # transforming parameters to natural
  phi = plogis(theta.star[1])
  sigma = exp(theta.star[2])
  beta = exp(theta.star[3])
  # defining intervals and gridpoints for numerical integration
  b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # approximating tpm resulting from midpoint quadrature
  Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
  Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
  delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
  # approximating state-dependent density based on midpoints
  allprobs = matrix(1, length(y), m)
  ind = which(!is.na(y))
  allprobs[ind] = t(sapply(y[ind], dnorm, mean = 0, sd = beta * exp(bstar/2)))
  # forward algorithm to calculate the approximate log-likelihood recursively
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo / sum(foo)
  for(t in 2:nrow(X)){
    foo = phi%*%Gamma%*%diag(allprobs[t,])
    l = l + log(sum(foo))
    phi = foo / sum(foo)
  }
  return(-l)
}

### short version using the package LaMa
library(LaMa)
# with the building-block functions from the R package LaMa, we can write the 
# negative log-likelihood function much shorter and make it fuch faster (by a factor of 10-20).
mllk_ssm_short = function(theta.star, y, bm, m){
  # transforming parameters to natural
  phi = plogis(theta.star[1])
  sigma = exp(theta.star[2])
  beta = exp(theta.star[3])
  mu = theta.star[4]
  # defining intervals and gridpoints for numerical integration
  b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # approximating tpm resulting from midpoint quadrature
  Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
  Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
  delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
  # approximating state-dependent density based on midpoints
  allprobs = matrix(1, length(y), m)
  ind = which(!is.na(y))
  allprobs[ind,] = t(sapply(y[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))
  # forward algorithm to calculate the approximate log-likelihood recursively
  -forward(delta, Gamma, allprobs)
}


# Model fitting -----------------------------------------------------------

# initial parameters for numerical optimisation
hist(data$return, xlim = c(-0.1,0.1), prob = TRUE, border = "white", breaks = 100)
curve(dnorm(x, 0, sd(data$return[-1])), add = TRUE)
curve(dnorm(x, 0, 0.005), add = TRUE)
curve(dnorm(x, 0, 0.05), add = TRUE)

phi0 = 0.9
sigma0 = 0.65
beta0 = 0.02
mu0 = 0.001
qnorm(0.025, 0, sigma0/sqrt(1-phi0^2)) # 0.025 quantiale of stationary distribution

# working scale initial parameter vector
theta.star0 = c(qlogis(phi0), log(sigma0), log(beta0), mu0)

mod_BC1 = nlm(mllk_ssm_short, theta.star0, y = data$return, bm = 3.5, m = 200,
                    print.level = 2, hessian = TRUE, stepmax = 10)

# obtaining the estimated parameters
(phi = plogis(mod_BC1$estimate[1]))
(sigma = exp(mod_BC1$estimate[2]))
(beta = exp(mod_BC1$estimate[3]))
(mu = mod_BC1$estimate[4]) # expected return slightly positive
qnorm(0.025, 0, sigma/sqrt(1-phi^2))

I = solve(mod_BC1$hessian)
# standard deviations with delta method
sd_phi = sqrt(I[1,1]*(
  exp(mod_BC1$estimate[1])/
    (1+exp(mod_BC1$estimate[1]))^2
  )^2)
sd_sigma = sqrt(I[2,2]*exp(mod_BC1$estimate[2])^2)
sd_beta = sqrt(I[3,3]*exp(mod_BC1$estimate[3])^2)
sd_mu = I[4,4]

round(sd_phi, 3)
round(sd_sigma, 3)
round(sd_beta, 3)
round(sd_mu, 4)


# Results -----------------------------------------------------------------

## state decoding
bm = 3.5; m = 200
# defining intervals and gridpoints for numerical integration
b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
h = b[2]-b[1] # interval width
bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
# approximating tpm resulting from midpoint quadrature
Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
# approximating state-dependent density based on midpoints
allprobs = matrix(1, length(data$return), m)
ind = which(!is.na(data$return))
allprobs[ind,] = t(sapply(data$return[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))

mod_BC1$rawstates = viterbi(delta, Gamma, allprobs)
mod_BC1$states = bstar[mod_BC1$rawstates]

pdf("./figs/ssm_decoded.pdf", width = 8, height = 4.5)
par(mfrow = c(1,1), mar = c(5,4,4,4)+0.1)
plot(as.POSIXct(data$Date), data$return, yaxt = "n", 
     type = "l", bty = "n", ylab = "", ylim = c(-0.5, 0.2), xlab = "date")
lines(as.POSIXct(data$Date), 2*beta*exp(mod_BC1$states/2)-0.45, type = "l", col = "orange")
axis(2, at = seq(-0.2, 0.2, by = 0.1), labels = seq(-0.2, 0.2, by = 0.1))
axis(4, at = seq(-0.5, -0.1, by = 0.1), labels = 0.5*(seq(-0.5, -0.1, by = 0.1)+0.45))
mtext("volatility", side=4, line=3, at = -0.3)
mtext("return", side=2, line=3, at = 0)
dev.off()

# Calculating VaR ---------------------------------------------------------

library(expm)

# calculating the scaled forward variables
# first unscaled, but on log scale
lalpha = LaMa:::logalpha(delta, array(Gamma, dim = c(m,m,nrow(allprobs)-1)), allprobs)

# function to evaluate the density of the forecast distribution

dforecast = function(x, lalpha, t, Gamma, step = 1, mu, beta, bm=3.5){
  # calculating the scaled forward variable at time point of interest
  lalphaT = lalpha[t,]
  ma = max(lalphaT)
  phiT = exp(lalphaT-ma)/sum(exp(lalphaT-ma))
  # details for numerical integration
  m = dim(Gamma)[1]
  b = seq(-bm, bm, length = m+1)
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # calculating state-dependent probs for all values in x
  allprobs = t(sapply(x, dnorm, mean = mu, sd = beta * exp(bstar/2)))
  out = rep(NA, length(x))
  # forecast distribution
  for(i in 1:length(x)){
    out[i] = sum(phiT%*%(Gamma%^%step)%*%diag(allprobs[i,]))
  }
  out
}

xseq = seq(-0.25, 0.25, length=400)
h2 = 0.5/400

for(t in (nrow(data)-20) : (nrow(data)-1)){
  pred = dforecast(xseq, lalpha, t=t, Gamma, mu=mu, beta=beta, bm=bm)
  pred = pred/sum(pred) # we have to normalize
  plot(xseq, pred/h, type = "l", main = t, bty = "n", lwd = 2)
  # deviding by h to plot the density
  abline(v=data$return[t+1])
  # plotting the observed data point
  Sys.sleep(0.4)
}

## actually calculating VaR based on forecast distribution
alpha = 0.025 # confidence level
pred = dforecast(xseq, lalpha, t=nrow(data), Gamma, mu=mu, beta=beta, bm=bm)
pred = pred / sum(pred)
smaller = which(cumsum(pred)<alpha) # find the values that have smaller cum prob
VaR = -xseq[smaller[length(smaller)]]

library(scales)

pdf("./figs/VaR.pdf", width = 7, height = 4)
dens = pred/h2
ind = 1:which.min((xseq+VaR)^2)
plot(NA, main = paste("Value at Risk =",round(VaR,3)), xlim = c(-0.2,0.2), ylim = c(0, 15),
     xlab = "return", ylab = "density", bty = "n")
polygon(c(xseq[ind], rev(xseq[ind])), c(rep(0, length(ind)), rev(dens[ind])), 
        col = alpha("orange",0.4), border = F)
lines(xseq, dens, lwd = 2)
abline(v = -VaR, col = "orange", lwd = 2)
dev.off()
