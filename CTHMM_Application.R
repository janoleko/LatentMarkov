
# Loading packages --------------------------------------------------------

# install.packages("msm")
library(msm)
# install.packages("LaMa")
library(LaMa)


# Data --------------------------------------------------------------------

data = fev # data from the msm package

## little EDA
nrow(fev) 
# 5896 observations
colnames(fev)
length(unique(fev$ptnum))
# 203 patients
datalist = split(fev, fev$ptnum)
sum(sapply(datalist, function(X) sum(X$fev == 999)))
# 96 patients died

plot(data$days, data$fev, xlab = "time", ylab = "fav", bty = "n", pch = 20)
hist(data$fev, prob = T, border = "white", main = "Histogram of fev", xlab = "fev", breaks = 100)

# 999 codes death: we will just use it as is, as it is so far from all other values


# Model with state transition structure -----------------------------------

# Loading the negative log-likelihood function ----------------------------

source("likelihood_functions.R")

# ### long version
# mllk_ct_slow = function(theta.star, X){
#   beta = matrix(theta.star[1:4], ncol = 2)
#   sigma = exp(theta.star[4+1:2])
#   # structured generator matrix
#   Q = matrix(0, 3, 3)
#   Q[1,2:3] = exp(theta.star[6+1:2])
#   Q[2,3] = exp(theta.star[9])
#   diag(Q) = -rowSums(Q)
#   delta = c(exp(theta.star[10:11]),0)
#   delta = delta/sum(delta)
#   # for the rest: loop over patients
#   l = 0 # initialize log likelihood
#   ptnums = unique(X$ptnum) # patient numbers
#   for(i in ptnums){
#     X_p = X[which(X$ptnum==i),]
#     n_p = nrow(X_p)
#     timediff = diff(X_p$days)
#     Qube = array(dim = c(3,3,n_p-1))
#     for(t in 1:(n_p-1)){
#       Qube[,,t] = expm::expm(Q*timediff[t])
#     }
#     allprobs = matrix(1, nrow = n_p, ncol = 3)
#     ind = which(!is.na(X_p$fev))
#     for(j in 1:2){
#       allprobs[ind,j] = dnorm(X_p$fev[ind], beta[j,1]+beta[j,2]*X_p$acute[ind], sigma[j])
#     }
#     allprobs[,3] = 0
#     allprobs[which(X_p$fev==999),] = c(0,0,1)
#     # forward algorithm to calculate the log-likelihood recursively
#     foo = delta%*%diag(allprobs[1,])
#     phi = foo/sum(foo)
#     l_p = log(sum(foo))
#     for(t in 2:n_p){
#       foo = phi%*%Qube[,,t-1]%*%diag(allprobs[t,])
#       phi = foo/sum(foo)
#       l_p = l_p + log(sum(foo))
#     }
#     l = l + l_p
#   }
#   return(-l)
# }
# 
# ### short version using the package LaMa
# mllk_ct_fast = function(theta.star, X){
#   beta = matrix(theta.star[1:4], ncol = 2)
#   sigma = exp(theta.star[4+1:2])
#   # structured generator matrix
#   Q = matrix(0, 3, 3)
#   Q[1,2:3] = exp(theta.star[6+1:2])
#   Q[2,3] = exp(theta.star[9])
#   diag(Q) = -rowSums(Q)
#   delta = c(exp(theta.star[10:11]),1)
#   delta = delta/sum(delta)
#   # for the rest: loop over patients
#   l = 0 # initialize log likelihood
#   ptnums = unique(X$ptnum) # patient numbers
#   for(i in ptnums){
#     X_p = X[which(X$ptnum==i),]
#     n_p = nrow(X_p)
#     timediff = diff(X_p$days)
#     Qube = LaMa::tpm_cont(Q, timediff) # exp(Q*dt)
#     allprobs = matrix(1, nrow = n_p, ncol = 3)
#     ind = which(!is.na(X_p$fev))
#     for(j in 1:2){
#       allprobs[ind,j] = dnorm(X_p$fev[ind], beta[j,1]+beta[j,2]*X_p$acute[ind], sigma[j])
#     }
#     allprobs[,3] = 0
#     allprobs[which(X_p$fev==999),] = c(0,0,1)
#     # forward algorithm to calculate the log-likelihood recursively
#     l = l + LaMa::forward_g(delta, Qube, allprobs)
#   }
#   return(-l)
# }


# Model fitting -----------------------------------------------------------

# initial parameters for numerical optimisation
beta0 = c(100, 50) # intercept
sd0 = c(15,15)
qs0 = rep(1/100, 3)

# working scale initial parameter vector
theta.star0 = c(beta0, rep(0,2), log(sd0), log(qs0), rep(0,2))

# fitting the model by numerically optimising the log-likelihood with nlm
mod_fev5 = nlm(mllk_ct_fast, theta.star0, X = data,
               iterlim = 1000, print.level = 2, stepmax = 10, hessian = TRUE)
I_fev5 = solve(mod_fev5$hessian)

theta.star = mod_fev5$estimate
beta = matrix(theta.star[1:4], ncol = 2)

beta1CI = cbind(beta[,2] - 1.96 * sqrt(diag(I_fev5)[3:4]),
                beta[,2] + 1.96 * sqrt(diag(I_fev5)[3:4]))

round(beta1CI, 3)

colnames(beta) = c("intercept", "accute inf.")
rownames(beta) = c("healthy state", "sick state")
sqrt(diag(I_fev5))[1:4]
sigma = exp(theta.star[4+1:2])
# structured generator matrix
Q = matrix(0, 3, 3)
Q[1,2:3] = exp(theta.star[6+1:2])
Q[2,3] = exp(theta.star[9])
diag(Q) = -rowSums(Q)
rownames(Q) = colnames(Q) = c("healthy", "sick", "dead")
# average number of years it takes to get to next diseasse stage/ death:
round((1/Q)/365,2)

# initial distribution (not too interesting here)
delta = c(exp(theta.star[10:11]),0)
delta = delta/sum(delta)

# decoding states
ptnums = unique(data$ptnum) # patient numbers
states = rep(NA, nrow(data))
for(i in ptnums){
  X_p = data[which(data$ptnum==i),]
  n_p = nrow(X_p)
  timediff = diff(X_p$days)
  Qube = LaMa::tpm_cont(Q, timediff) # exp(Q*dt)
  allprobs = matrix(1, nrow = n_p, ncol = 3)
  ind = which(!is.na(X_p$fev))
  for(j in 1:2){
    allprobs[ind,j] = dnorm(X_p$fev[ind], beta[j,1]+beta[j,2]*X_p$acute[ind], sigma[j])
  }
  allprobs[,3] = 0
  allprobs[which(X_p$fev==999),] = c(0,0,1)
  # forward algorithm to calculate the log-likelihood recursively
  states[which(data$ptnum==i)] = viterbi_g(delta, Qube, allprobs)
}

# empirical state frequencies for accute = 0 and accute = 1
(delta0 = c(sum(states[which(data$acute==0)]==1), 
           sum(states[which(data$acute==0)]==2))/nrow(data[which(data$acute==0),]))
(delta1 = c(sum(states[which(data$acute==1)]==1), 
            sum(states[which(data$acute==1)]==2))/nrow(data[which(data$acute==1),]))
# probabilty of being in the diseased state is much higher after infection



# Figure as shown in the paper --------------------------------------------

color = c("orange", "deepskyblue")
pdf("./figs/cthmm_marginal.pdf", width = 8, height = 4)
# marginal distribution without and with accute infection in the last 14 days
par(mfrow = c(1,2))
hist(data$fev[which(data$acute==0)], prob = T, border = "white", ylab = "density",
     main = "accute = 0", xlab = "fev", breaks = 150, ylim = c(0,0.02), xlim = c(0,150))
for(j in 1:2) curve(delta0[j]*dnorm(x, beta[j,1], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta0[1]*dnorm(x, beta[1,1], sigma[1])+delta0[2]*dnorm(x, beta[2,1], sigma[2]),
      add = T, lwd = 2, lty = 2)

hist(data$fev[which(data$acute==1)], prob = T, border = "white",  ylab = "density",
     main = "acute = 1", xlab = "fev", breaks = 50, ylim = c(0,0.02), xlim = c(0,150))
for(j in 1:2) curve(delta1[j]*dnorm(x, beta[j,1]+beta[j,2], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta1[1]*dnorm(x, beta[1,1]+beta[1,2], sigma[1])+delta1[2]*dnorm(x, beta[2,1]+beta[2,2], sigma[2]),
      add = T, lwd = 2, lty = 2)
legend("topright", lwd = 2, col = color, legend = paste("state", 1:2), bty = "n") 
dev.off()


