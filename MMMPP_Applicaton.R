# load("./data/bluewhalelist.RData")
# 
# data = datalist0[[24]]
# 
# head(data)
# 
# mllk = function(theta.star, divetime, N=2){
#   lambda = exp(theta.star[1:N]) # state specific rates
#   Q = diag(N) # generator matrix
#   Q[!Q] = exp(theta.star[N+1:(N*(N-1))])
#   diag(Q) = 0
#   diag(Q) = -rowSums(Q)
#   # delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
#   
#   l = 0
#   for(i in 1:length(divetime)){
#     Qube = LaMa::tpm_cont(Q-diag(lambda), divetime[[i]])
#     allprobs = matrix(lambda, nrow = length(divetime[[i]]+1), ncol = N, byrow = T)
#     allprobs[1,] = 1
#     l = l + LaMa::forward_g(delta, Qube, allprobs)
#   }
#   -l
# }
# 
# 
# 
# load("./data/minkewhales.RData")
# # divetime = dive2
# divetime = lapply(dive2, function(x) x/ 60)
# which(unlist(lapply(dive2, length))>=60)
# 
# lambda = 1 / c(quantile(divetime[[1]], 0.1), quantile(divetime[[1]], 0.9))
# qs = mean(lambda)/c(10,10)
# 
# theta.star = c(log(lambda), log(qs))
# 
# mod = nlm(mllk, theta.star, divetime = divetime, 
#           print.level = 2, iterlim = 1000, stepmax = 10, hessian = TRUE)
# lambda_hat = exp(mod$estimate[1:2])
# q_hat = exp(mod$estimate[3:4])
# 
# I = solve(mod$hessian)
# se = sqrt(diag(I))
# 
# round(1/lambda_hat, 1) 
# round(1/q_hat, 1)
# 
# 
# 
# 
# 
# 
# # Handball ----------------------------------------------------------------
# 
# # load data
# handball_data = read.csv("./data/handball_data.csv")
# 
# handball_data$time2 = as.POSIXct(handball_data$Time, format = "%M:%S")
# 
# library(tidyverse)
# handball_data = handball_data %>% fill(time2, .direction = "down")
# 
# # within each match, sort by time
# datalist = split(handball_data, handball_data$matchID)
# for(i in 1:length(unique(handball_data$matchID))){
#   datalist[[i]] = datalist[[i]] %>% arrange(time2)
# }
# 
# handball_data = bind_rows(datalist)
# 
# timediff = rep(NA, nrow(handball_data))
# for(i in 1:length(unique(handball_data$matchID))){
#   ind = which(handball_data$matchID == unique(handball_data$matchID)[i])
#   timediff[ind] = c(NA, diff(handball_data$time2[ind]))
# }
# 
# library(LaMa)
# trackInd = calc_trackInd(handball_data$matchID)
# 
# timediff[which(is.na(timediff))] = 0
# 
# mllk = function(theta.star, timediff, N=2, trackInd){
#   lambda = exp(theta.star[1:N]) # state specific rates
#   Q = diag(N) # generator matrix
#   Q[!Q] = exp(theta.star[N + 1:(N*(N-1))])
#   diag(Q) = 0
#   diag(Q) = -rowSums(Q)
#   delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
#   
#   Qube = LaMa::tpm_cont(Q-diag(lambda), timediff) 
#   # Qube = array(dim = c(N,N,length(timediff)))
#   # for(t in 1:length(timediff)){
#   #   Qube[,,t] = expm::expm(Q-diag(lambda)*timediff[t])
#   # }
#   
#   # first slice per game is identity, but will be ignored by forward_g()
#   allprobs = matrix(lambda, nrow = length(timediff), ncol = N, byrow = TRUE)
#   allprobs[1,] = 1
#   LaMa::forward_g(delta, Qube, allprobs, trackInd)
# 
# }
# 
# lambda = 1 / c(quantile(timediff, 0.2, na.rm = T), quantile(timediff, 0.5, na.rm = T), quantile(timediff, 0.9, na.rm = T))
# qs = rep(mean(lambda) / 10, 6)
# 
# theta.star = as.numeric(c(log(lambda), log(qs)))
# N = 3
# 
# mod = nlm(mllk, theta.star, timediff = timediff, N = 3, trackInd = trackInd,
#           print.level = 2, iterlim = 1000, stepmax = 5, hessian = TRUE)
# 
# 
# 
# 
# 
# 
# # Other whale data set ----------------------------------------------------
# 
# whale = read.csv("/Users/jan-ole/Downloads/beaked_whale_baird.csv")
# 
# dim(whale)
# head(whale)
# 
# depth = whale[,2]
# plot(depth[1:10000], type = "l", ylim = c(0,50))
# 
# surface = which(depth < 10)
# whale$surfacing = 0
# whale$surfacing[surface] = 1
# 
# 
# dives = rle(whale$surfacing)
# divelengths = dives$lengths[which(dives$values == 0)]
# surfacetimes = dives$lengths[which(dives$values == 1)]
# 
# divelengths
# 
# plot(depth[sum(divelengths[1:7]) + sum(surfacetimes[1:7]) - 50 + 1:500], type = "l", ylim = c(0,20))
# abline(h = 10)
# 
# 
# ind = which(divelengths > 5)
# data = data.frame(divetime = divelengths, surfacetime = surfacetimes)[ind,]
# data$divetime = c(NA, data$divetime[-nrow(data)])
# 
# mllk = function(theta.star, y, timediff, N){
#   lambda = exp(theta.star[1:N]) # state specific rates
#   rate = exp(theta.star[N+1:N])
# 
#   Q = diag(N) # generator matrix
#   Q[!Q] = exp(theta.star[2*N+1:(N*(N-1))])
#   diag(Q) = 0
#   diag(Q) = -rowSums(Q)
#   delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
#   Qube = LaMa::tpm_cont(Q-diag(lambda), timediff) # exp((Q-Lambda)*deltat)
#   
#   allprobs = matrix(1, length(y), N)
#   for(j in 1:N){
#     allprobs[,j] = dexp(y, rate = rate[j])
#   }
#   allprobs[-1,] = allprobs[-1,] * matrix(lambda, length(y)-1, N, byrow = T)
#   -LaMa::forward_g(delta, Qube, allprobs)
# }
# 
# y = data$surfacetime
# timediff = data$divetime[-1]
# 
# lambda = 1 / c(quantile(data$divetime[-1], 0.2), quantile(data$divetime[-1], 0.8))
# rate = 1 / c(quantile(data$surfacetime, 0.2), quantile(data$surfacetime, 0.8))
# qs = rep(mean(lambda)/10, 2)
# 
# theta.star = as.numeric(c(log(lambda), log(rate), log(qs)))
# 
# mod = nlm(mllk, theta.star, y = data$surfacetime, timediff = data$divetime[-1], N = 2,
#           print.level = 2, iterlim = 1000, stepmax = 10, hessian = TRUE)
# 
# sqrt(diag(solve(mod$hessian)))
# exp(mod$estimate)
# 



## loading Minkewhale data
mink = scan("/Users/jan-ole/Downloads/divepatterns.txt")

splitInd = which(mink == 1)
# split into list by splitInd
dive = split(mink, cumsum(mink == 1))

# filter for more than 60 observations
surftimes = dive[which(unlist(lapply(dive2, length))>=60)]

# 7 track have more than 60 observations
# we only consider whales 2, 3, 4 and 5 for the application, as the other whales
# lead to the special case of a homogeneous Poisson process -> absorbing state

surftimes = surftimes[2:5]


# little EDA
pdf("./figs/mmpp_divetimes.pdf", width = 7, height = 4)

par(mfrow = c(4,1), mar = c(4.1,2,0,2) + 0.1)
for(i in 1:4){
  plot(surftimes[[i]], rep(1, length(surftimes[[i]])), ylim = c(0,1.5), type = "h", 
       yaxt = "n", xaxt = "n", xlab = "Dive time in minutes", ylab = "", bty = "n", xlim = c(0, 3500))
  axis(1, at = seq(0, 3500, 500), labels = seq(0, 3500, 500))
}
# plot(surftimes[[4]], rep(1, length(surftimes[[4]])), ylim = c(0,1.5), type = "h",
#      yaxt = "n", xaxt = "n", xlab = "Dive time in minutes", ylab = "", bty = "n", xlim = c(0, 3500))
# axis(1, at = seq(0, 3500, 500), labels = seq(0, 3500, 500))
dev.off()


# writing the likelihood function using LaMA
mllk = function(theta.star, timediff, N=2){
  lambda = exp(theta.star[1:N]) # state specific rates
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
  
  Qube = LaMa::tpm_cont(Q-diag(lambda), timediff)
  allprobs = matrix(lambda, nrow = length(timediff+1), ncol = N, byrow = T)
  allprobs[1,] = 1
  
  -LaMa::forward_g(delta, Qube, allprobs)
}


## fitting the models
lambda = 1 / c(quantile(unlist(surftimes), 0.4), quantile(unlist(surftimes), 0.6))
qs = rep(mean(lambda)/5, 2)

# initial values
theta.star = log(c(lambda, qs))

mods = list()
for(i in 1:length(surftimes)){
  mods[[i]] = nlm(mllk, theta.star, timediff = diff(surftimes[[i]]), N = 2,
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

round(colMeans(1/Lambdas), 2)
round(colMeans(1/Qs), 2)

## nur 4 gute Wale nehmen für Application.
# Verweis auf JASA paper
