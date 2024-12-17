# This case study illustrates continuous-time SSMs by analyzing seven-metre throws in Handball to find evidence for a hot hand effect.

# Loading packages --------------------------------------------------------

# install.packages("tidyverse")
library(tidyverse)
# install.packages("LaMa")
library(LaMa)


# Loading data ------------------------------------------------------------

# read in the data
data = read.csv("./data/handball_data_7m.csv")

# some data cleaning
n = nrow(data)
data = data[n:1, ] # reverse the order of the data

unique_players = unique(data$player)
data$player = as.factor(data$player)
data = split(data, data$player)
data = bind_rows(data[which(lapply(data, nrow) > 10)]) # only consider players with more than 10 throws


# Summary statistics ------------------------------------------------------

# number of matches
length(unique(data$matchID))

# check how many throws in a game
view = data %>% group_by(player, matchID) %>% summarise(n = n())
# View(view)
max(view$n) # maximum of 9 throws by one player in one game
summary(view$n) # min of 1 throw, median 2 and average 2.4 throws per game


# Preparing data for model fitting ----------------------------------------

data$uID = paste0(data$matchID, "_", data$player) # creating a unique ID for each player and game

# We have a large amount of short tracks - one track is a series of throws by one player in one game
# we have to create the indices at which tracks starts for LaMa to sum the individual log likelihoods later
trackInds = list()
data$time2 = as.POSIXct(data$time2)
data$player = as.factor(as.character(data$player))
players = unique(data$player)
data = split(data, data$player)
for(i in 1:length(players)){
  if(players[i] == "Uwe Gensheimer"){
    data[[i]][1:2, ] = data[[i]][2:1, ] # throws by Uwe Gensheimer were in wrong order
  }
  data[[i]]$timediff = c(NA, diff(data[[i]]$time2))
  trackInds[[i]] = calc_trackInd(data[[i]]$uID)
  data[[i]]$timediff[trackInds[[i]]] = NA
}
data = bind_rows(data)

data$timediff = data$timediff / 60 # timediff in minutes not seconds

# View(data)

# filter out games with only one throw (as we need a minumu of 2 for a "track" -> at least one transition)
uIDs = unique(data$uID)
data$uID = as.factor(data$uID)
data = split(data, data$uID)
data = bind_rows(data[which(lapply(data, nrow) > 2)])

# finally filter out tracks for which time diffs are missing
data = split(data, data$uID)
data = bind_rows(data[which(unlist(lapply(data, function(X) sum(is.na(X[, 12])) <= 1)))])

nrow(data) # 1131 total throws
length(unique(data$player)) # by 46 players
length(unique(data$matchID)) # in 205 matches


# Loading the negative log-likelihood function ----------------------------

source("likelihood_functions.R")


# To avoid redundant calculations, we calculate the unique time differences and match them to the corresponding rows
deltat = sort(unique(data$timediff))
# we do this, as there are less unique time differences than overall transitions in the data set, this way saving computation time.

# creating a column in the data set to match the unique time differences to the corresponding row
data$match2Array = NA
for(i in 1:nrow(data)){
  if(!is.na(data$timediff[i])){
    data$match2Array[i] = which(deltat == data$timediff[i])
  }
}

# this vector defines where each single track starts and is necessary for LaMa when summing the individual log-likelihoods
trackInds = calc_trackInd(as.character(data$uID))


# Fitting the model -------------------------------------------------------

# initial parameter values for numerical optimization
theta0 = 0.2
sigma0 = 0.25
beta0 = qlogis(0.8)

theta.star = c(log(theta0), log(sigma0), plogis(beta0))

modOU = nlm(mllk_ctSSM_fast, theta.star, X = data, deltat = deltat, bm = 3.5, m = 250,
            print.level = 2, iterlim = 1000, hessian = TRUE)

# this is still rather slow, as calculating the transition matrices takes the most time here, 
# therefore running the forward algorithm in C++ does not help much

(theta = exp(modOU$estimate[1])) # mean reversion strength of the OU process
(sigma = exp(modOU$estimate[2])) # variability of the OU process
(beta = modOU$estimate[3]) # Intercept of the state-dependent process on logit-scale
plogis(beta) # scoring probability when state process in equilibrium

# stationary distribution
qnorm(0.01, 0, sqrt(sigma^2 / (2 * theta))) # initial distribution = stationary distribution of OU process

# standard errors
I = solve(modOU$hessian)
se = sqrt(diag(I))

(thetaCI = exp(modOU$estimate[1] + 1.96 * se[1]*c(-1,1)))
(sigmaCI = exp(modOU$estimate[2] + 1.96 * se[2]*c(-1,1)))
(betaCI = modOU$estimate[3] + 1.96 * se[3]*c(-1,1))
(probCI = plogis(modOU$estimate[3] + 1.96 * se[3]*c(-1,1)))

# estimated parameters and confidence intervals
round(theta, 3)
round(thetaCI, 3)

round(sigma, 3)
round(sigmaCI, 3)

round(beta, 3)
round(betaCI, 3)

round(plogis(beta), 4)


# ACF deduced from OU process
exp(-theta * 5) # 92% autocorrelation after one minute
exp(-theta * 60) # 35% autocorrelation after 60 minutes (entire game)


# sampling from fitted OU process

dt = 1/60 # one second intervals
n = 60 / dt

color = c("orange", "deepskyblue", "seagreen2", "plum", "navy", "firebrick3")

pdf("./figs/simOU_process.pdf", width = 8, height = 5)

set.seed(13)
# one realization
plot(NA, xlim = c(0, 60), ylim = c(-2.5, 2.5), xlab = "time (minutes)", ylab = "hotness", bty = "n")
x = numeric(n)
time = 1:n * dt
x[1] = rnorm(1, 0, sqrt(sigma^2 / (2 * theta)))
for(t in 2:n){
  x[t] = rnorm(1, mean = exp(-theta * dt) * x[t-1], 
               sd = sqrt((1 - exp(-2 * theta * dt)) * sigma^2 / (2 * theta)))
}
lines(time, x, type = "l", col = "black", lwd = 1)

# adding a few realizations
k = 5
for(i in 1:k){
  x = numeric(n)
  time = 1:n *dt
  x[1] = rnorm(1, 0, sqrt(sigma^2 / (2 * theta)))
  for(t in 2:n){
    x[t] = rnorm(1, mean = exp(-theta * dt) * x[t-1], 
                 sd = sqrt((1 - exp(-2 * theta * dt)) * sigma^2 / (2 * theta)))
  }
  lines(time, x, type = "l", col = color[i], lwd = 1)
}

dev.off()

