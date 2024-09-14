rm(list = ls()) # cleans environment

time_steps = 100 # 100 time steps 

f = 10 # number of fishes 

coordinates = 2*f 

mu = 50

alpha = 0.1

sigma = 1

fish_data <- read.csv("fishdata.csv")

head(fish_data)

fish0 = fish_data[fish_data$who == 0, ]

xcors = fish0$xcor # x coordinates of fish 0
ycors = fish0$ycor # y coordinates of fish 0

head(fish0)

log_likelihood_func_one_fish <- function(xcors, ycors, alpha, sigma, mu, f) {
  
  likelihood = (1/2)^(2*f)
  log_likelihood = log(likelihood) # first coordinate of fish 0 Bernoulli
  
  for (t in 2:time_steps) {
    xlikelihood = dnorm(xcors[t], xcors[t-1]+alpha*(mu-xcors[t-1]), sigma)
    ylikelihood = dnorm(ycors[t], ycors[t-1]+alpha*(mu-ycors[t-1]), sigma)
    log_likelihood = log_likelihood + log(xlikelihood) + log(ylikelihood)
    likelihood = likelihood * xlikelihood * ylikelihood
  }
  
  return (c(likelihood, log_likelihood))
  
}

log_likelihood_func <- function(fish_data, alpha, sigma, mu, fishes) {
  
  likelihood = (1/2)^(2*fishes) # initial x and y coordinate of each fish is a bernoulli(1/2) , there are 10 fishes so 1/2^20
  log_likelihood = log(likelihood) 
  
  for (f in 0:fishes) { # for each fish
    
    fish = fish_data[fish_data$who == f, ] # get only their data
    
    xcors = fish$xcor # get their x coordinates for each time step
    ycors = fish$ycor # get their y coordinates for each time step
    
    for (t in 2:time_steps) { # update likelihood for each time step for the given fish
      xlikelihood = dnorm(xcors[t], xcors[t-1]+alpha*(mu-xcors[t-1]), sigma)
      ylikelihood = dnorm(ycors[t], ycors[t-1]+alpha*(mu-ycors[t-1]), sigma)
      log_likelihood = log_likelihood + log(xlikelihood) + log(ylikelihood)
      likelihood = likelihood * xlikelihood * ylikelihood
    }
    
  }
  return (c(likelihood, log_likelihood))
  
}

log_likelihood_func(fish_data, alpha, sigma, 50, 10)

# alphas = seq(0, 0.5, by = 0.001)
# sigmas = seq(0, 2, by = 0.001)
# likelihoods = integer(length(alphas)*length(sigmas))
# log_likelihoods = integer(length(alphas)*length(sigmas))
# count = 1

# for (alpha in alphas) {
  # for (sigma in sigmas) {
    # a = log_likelihood_func(xcors, ycors, alpha, sigma, 50)
    # likelihoods[count] = a[1]
    # log_likelihoods[count] = a[2]
    # count = count + 1
  # }
# }

# print(likelihoods)
#print(log_likelihoods)

