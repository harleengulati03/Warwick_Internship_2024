rm(list = ls()) # cleans environment

fish_data <- read.csv("fishdata.csv")

time_steps = 100 # 100 time steps 

mu = 50

# alpha = 0.250

# sigma = 1.9

fishes = 1

log_likelihood_func = function(fish_data, time_steps, mu, alpha, sigma, fishes) {
  
  likelihood = (1/2)^(2*fishes) # initial x and y coordinate of each fish is a bernoulli(1/2) , there are 10 fishes so 1/2^20
  log_likelihood = log(likelihood) 
  
  fish = fish_data[fish_data$who == 0, ] # get only their data
  
  
  xcors = fish$xcor # get their x coordinates for each time step
  ycors = fish$ycor # get their y coordinates for each time step
  
  for (t in 2:time_steps) { # update likelihood for each time step for the given fish
    xlikelihood = dnorm(xcors[t], xcors[t-1]+alpha*(mu-xcors[t-1]), sigma)
    ylikelihood = dnorm(ycors[t], ycors[t-1]+alpha*(mu-ycors[t-1]), sigma)
    log_likelihood = log_likelihood + log(xlikelihood) + log(ylikelihood)
    likelihood = likelihood * xlikelihood * ylikelihood
  }
  
  return (log_likelihood)
}

# likelihood = (1/2)^(2*fishes) # initial x and y coordinate of each fish is a bernoulli(1/2) , there are 10 fishes so 1/2^20
# log_likelihood = log(likelihood) 
# 
# fish = fish_data[fish_data$who == 0, ] # get only their data
#   
# xcors = fish$xcor # get their x coordinates for each time step
# ycors = fish$ycor # get their y coordinates for each time step
# 
# for (t in 2:time_steps) { # update likelihood for each time step for the given fish
#   xlikelihood = dnorm(xcors[t], xcors[t-1]+alpha*(mu-xcors[t-1]), sigma)
#   ylikelihood = dnorm(ycors[t], ycors[t-1]+alpha*(mu-ycors[t-1]), sigma)
#   log_likelihood = log_likelihood + log(xlikelihood) + log(ylikelihood)
#   likelihood = likelihood * xlikelihood * ylikelihood
# }

alphas = seq(0, 0.5, length.out = 100)

sigmas = seq(0, 2, length.out = 100)

log_likelihoods_store = integer(length(alphas) * length(sigmas))
count = 0

for (alpha in alphas) {
  for (sigma in sigmas) {
    a = log_likelihood_func(fish_data, time_steps, mu, alpha, sigma, fishes)
    if (a > 0) {
      log_likelihoods_store[count] = 1
      count = count + 1
      }
  }
}

print("DONE")

print(max(log_likelihoods_store))


