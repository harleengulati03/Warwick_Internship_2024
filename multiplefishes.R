rm(list = ls()) # cleans environment

log_likelihood_func <- function(fish_data, time_steps, mu, alpha, sigma, fishes) {
  likelihood = (1/2)^(2*fishes) # initial x and y coordinate of each fish is a bernoulli(1/2) , there are 10 fishes so 1/2^20
  log_likelihood = log(likelihood)
  
  for (f in 0:fishes) { # for each fish
    
    fish = fish_data[fish_data$who == f, ] # get only their data over 100 time steps
    
    xcors = fish$xcor # get their x coordinates for each time step
    ycors = fish$ycor # get their y coordinates for each time step
    
    for (t in 2:time_steps) { # update likelihood for each time step for the given fish
      xlikelihood = dnorm(xcors[t], xcors[t-1]+alpha*(mu-xcors[t-1]), sigma)
      ylikelihood = dnorm(ycors[t], ycors[t-1]+alpha*(mu-ycors[t-1]), sigma)
      log_likelihood = log_likelihood + log(xlikelihood) + log(ylikelihood)
      likelihood = likelihood * xlikelihood * ylikelihood
    }
  }
  
  return (log_likelihood)
}

fish_data <- read.csv("fishdata.csv") # gets csv file of fish data 

time_steps = 100 # 100 time steps 

mu = 50 # fixed value

# alpha = 0.0250 # initial value of alpha

# sigma = 1.9 # initial value of sigma

fishes = 10 # 10 fishes

# likelihood = (1/2)^(2*fishes) # initial x and y coordinate of each fish is a bernoulli(1/2) , there are 10 fishes so 1/2^20
# log_likelihood = log(likelihood) 
# 
# for (f in 0:fishes) { # for each fish
#   
#   fish = fish_data[fish_data$who == f, ] # get only their data over 100 time steps
#   
#   xcors = fish$xcor # get their x coordinates for each time step
#   ycors = fish$ycor # get their y coordinates for each time step
#   
#   for (t in 2:time_steps) { # update likelihood for each time step for the given fish
#     xlikelihood = dnorm(xcors[t], xcors[t-1]+alpha*(mu-xcors[t-1]), sigma)
#     ylikelihood = dnorm(ycors[t], ycors[t-1]+alpha*(mu-ycors[t-1]), sigma)
#     log_likelihood = log_likelihood + log(xlikelihood) + log(ylikelihood)
#     likelihood = likelihood * xlikelihood * ylikelihood
#   }
# }


# wrap in function

alphas = seq(0, 0.5, length.out = 100)

sigmas = seq(0, 2, length.out = 100)

log_likelihoods_store = integer(length(alphas) * length(sigmas))
count = 0

for (alpha in alphas) {
  for (sigma in sigmas) {
    log_likelihoods_store[count] = log_likelihood_func(fish_data, time_steps, mu, alpha, sigma, fishes)
    count = count + 1
  } 
}

prin(max(log_likelihoods_store))



