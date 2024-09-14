rm(list = ls()) # cleans environment 

library(tmvtnorm) # multivariate truncated norm library

alpha_current = runif(1, 0, 0.5) # prior on alpha uniform in [0, 0.5]

sigma_current = runif(1, 0, 2) # prior on sigma uniform in [0, 2]

log_likelihood_func <- function(fish_data, time_steps, mu, alpha, sigma, fishes) { # log(f(xobs | alpha, sigma)) 
  
  likelihood = (1/2)^(2*fishes) # initial x and y coordinate of each fish is a bernoulli(1/2) , there are 10 fishes so 1/2^20
  log_likelihood = log(likelihood)
  

  fish = fish_data[fish_data$who == 0, ] # get only their data over 100 time steps
    
  xcors = fish$xcor # get their x coordinates for each time step
  ycors = fish$ycor # get their y coordinates for each time step
    
  for (t in 2:time_steps) { # update likelihood for each time step for the given fish
    xloglikelihood = dnorm(xcors[t], xcors[t-1]+alpha*(mu-xcors[t-1]), sigma, log=TRUE)
    yloglikelihood = dnorm(ycors[t], ycors[t-1]+alpha*(mu-ycors[t-1]), sigma, log=TRUE)
    log_likelihood = log_likelihood + xloglikelihood + yloglikelihood
    # likelihood = likelihood * xlikelihood * ylikelihood
  }
  return (log_likelihood)
}

proposal_density_sampling <- function(alpha_current, sigma_current) { # samples proposed alpha and sigma from given alpha and sigma
  mean_vector = c(alpha_current, sigma_current) # current mean of multivariate norm
  covariance_matrix = matrix(c(0.3, 0, 0, 0.3), nrow = 2) # covariance matrix of multivariate norm
  lowerbounds = c(0, 0) # norm is truncated because alpha and sigma bounded
  upperbounds = c(0.5, 2)
  proposed_vals = rtmvnorm(1, mean_vector, covariance_matrix, lowerbounds, upperbounds) # proposal values
  return (proposed_vals) 
}

proposal_density_log = function(given_alpha, given_sigma, proposal_alpha, proposal_sigma) { # q(alpha', sigma' | alpha, sigma)
  mean_vector = c(given_alpha, given_sigma) # given values as a mean 
  covariance_matrix = matrix(c(0.3, 0, 0, 0.3), nrow = 2) # covariance matrix as before
  lowerbounds = c(0, 0) # bounds as before 
  upperbounds = c(0.5, 2)
  proposed_vector = c(proposal_alpha, proposal_sigma) # density of values wanting to find 
  density_val = dtmvnorm(proposed_vector, mean_vector, covariance_matrix, lowerbounds, upperbounds, log=TRUE) # density 
  return (density_val)
}

prior_density_alpha_log = function(alpha) { # prior density of alpha uniform bounded [0, 0.5]
  return (dunif(alpha, 0, 0.5, log=TRUE))
}

prior_density_sigma_log = function(sigma) { # prior density of sigma uniform bounded [0, 2]
  return (dunif(sigma, 0, 2, log=TRUE))
}

iterations = 10000 # number of iterations 

fish_data <- read.csv("fishdata.csv") # gets csv file of fish data 

time_steps = 100 # 100 time steps 

mu = 50 # fixed value

fishes = 1 # 10 fishes

accepted_alphas = list() # list of accepted alphas
accepted_sigmas = list() # list of accepted sigmas 

for (t in 1:iterations) { # for each iteration
  proposed_vals = proposal_density_sampling(alpha_current, sigma_current) # propose values from current values
  proposal_alpha = proposed_vals[1] # store
  proposal_sigma = proposed_vals[2]

  # print("-----------------------------------")
  # print((prior_density_alpha(proposal_alpha)))
  # print(log(prior_density_sigma(proposal_sigma)))
  # print(log_likelihood_func(fish_data, time_steps, mu, proposal_alpha, proposal_sigma, fishes))
  # print(log(proposal_density(proposal_alpha, proposal_sigma, alpha_current, sigma_current)) )

  numerator = prior_density_alpha_log(proposal_alpha) +
              prior_density_sigma_log(proposal_sigma) +
              log_likelihood_func(fish_data, time_steps, mu, proposal_alpha, proposal_sigma, fishes) +
              proposal_density_log(proposal_alpha, proposal_sigma, alpha_current, sigma_current)

  # print(log(prior_density_alpha(alpha_current)))
  # print(log(prior_density_sigma(sigma_current)))
  # print(log_likelihood_func(fish_data, time_steps, mu, alpha_current, sigma_current, fishes))
  # print(log(proposal_density(alpha_current, sigma_current, proposal_alpha, proposal_sigma)))

  denominator = prior_density_alpha_log(alpha_current) +
                prior_density_sigma_log(sigma_current) +
                log_likelihood_func(fish_data, time_steps, mu, alpha_current, sigma_current, fishes) +
                proposal_density_log(alpha_current, sigma_current, proposal_alpha, proposal_sigma)

  log_acceptance_prob = numerator - denominator # log acceptance

  r = runif(1, 0, 1)

  if (log_acceptance_prob >= log(r)) { # if log acceptance is high
    accepted_alphas = c(accepted_alphas, proposal_alpha) # add alpha to samples
    accepted_sigmas = c(accepted_sigmas, proposal_sigma) # add sigma to samples
    alpha_current = proposal_alpha # update current alpha
    sigma_current = proposal_sigma # update current sigma
  }
  
  else {
    accepted_alphas = c(accepted_alphas, alpha_current)
    accepted_sigmas = c(accepted_sigmas, sigma_current)
    }
}

alphas = list(mean=mean(unlist(accepted_alphas)), std=sd(unlist(accepted_alphas)))

sigmas = list(mean=mean(unlist(accepted_sigmas)), std=sd(unlist(accepted_sigmas)))

sequence_samples = seq(0, 0.5, length=1000)
hist(unlist(accepted_alphas), freq=FALSE, xlab="Alphas", main="Histogram of samples")

sequence_samples = seq(0, 2, length=1000)
hist(unlist(accepted_sigmas), freq=FALSE, xlab="Sigmas", main="Histogram of samples")

