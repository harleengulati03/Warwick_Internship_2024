rm(list = ls()) 
library(mcmcse)
library(tmvtnorm) 

proposal_density_sampling <- function(alpha_current, sigma_current) { # samples proposed alpha and sigma from given alpha and sigma
  mean_vector = c(alpha_current, sigma_current) # current mean of multivariate norm
  covariance_matrix = matrix(c((2.38/sqrt(2))*0.001898841, 0, 0, (2.38/sqrt(2))*0.01819092), nrow = 2) # covariance matrix of multivariate norm
  lowerbounds = c(0, 0) # norm is truncated because alpha and sigma bounded
  upperbounds = c(0.5, 2)
  proposed_vals = rtmvnorm(1, mean_vector, covariance_matrix, lowerbounds, upperbounds) # proposal values
  return (proposed_vals) 
} # simulated alpha and sigma 

proposal_density_log = function(given_alpha, given_sigma, proposal_alpha, proposal_sigma) { # q(alpha', sigma' | alpha, sigma)
  mean_vector = c(given_alpha, given_sigma) # given values as a mean 
  covariance_matrix = matrix(c((2.38/sqrt(2))*0.001898841, 0, 0, (2.38/sqrt(2))*0.01819092), nrow = 2) # covariance matrix of multivariate norm
  lowerbounds = c(0, 0) # bounds as before 
  upperbounds = c(0.5, 2)
  proposed_vector = c(proposal_alpha, proposal_sigma) # density of values wanting to find 
  density_val = dtmvnorm(proposed_vector, mean_vector, covariance_matrix, lowerbounds, upperbounds, log = TRUE) # density 
  return (density_val)
} # proposal density log 

prior_density_alpha_log = function(alpha) { # prior density of alpha uniform bounded [0, 0.5]
  return (dunif(alpha, 0, 0.5, log=TRUE))
} # prior for alpha density 

prior_density_sigma_log = function(sigma) { # prior density of sigma uniform bounded [0, 2]
  return (dunif(sigma, 0, 2, log=TRUE))
} # prior for sigma density

euclidean_distance = function(x, y) {sqrt(x^2 + y^2)} # euclidean distance for (x,y) from (0,0)

euclidean_distances = function(fish_cor) {
  distances = matrix(0, nrow=101, ncol=1)
  for (count in 1:nrow(fish_cor)) {
    distances[count, 1] = euclidean_distance(fish_cor[count, 1], fish_cor[count, 2])
  }
  return (distances)
} # for a fish, vector of euclidean distances at each time step

distance_measure = function(observed_euclidean_dists, simulated_dists) {
  closeness_measure = abs(sum(observed_euclidean_dists - simulated_dists))
  return (closeness_measure)
  # return (1/closeness_measure)
} # for [x,y,z] and [xobs, yobs, zobs] we return ((xobs-x) + (yobs-y) + (zobs-z))^-1 for x,y,z,xobs,yobs,zobs being euclidean distances of each time step 

# kernel_function_log = function(u, h) {
#   # z = 1/(h*sqrt(2 * pi))
#   return ((-1/2 * u^2))
# } # gaussian kernel with u being measure of closeness 

# kernel_function_triangular = function(u, h) {
#   if (u <= h) {return (0)}
#   else {return (log (1/h) + log(h-u))}
# }

kernel_func_uniform = function(u, h) {
  if (u <= h) {return (1)}
  else {return (0)}
}

kernel_finding_function_unif = function(simulated_positions, observed_dists) {
  simulated_dists = euclidean_distances(simulated_positions)
  print(distance_measure(observed_dists, simulated_dists))
  kernel_val = kernel_function_uniform(distance_measure(observed_dists, simulated_dists), 0.0005)
  if (kernel_val == 1) {return (0)}
  else {return (-1)}
} # call to other functions to find kernel 

# kernel_finding_function = function(simulated_positions, observed_dists) {
#   simulated_dists = euclidean_distances(simulated_positions)
#   print(distance_measure(observed_dists, simulated_dists))
#   return (kernel_function_uniform(distance_measure(observed_dists, simulated_dists), 1))
# } # call to other functions to find kernel 

simulate_ou_step = function(x,y,alpha,sigma,mu)
{
  new_x = x + alpha*(mu - x) + rnorm(1,0,sigma)
  new_y = y + alpha*(mu - y) + rnorm(1,0,sigma)
  return(list(new_x,new_y))
}

multiple_simulate_ou_step = function(x,y,alpha,sigma,mu)
{
  new_x = x
  new_y = y
  for (i in 1:length(x)) # for each turtles x and y position
  {
    new_position = simulate_ou_step(x[i],y[i],alpha,sigma,mu)
    new_x[i] = new_position[[1]]
    new_y[i] = new_position[[2]]
  }
  return(list(new_x,new_y))
}

simulate_turtles = function(alpha,sigma,mu,T,n)
{
  positions = matrix(0,(T+1)*n,2) 
  for (i in 1:n) # for each turtle , set their initial position to either 0 or 1
  {
    positions[i,1] = round(runif(1)) 
    positions[i,2] = round(runif(1)) 
  }
  
  for (t in 1:T) # for each time-step 
  { 
    new_positions = multiple_simulate_ou_step(positions[(t-1)*n+(1:n),1],positions[(t-1)*n+(1:n),2],alpha,sigma,mu)
    positions[t*n+(1:n),1] = new_positions[[1]] 
    positions[t*n+(1:n),2] = new_positions[[2]]
  }
  
  return(positions)
}

## ------------------------------------------------------------------------- ##

fish_data <- read.csv("fishdata.csv") # gets csv file of fish data 
time_steps = 100 # 100 time steps 
mu = 50 # fixed value
fishes = 1 # 10 fishes
fish = fish_data[fish_data$who == 0, ] # data for one fish 
fish_cor = data.frame("xcoordinate" = fish$xcor, "ycoordinate" = fish$ycor) # observed data
observed_distances = euclidean_distances(fish_cor) # euclidean distances of each time point in observed data
iterations = 100 # number of iterations 
burn_in = 0
alpha_current = runif(1, 0, 0.5) # prior on alpha uniform in [0, 0.5] 0.0279
sigma_current = runif(1, 0, 2) # prior on sigma uniform in [0, 2]
simulated_positions_current = simulate_turtles(alpha_current, sigma_current, mu, time_steps, 1) # simulated distances given alpha and sigma

while (kernel_finding_function_unif(simulated_positions_current, observed_distances) == -1) {
  alpha_current = runif(1, 0, 0.5) # prior on alpha uniform in [0, 0.5]
  sigma_current = runif(1, 0, 2) # prior on sigma uniform in [0, 2]
  simulated_positions_current = simulate_turtles(alpha_current, sigma_current, mu, time_steps, 1)
}
accepted_vals = matrix(0, nrow=iterations, ncol=2) 

## ------------------------------------------------------------------------- ##

# proposed_vals = proposal_density_sampling(alpha_current, sigma_current) # gets alpha and sigma proposal
# proposal_alpha = proposed_vals[1] 
# proposal_sigma = proposed_vals[2]
# 
# simulated_positions_new = simulate_turtles(proposal_alpha, proposal_sigma, mu, time_steps, 1) # corresponding simulated data
# kernel_log_new = kernel_finding_function(simulated_positions_new, observed_distances) # kernel of new simulated points
# kernel_log_old = kernel_finding_function(simulated_positions_current, observed_distances) # kernel of old simulated points

## ------------------------------------------------------------------------- ##
for (t in 1:iterations) { # for each iteration
  proposed_vals = proposal_density_sampling(alpha_current, sigma_current) # propose values from current values
  proposal_alpha = proposed_vals[1] 
  proposal_sigma = proposed_vals[2]
  simulated_positions_new = simulate_turtles(proposal_alpha, proposal_sigma, mu, time_steps, 1)
  
  kernel_log_new = kernel_finding_function_unif(simulated_positions_new, observed_distances)
  
  kernel_log_old = kernel_finding_function_unif(simulated_positions_current, observed_distances)
  
  if(kernel_log_new == -1) {
    log_acceptance_prob = -1000000000000000000000
  }
  
  else {
    numerator = prior_density_alpha_log(proposal_alpha) +
      prior_density_sigma_log(proposal_sigma) +
      proposal_density_log(proposal_alpha, proposal_sigma, alpha_current, sigma_current) +
      kernel_log_new
    
    denominator = prior_density_alpha_log(alpha_current) +
      prior_density_sigma_log(sigma_current) +
      proposal_density_log(alpha_current, sigma_current, proposal_alpha, proposal_sigma) +
      kernel_log_old
    
    log_acceptance_prob = numerator - denominator # log acceptance
  }
  
  r = runif(1, 0, 1)
  
  if (log_acceptance_prob >= log(r)) { # if log acceptance is high
    accepted_vals[t, 1] = proposal_alpha
    accepted_vals[t, 2] = proposal_sigma
    # accepted_alphas = c(accepted_alphas, proposal_alpha) # add alpha to samples
    # accepted_sigmas = c(accepted_sigmas, proposal_sigma) # add sigma to samples
    alpha_current = proposal_alpha # update current alpha
    sigma_current = proposal_sigma # update current sigma
    simulated_positions_current = simulated_positions_new
  }
  
  else {
    accepted_vals[t, 1] = alpha_current
    accepted_vals[t, 2] = sigma_current
  }
}

alphas_efficent = (accepted_vals[, 1])[burn_in+1:iterations]
sigmas_efficent = (accepted_vals[, 2])[burn_in+1:iterations]

alphas = list(mean=mean(unlist(alphas_efficent)), std=sd(unlist(alphas_efficent)))
sigmas = list(mean=mean(unlist(sigmas_efficent)), std=sd(unlist(sigmas_efficent)))

hist(unlist(alphas_efficent), freq=FALSE, xlab="Alphas", main="Histogram of samples")

hist(unlist(sigmas_efficent), freq=FALSE, xlab="Sigmas", main="Histogram of samples")


plot(as.matrix(unlist(alphas_efficent)))

plot(as.matrix(unlist(sigmas_efficent)))

multiESS(as.matrix(unlist(alphas_efficent)))

multiESS(as.matrix(unlist(sigmas_efficent)))
