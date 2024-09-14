rm(list = ls()) 
library(mcmcse)
library(tmvtnorm) 
library(truncnorm)

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

summary_stats_sigma = function(data) {
  last_time_step = data[(nrow(data)-9):nrow(data), ]
  all_values <- as.vector(last_time_step)
  std <- sd(all_values)
  return (std)
}

prior_density_log = function(sigma) {
  return (dunif(sigma, 0, 2, log=TRUE))
}

proposal_density_sampling <- function(sigma_current) {
  proposed_val = rtruncnorm(1, a =0, b=2, mean=sigma_current, sd=(2.38/sqrt(2))*0.01819092)
  return (proposed_val) 
} # simulated sigma 

proposal_density_log = function(given_sigma, proposal_sigma) { # q(alpha', sigma' | alpha, sigma)
  density_val = log(dtruncnorm(proposal_sigma, a = 0, b= 2, mean=given_sigma, sd=(2.38/sqrt(2))*0.01819092))
  return (density_val)
} # proposal density log

kernel_gaussian_log = function(s, sobs) {
 u = abs(s - sobs)
 return (-1/2 * u^2)
}

iterations = 300000
mu = 50
time_steps = 100
fishes = 10
true_alpha = 0.030
sigma_current = runif(1, 0, 2)
simulated_data_current = simulate_turtles(true_alpha, sigma_current, mu, time_steps, fishes)
summary_stat_sigma_current = summary_stats_sigma(simulated_data_current)
fish_data <- read.csv("fishdata.csv") # gets csv file of fish data 
xcors = fish_data$xcor
ycors = fish_data$ycor
observed_data = cbind(xcors, ycors)
sigmaobs = summary_stats_sigma(observed_data)
accepted_vals = matrix(0, nrow=iterations, ncol=2) 

for (t in 1:iterations) {
  sigma_proposed = proposal_density_sampling(sigma_current)
  simulated_data_proposed = simulate_turtles(true_alpha, sigma_proposed, mu, time_steps, fishes)
  summary_stat_sigma_proposed = summary_stats_sigma(simulated_data_proposed)
  
  numerator = prior_density_log(sigma_proposed) + kernel_gaussian_log(summary_stat_sigma_proposed, sigmaobs) + proposal_density_log(sigma_proposed, sigma_current)
  denominator = prior_density_log(sigma_current) + kernel_gaussian_log(summary_stat_sigma_current, sigmaobs) + proposal_density_log(sigma_current, sigma_proposed)
  
  logaccept = numerator - denominator 
  
  r = runif(1, 0, 1)
  
  if (logaccept >= log(r)) {
    sigma_current = sigma_proposed 
    simulated_data_current = simulated_data_proposed
    summary_stat_sigma_current = summary_stat_sigma_proposed
    accepted_vals[t, 1] = sigma_proposed
  }
  
  else {accepted_vals[t, 1] = sigma_current}
}

sigmas_efficent = (accepted_vals[, 1])[1:iterations]

sigmas = list(mean=mean(unlist(sigmas_efficent)), std=sd(unlist(sigmas_efficent)))

hist(unlist(sigmas_efficent), freq=FALSE, xlab="Sigmas", main="Histogram of samples")




