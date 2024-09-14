rm(list = ls()) # cleans environment 

p = 0.5
n = 10
prior_alpha = 1
prior_beta = 1

set.seed(2)
cointosses = rbinom(n, 1, p)
yobs = sum(cointosses == 1) # number of successes = 5 currently

posterior_alpha = prior_alpha + yobs # 6 
posterior_beta = prior_beta + n - yobs # 6 

posterior_mean = posterior_alpha / (posterior_alpha + posterior_beta) # 0.5
denominator_intermediate_step = (posterior_alpha + posterior_beta)^2 * (posterior_alpha + posterior_beta + 1)
posterior_std = sqrt((posterior_alpha * posterior_beta) / denominator_intermediate_step) # 0.139

library(truncnorm)

abc_mcmc = function(h, yobs, n) {
  
  proposal_density = function(x, mean, std) { # gets values from the proposal density 
  return (dtruncnorm(x, a =0, b=1, mean=mean, sd=std))
}

  proposal_density_sampling = function(mean, std) { # gets random samples from the proposal density 
    return (rtruncnorm(1, a =0, b=1, mean=mean, sd=std))
  }
  
  N = 10000 # number of iterations 
  var_chosen = 0.01 # some value for variance chosen
  abc_mcmc_samples = list() # samples accepted 
  
  uniform_kernel = function(h, y, yobs) {
    if (abs(y-yobs) <= h) {return (1)}
    else {return (0)}
  }
  
  initial_p = runif(1, 0, 1) # initial p
  
  initial_y = sum(rbinom(n, 1, initial_p) == 1) # number of successes from initial p
  
  while (uniform_kernel(h, initial_y, yobs) == 0) { # if initial y not within bandwidth of h , then re-sample a new p and a new y
    initial_p = runif(1, 0, 1)
    
    initial_y = sum(rbinom(n, 1, initial_p) == 1)
  } # continue doing this until get an initial p and initial y which is within bandwidth h
  
  current_p = initial_p # once out of while loop, set the initial values to be the current ones
  
  current_y = initial_y 
  
  for (i in 1:N) {
    
    proposed_p = proposal_density_sampling(current_p, sqrt(var_chosen)) # gets a proposed p from a truncated normal mean = current p, var = 0.01
    
    proposed_y = sum(rbinom(n, 1, proposed_p) == 1) # gets a proposed y from the proposed p
    
    old_given_new = proposal_density(current_p, proposed_p, sqrt(var_chosen))
    
    new_given_old = proposal_density(proposed_p, current_p, sqrt(var_chosen))
    
    acceptanceprob = (1/(2*h) *uniform_kernel(h, proposed_y, yobs) * dunif(proposed_p, 0, 1) * old_given_new) / 
      (1/(2*h) * uniform_kernel(h, current_y, yobs) * dunif(current_p, 0, 1) * new_given_old)
    
    acceptanceprob = min(1, acceptanceprob)
    
    if ((acceptanceprob >= runif(1,0,1)) && (i <= 100)) {
      current_p = proposed_p # the proposed values becomes the current values 
      current_y = proposed_y
    }
    
    else if (acceptanceprob >= runif(1,0,1)) {
      abc_mcmc_samples = c(abc_mcmc_samples, list(proposed_p))
      current_p = proposed_p # the proposed values becomes the current values 
      current_y = proposed_y
    }
    
    else if (i > 100) {abc_mcmc_samples = c(abc_mcmc_samples, list(current_p))}
    
    # else {abc_mcmc_samples = c(abc_mcmc_samples, list(current_p))} # add the current sample to the list of chosen samples 
    
  }
  
  return (abc_mcmc_samples)
  
}

abc_mcmc_h_1 = abc_mcmc(1, yobs, n)

abc_mcmc_h_2 = abc_mcmc(1, yobs, n)

abc_mcmc_h_3 = abc_mcmc(3, yobs, n)

abc_mcmc_h_1_summary = list(mean=mean(unlist(abc_mcmc_h_1)), std = sqrt(var(unlist(abc_mcmc_h_1))))

abc_mcmc_h_2_summary = list(mean=mean(unlist(abc_mcmc_h_2)), std = sqrt(var(unlist(abc_mcmc_h_2))))

abc_mcmc_h_3_summary = list(mean=mean(unlist(abc_mcmc_h_3)), std = sqrt(var(unlist(abc_mcmc_h_3))))

par(mfrow=c(1,3))

sequence_samples = seq(0, 1, length=1000)

hist(unlist(abc_mcmc_h_1), freq=FALSE, xlab="Samples", main="ABC MCMC")
lines(sequence_samples, dbeta(sequence_samples, posterior_alpha, posterior_beta), type ='l')

hist(unlist(abc_mcmc_h_2), freq=FALSE, xlab="Samples", main="ABC MCMC")
lines(sequence_samples, dbeta(sequence_samples, posterior_alpha, posterior_beta), type ='l')

hist(unlist(abc_mcmc_h_3), freq=FALSE, xlab="Samples", main="ABC MCMC")
lines(sequence_samples, dbeta(sequence_samples, posterior_alpha, posterior_beta), type ='l')



