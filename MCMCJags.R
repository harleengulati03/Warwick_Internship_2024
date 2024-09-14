

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

library(rjags) # call library

# define the model
modelstring = "
model {
  p ~ dunif(0,1)
  y ~ dbinom(p, n)

}" 

bin_dta = list("y" = yobs, "n" = n) # data observed

model = jags.model(textConnection(modelstring), data=bin_dta) # initialize model

update(model, n.iter=1000) # burn-in for 1000 iterations

samples = coda.samples(model, variable.names = c("p"), n.iter =10000) # 10000 iterations of samples collected

traceplot(samples) # diagnostics
output_matrix = as.matrix(samples)
acf(output_matrix) # diagnostics

# from diagnostics find thinning needed, repeat with thinning
modelstring = "
model {
  p ~ dunif(0,1)
  y ~ dbinom(p, n)

}"

bin_dta = list("y" = yobs, "n" = n)

model = jags.model(textConnection(modelstring), data=bin_dta)

update(model, n.iter=1000)

samples = coda.samples(model, variable.names = c("p"), n.iter =10000, thin=5)

traceplot(samples)
output_matrix = as.matrix(samples)
acf(output_matrix)

summary(samples) # summary of samples

# graph plotting
sequence_samples = seq(0, 1, length=1000)
hist(unlist(samples), freq=FALSE, xlab="Samples", main="Histogram of samples")
lines(sequence_samples, dbeta(sequence_samples, posterior_alpha, posterior_beta), type ='l')
