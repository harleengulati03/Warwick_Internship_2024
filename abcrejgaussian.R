rm(list = ls()) # cleans environment 

p = 0.5 
n = 10
prior_alpha = 1
prior_beta = 1

set.seed(2)
cointosses = rbinom(n, 1, p)
yobs = sum(cointosses == 1) # number of successes = 5 currently 
set.seed(NULL)

posterior_alpha = prior_alpha + yobs # 6
posterior_beta = prior_beta + n - yobs # 6

posterior_mean = posterior_alpha / (posterior_alpha + posterior_beta) # 0.5
denominator_intermediate_step = (posterior_alpha + posterior_beta)^2 * (posterior_alpha + posterior_beta + 1)
posterior_std = sqrt((posterior_alpha * posterior_beta) / denominator_intermediate_step) # 0.139

N = 10000 # number of iterations 

abc = function(h, yobs) {
  
  K = 1 / (h * sqrt(2 * pi)) # as K>= Kh(0) * prior/g = Kh(0) = 1/h * 1/ sqrt(2pi) * e^(-1/2 *0^2) = 1/(h*sqrt(2pi))
  abc_gaussianker = list() # samples accepted 
  
  gaussian_ker = function(y, yobs, h) {
    u = abs(y-yobs)
    return (1/(h * sqrt(2 * pi)) * exp( (-1/2) * u^2))
  }
  
  for (i in 1:N) {
    proposed_p = runif(1, 0, 1) # proposal density 
    proposed_sucs = rbinom(1, n, proposed_p) # how many successes are observed under proposed_p
    acceptprob = (gaussian_ker(proposed_sucs, yobs, h) * dunif(proposed_p, 0, 1)) / (K * dunif(proposed_p, 0, 1)) # for now h = 1
    randomunif = runif(1, 0, 1)
    if (randomunif <= acceptprob) {
      abc_gaussianker <- c(abc_gaussianker, list(proposed_p))
    }
  }
  
  return (abc_gaussianker)
} 

abcgaussiankerh1 = abc(1, yobs)
abcgaussiankerh2 = abc(2, yobs)
abcgaussiankerh3 = abc(3, yobs)

summarygaussh1 = list(mean = mean(unlist(abcgaussiankerh1)), std=sd(unlist(abcgaussiankerh1)))

summarygaussh2 = list(mean = mean(unlist(abcgaussiankerh2)), std=sd(unlist(abcgaussiankerh2)))

summarygaussh3 = list(mean = mean(unlist(abcgaussiankerh3)), std=sd(unlist(abcgaussiankerh3)))

par(mfrow=c(1,3))

sequence_samples = seq(0, 1, length=1000)

hist(unlist(abcgaussiankerh1), freq=FALSE, xlab="Samples", main="ABC Gaussian Kernel , h = 1")
lines(sequence_samples, dbeta(sequence_samples, posterior_alpha, posterior_beta), type ='l')

hist(unlist(abcgaussiankerh2), freq=FALSE, xlab="Samples", main="ABC Gaussian Kernel , h = 2")
lines(sequence_samples, dbeta(sequence_samples, posterior_alpha, posterior_beta), type ='l')

hist(unlist(abcgaussiankerh3), freq=FALSE, xlab="Samples", main="ABC Gaussian Kernel , h = 3")
lines(sequence_samples, dbeta(sequence_samples, posterior_alpha, posterior_beta), type ='l')





