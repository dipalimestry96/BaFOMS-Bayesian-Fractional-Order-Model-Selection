## R code for simulation study 
## We have total six fractional logistic models under six different fractional derivative (FD) definitions
##   M1: Under Caputo FD definition
##   M2: Under Caputo-Fabrizio FD definition
##   M3: Under ABC FD definition
##   M4: Under Riemann-Liouville FD definition
##   M5: Under CF-RL FD definition
##   M6: Under Grunwal-Letnikov definition

## Simulated data: we have simulated data from the model M1 and then using this
## dataset we have performed the model selection


## sequence of the code is as follows:
#### 1. Simulated data Generation
#### 2. model fitting and getting posterior distribution of model parameters using JAGS
#### 3. final step: Multimodel inference using RJMCMC
#### 4. visualization of outputs by plotting them

##################################################################################
##################################################################################
#### 1. Simulated data Generation
set.seed(1222)  # for reproducibility of results
setwd("specify path of local directory here") # path to save outputs

## data simulation step from the model M1: Caputo logistic model
K_sim = 30000         # Carrying capacity
alpha_sim = 0.25       # order of FDE
r_sim = 0.9           # growth rate
sigma_sim = 300       # standard deviation
tau_sim = 1/sigma_sim^2   # precision = 1/ variance
n = 40     # length of data to be generated
t0 = 0     # initial time point
tn = 20    # end time point

t = seq(t0, tn, length.out = n)    # time vector
h = t[2]-t[1]                      # time step

M_alpha = 1                      # normalising function
ABC_alpha = 1                    # normalising function


# mean function for model M1
mean_M1 = function(N_p, N_pp, t_c, t_p, t_pp, h, alpha, r, K){          
  A = ((t_c^(alpha+1))/(alpha*(alpha+1)))- (((t_c^alpha)*t_pp)/alpha)-
    ((t_p^(alpha+1))/(alpha*(alpha+1))) + (((t_p^alpha)*t_pp)/alpha)
  B = ((t_c^(alpha+1))/(alpha*(alpha+1))) - (((t_c^alpha)*t_p)/alpha) + 
    ((t_p^(alpha+1))/(alpha+1))
  C = N_p + ((r*N_p*(1-(N_p/K)))/(h*gamma(alpha)))*A - 
    (((r*N_pp*(1-(N_pp/K))))/(h*gamma(alpha)))*B
  return(C)
} 

N = numeric(n)   # array to store the data
mu_0 = 5000      # mean for initial population size
N[1] = rnorm(1,mean = mu_0, sd = sigma_sim)
mu_star = N[1] + h*r_sim*N[1]*(1-(N[1]/K_sim))
N[2] = rnorm(1, mean = mu_star, sd = sigma_sim )
for (i in 3:n) {
  mean = mean_M1(N[i-1], N[i-2], t[i], t[i-1],t[i-2], h, alpha_sim, r_sim, K_sim)
  N[i] = rnorm(1, mean = mean, sd = sigma_sim)
}

par(mfrow = c(1,1))
plot(t, N, type = "b", col= "red", pch = 16, main = "simulated data from M1")


### this simulated data we treat as a true data

## for model selection, we consider the RJMCMC approach which uses 
## the posterior model probabilities of models.
## For that we require the posterior distribution of model parameters of each model
## for getting the posterior distribution we are plane to use the jags() function from the package R2jags

################################################################################
################################################################################
#### 2. model fitting and getting posterior distribution of model parameters using JAGS

library(R2jags)
library(coda)
library(madness)

## theta = c(alpha, r, K, tau, mu_0) where tau = 1/sigma^2
## prior specification of all model parameters

# alpha ~ unif(0,1)
# r ~ gamma(2,1/2)
# K ~ normal(30000, 10000)
# tau ~ gamma(0.01, 0.01)
# mu_0 ~ norm(5000, 1000)


# log-prior for all model parameters
log_prior = function(theta){
  sum(dunif(theta[1], min = 0, max = 1, log = TRUE),
      dgamma(theta[2], 2, 1/2, log = TRUE) ,
      dnorm(theta[3], mean = 30000, sd = 100, log = TRUE),
      dgamma(theta[4], 0.01, 0.01, log = TRUE ),
      dnorm(theta[5], 5000, sqrt(1000), log = TRUE))
}

## common initial values of parameters for all the models fro passing to JAGS
inits = function(){
  list(alpha = 0.7, r = 0.7, K = 1000, tau = 2, mu_0 = 1000)
}

params = c("alpha", "r", "K", "tau", "mu_0")   # parameters
 
## bijection and its inverse for RJMCMC
g = function(psi){ psi }
ginv = function(theta){ theta }


n_iter = 2e5  # no of iteration for jags

#################################################################################################################

## for model M1: Caputo logistic model

## mean function

mean_M1 = function(N_p, N_pp, t_c, t_p, t_pp, h, alpha, r, K){         
  A = ((t_c^(alpha+1))/(alpha*(alpha+1)))- (((t_c^alpha)*t_pp)/alpha)-
    ((t_p^(alpha+1))/(alpha*(alpha+1))) + (((t_p^alpha)*t_pp)/alpha)
  B = ((t_c^(alpha+1))/(alpha*(alpha+1))) - (((t_c^alpha)*t_p)/alpha) + 
    ((t_p^(alpha+1))/(alpha+1))
  C = N_p + ((r*N_p*(1-(N_p/K)))/(h*gamma(alpha)))*A - 
    (((r*N_pp*(1-(N_pp/K))))/(h*gamma(alpha)))*B
  return(C)
}

# log_likelihood for M1
log_L_M1 = function(theta){
  s1 = dnorm(N[1], mean = theta[5], sd = sigma_sim, log = TRUE)
  s2 = dnorm(N[2],mean = N[1]+ h*r_sim*N[1]*(1-(N[1]/K_sim)),sd = sigma_sim, log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M1(N[i], N[i-1],t[i+1], t[i],t[i-1],h,theta[1], 
                                           theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE) 
  }
  z = s1+ s2 + s3
  return(z)
}


# model txt file of M1
cat("model{
  N[1] ~ dnorm(mu_0,1/sigma_sim^2)
  N[2] ~ dnorm(N[1] + h*r_sim*N[1]*(1-N[1]/K_sim), 1/sigma_sim^2)
  for (i in 3:n) {
       mu[i] = N[i-1] + ((r*N[i-1]*(1-(N[i-1]/K)))/(h*exp(loggam(alpha))))*(((t[i]^(alpha+1))/(alpha*(alpha+1)))-
                (((t[i]^alpha)*t[i-2])/alpha) - ((t[i-1]^(alpha+1))/(alpha*(alpha+1))) + (((t[i-1]^alpha)*t[i-2])/alpha) ) -
                ((r*N[i-2]*(1-(N[i-2]/K)))/(h*exp(loggam(alpha))))*( ((t[i]^(alpha+1))/(alpha*(alpha+1))) -
                 (((t[i]^alpha)*t[i-1])/alpha) +((t[i-1]^(alpha+1))/(alpha+1)) )
  }
   for (j in 3:n){    
      N[j] ~ dnorm(mu[j], tau)
  }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dnorm(30000, 1/10000)
  tau ~ dgamma(0.01, 0.01)
  mu_0 ~ dnorm(5000, 1/1000)
 }", file = "M1_Caputo.txt")


## Posterior distribution using JAGS
jagsfit_M1 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim'), inits, params,
                  n.iter = n_iter, model.file = "M1_Caputo.txt", n.chains = 3, n.thin = 1)
fit_M1 = as.mcmc(jagsfit_M1) # mcmc after burnin period=iter/2

# checking posterior distribution convergence 
gelman.diag(fit_M1)
gelman.plot(fit_M1)
traceplot(fit_M1)

# storing the posterior distribution of parameters of model M1
C1 = as.matrix(fit_M1)
write.table(C1, "M1_output.txt", row.names = FALSE)

## function for sampling from posterior distribution
draw1 = function(){
  C1[sample(dim(C1)[1], 1, replace = TRUE), -which(colnames(C1) == "deviance")][c(1,4,2,5,3)]
}


########################################################################################################################

## for model M2: Caputo-Fabrizio logistic model

# mean function
mean_M2 = function(N_p, N_pp, M_alpha, h, alpha, r, K){
  A = ((1-alpha)/M_alpha + (3*alpha*h)/(2*M_alpha))*r*N_p*(1-(N_p/K))
  B = ((1-alpha)/M_alpha + (alpha*h)/(2*M_alpha))*r*N_pp*(1-(N_pp/K))
  return(N_p + A - B)
}

# log_likelihood for M2
log_L_M2 = function(theta){
  s1 = dnorm(N[1], mean = theta[5], sd = sigma_sim, log = TRUE)
  s2 = dnorm(N[2],mean = N[1]+ h*r_sim*N[1]*(1-(N[1]/K_sim)),sd = sigma_sim, log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M2(N[i], N[i-1],M_alpha,h,theta[1], theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE) 
  }
  z = s1+ s2 + s3
  return(z)
}

# model txt file of M2
cat("model{
  N[1] ~ dnorm(mu_0,1/sigma_sim^2)
  N[2] ~ dnorm(N[1] + h*r_sim*N[1]*(1-N[1]/K_sim), 1/sigma_sim^2)
  for(i in 3:n){
    mu[i] = N[i-1] + ((1-alpha)/M_alpha + (3*alpha*h)/(2*M_alpha))*r*N[i-1]*(1-(N[i-1]/K)) 
               - ((1-alpha)/M_alpha + (alpha*h)/(2*M_alpha))*r*N[i-2]*(1-(N[i-2]/K))
    }
   for (j in 3:n){    
      N[j] ~ dnorm(mu[j], tau)
  }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dnorm(30000, 1/10000)
  tau ~ dgamma(0.01, 0.01)
  mu_0 ~ dnorm(5000, 1/1000)
 }", file = "M2_Caputo_Fabrizio.txt")


## Posterior distribution using JAGS
jagsfit_M2 = jags(data = c('N', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim','M_alpha'), inits, params,
                  n.iter = n_iter, model.file = "M2_Caputo_Fabrizio.txt", n.chains = 3, n.thin = 1)
fit_M2 = as.mcmc(jagsfit_M2)
class(fit_M2)

# checking posterior distribution convergence
gelman.diag(fit_M2)
gelman.plot(fit_M2)
traceplot(fit_M2)

# storing the posterior distribution of parameters of model M2
C2 = as.matrix(fit_M2)
write.table(C2, "M2_output.txt", row.names = FALSE)

## function for sampling from posterior distribution
draw2 = function(){
  C2[sample(dim(C2)[1], 1, replace = TRUE), -which(colnames(C2) == "deviance")][c(1,4,2,5,3)]
}

#############################################################################################################################

## for model M3: ABC logistic model

## mean function
mean_M3 = function(N_p, N_pp, t_c, t_p,t_pp, h, ABC_alpha, alpha, r, K){
  A = ((1-alpha)/ABC_alpha)*(r*N_p*(1-(N_p/K)) - r*N_pp*(1-(N_pp/K)))
  B = (t_c^(alpha+1)/(alpha*(alpha+1)))-(((t_c^alpha)*t_pp)/alpha)
  -(t_p^(alpha+1)/(alpha*(alpha+1))) + (((t_p^alpha)*t_pp)/alpha)
  C = (t_c^(alpha+1)/(alpha*(alpha+1)))- (((t_c^alpha)*t_p)/alpha) 
  + (t_p^(alpha+1)/(alpha+1))
  D = N_p + A + ((alpha*r*N_p*(1-(N_p/K)))/(h*ABC_alpha*gamma(alpha)))*B 
  - ((alpha*r*N_pp*(1-(N_pp/K)))/(h*ABC_alpha*gamma(alpha)))*C
  return(D)
}

##log-likelihood for M3
log_L_M3 = function(theta){
  s1 = dnorm(N[1], mean = theta[5], sd = sigma_sim, log = TRUE)
  s2 = dnorm(N[2],mean = N[1]+ h*r_sim*N[1]*(1-(N[1]/K_sim)),sd = sigma_sim, log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M3(N[i], N[i-1],t[i+1],t[i],t[i-1],h,ABC_alpha,theta[1], theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE) 
  }
  z = s1+ s2 + s3
  return(z)
}


# model txt file of M3
cat("model{
  N[1] ~ dnorm(mu_0,1/sigma_sim^2)
  N[2] ~ dnorm(N[1] + h*r_sim*N[1]*(1-N[1]/K_sim), 1/sigma_sim^2)
  for(i in 3:n){
    mu[i] = N[i-1] + ((1-alpha)/ABC_alpha)*(r*N[i-1]*(1-(N[i-1]/K)) - r*N[i-2]*(1-(N[i-2]/K))) + 
            ((alpha*r*N[i-1]*(1-(N[i-1]/K)))/(ABC_alpha*h*exp(loggam(alpha))))*((t[i]^(alpha+1)/(alpha*(alpha+1)))-
            (((t[i]^alpha)*t[i-2])/alpha) -(t[i-1]^(alpha+1)/(alpha*(alpha+1))) + (((t[i-1]^alpha)*t[i-2])/alpha)) -
            ((alpha*r*N[i-2]*(1-(N[i-2]/K)))/(ABC_alpha*h*exp(loggam(alpha))))*((t[i]^(alpha+1)/(alpha*(alpha+1)))- 
            (((t[i]^alpha)*t[i-1])/alpha)  + (t[i-1]^(alpha+1)/(alpha+1)))
  }
  for (j in 3:n){    
      N[j] ~ dnorm(mu[j], tau)
  }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dnorm(30000, 1/10000)
  tau ~ dgamma(0.01, 0.01)
  mu_0 ~ dnorm(5000, 1/1000)
 }", file = "M3_ABC.txt")

## Posterior distribution using JAGS
jagsfit_M3 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim','ABC_alpha'), inits, params,
                  n.iter = n_iter, model.file = "M3_ABC.txt", n.chains = 3, n.thin = 1)
fit_M3 = as.mcmc(jagsfit_M3)
class(fit_M3)

# checking posterior distribution convergence
gelman.diag(fit_M3)
gelman.plot(fit_M3)
traceplot(fit_M3)

# storing the posterior distribution of parameters of model M3
C3 = as.matrix(fit_M3)
write.table(C3, "M3_output.txt", row.names = FALSE)

## function for sampling from posterior distribution
draw3 = function(){
  C3[sample(dim(C3)[1], 1, replace = TRUE), -which(colnames(C3) == "deviance")][c(1,4,2,5,3)]
}

############################################################################################################

## for model M4: RL logistic model

## Mean function
mean_M4 = function(N, t,i,h, alpha, r, K ){
  j = i
  s = 0
  for(k in 2:j){
    s = s + (N[k] + N[k-1])*((t[j+1]-t[k-1])^(1-alpha) - (t[j+1]-t[k])^(1-alpha)
                             - (t[j]-t[k-1])^(1-alpha) + (t[j]-t[k])^(1-alpha))
  }
  mu = -N[j] - (h^(alpha-1))*s + 2*(h^alpha)*gamma(2-alpha)*(r*N[j]*(1-(N[j]/K)))
  return(mu)
}

##log-likelihood for M4
log_L_M4 = function(theta){
  s1 = dnorm(N[1], mean = theta[5], sd = sigma_sim, log = TRUE)
  s2 = dnorm(N[2],mean = N[1]+ h*r_sim*N[1]*(1-(N[1]/K_sim)),
             sd = sigma_sim, log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M4(N, t,i, h, theta[1], 
                  theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE) 
  }
  z = s1+ s2 + s3
  return(z)
}

# model txt file of M4
cat("model{
  N[1] ~ dnorm(mu_0,1/sigma_sim^2)
  N[2] ~ dnorm(N[1] + h*r_sim*N[1]*(1-N[1]/K_sim), 1/sigma_sim^2)
  for(i in 3:n){
    for(k in 2:(i-1)){
      temp[i,k] = (N[k] + N[k-1])*((t[i]-t[k-1])^(1-alpha) - (t[i]-t[k])^(1-alpha)
                             - (t[i-1]-t[k-1])^(1-alpha) + (t[i-1]-t[k])^(1-alpha))
    }
    mu[i] = -N[i-1] - (h^(alpha-1))*sum(temp[i,2:(i-1)]) +
                  2*(h^alpha)*exp(loggam(2-alpha))*r*N[i-1]*(1-(N[i-1]/K))
  }
  for (j in 3:n) {
    N[j] ~ dnorm(mu[j], tau)
  }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dnorm(30000, 1/10000)
  tau ~ dgamma(0.01, 0.01)
  mu_0 ~ dnorm(5000, 1/1000)
 }", file = "M4_RL.txt")

## Posterior distribution using JAGS
jagsfit_M4 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim'), inits, params,
                  n.iter = n_iter, model.file = "M4_RL.txt", n.chains = 3, n.thin = 1)
fit_M4 = as.mcmc(jagsfit_M4)
class(fit_M4)

# checking posterior distribution convergence
gelman.diag(fit_M4)
gelman.plot(fit_M4)
traceplot(fit_M4)

# storing the posterior distribution of parameters of model M4
C4 = as.matrix(fit_M4)
write.table(C4, "M4_output.txt", row.names = FALSE)

## function for sampling from posterior distribution
draw4 = function(){
  C4[sample(dim(C4)[1], 1, replace = TRUE), -which(colnames(C4) == "deviance")][c(1,4,2,5,3)]
}

##################################################################################################################################

##for model M5: RL-CF logistic model

## mean function
mean_M5 = function(N,t,i,M_alpha, h, alpha, r,K){
  j = i
  s = 0
  for (k in 2:j){
    s = s + (N[k]+ N[k-1])*( exp((-alpha/(1-alpha))*(t[j]-t[k])) - exp((-alpha/(1-alpha))*(t[j]-t[k-1]))
                             -exp((-alpha/(1-alpha))*(t[j+1]-t[k])) + exp((-alpha/(1-alpha))*(t[j+1]-t[k-1]))  )
  }
  mu = - N[j] + (1/(1-exp((-alpha*h)/(1-alpha))))*s + ((2*alpha*h)/(M_alpha*(1-exp((-alpha*h)/(1-alpha)))))*r*N[j]*(1-(N[j]/K))
  return(mu)
  
}

# log_likelihood for M5
log_L_M5 = function(theta){
  s1 = dnorm(N[1], mean = theta[5], sd = sigma_sim, log = TRUE)
  s2 = dnorm(N[2],mean = N[1]+ h*r_sim*N[1]*(1-(N[1]/K_sim)),
             sd = sigma_sim, log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M5(N, t,i,M_alpha, h, theta[1],
                       theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE)
  }
  z = s1+ s2 + s3
  return(z)
}


# model txt file of M5
cat("model{
  N[1] ~ dnorm(mu_0,1/sigma_sim^2)
  N[2] ~ dnorm(N[1] + h*r_sim*N[1]*(1-N[1]/K_sim), 1/sigma_sim^2)
  for(i in 3:n){
    for(k in 2:(i-1)){
      temp[i,k] = (N[k]+ N[k-1])*( exp((-alpha/(1-alpha))*(t[i-1]-t[k])) - exp((-alpha/(1-alpha))*(t[i-1]-t[k-1]))
                             -exp((-alpha/(1-alpha))*(t[i]-t[k])) + exp((-alpha/(1-alpha))*(t[i]-t[k-1]))  )
    }
    mu[i] = -N[i-1] + (1/(1-exp((-alpha*h)/(1-alpha))))*sum(temp[i,2:(i-1)]) +
                  ((2*alpha*h)/(M_alpha*(1-exp((-alpha*h)/(1-alpha)))))*r*N[i-1]*(1-(N[i-1]/K))
  }
  for (j in 3:n) {
    N[j] ~ dnorm(mu[j], tau)
    }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dnorm(30000, 1/10000)
  tau ~ dgamma(0.01, 0.01)
  mu_0 ~ dnorm(5000, 1/1000)
 }", file = "M5_RL_CF.txt")


## Posterior distribution using JAGS
jagsfit_M5 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim', 'M_alpha'), inits, params,
                  n.iter = n_iter, model.file = "M5_RL_CF.txt", n.chains = 3, n.thin = 1)
fit_M5 = as.mcmc(jagsfit_M5)
class(fit_M5)

# checking posterior distribution convergence
gelman.diag(fit_M5)
gelman.plot(fit_M5)
traceplot(fit_M5)

# storing the posterior distribution of parameters of model M5
C5 = as.matrix(fit_M5)
write.table(C5, "M5_output.txt", row.names = FALSE)

## function for sampling from posterior distribution
draw5 = function(){
  C5[sample(dim(C5)[1], 1, replace = TRUE), -which(colnames(C5) == "deviance")][c(1,4,2,5,3)]
}

##################################################################################################################

## for model M6: Grunwald logistic model

## mean function
mean_M6 = function(N, t, i, h, alpha, r, K){
  j = i
  s = 0
  for (k in 1:j) {
    s = s + ((-1)^(k+1))*(factorial(alpha) / (factorial(alpha-k)* factorial(k)))*N[j-k+1]
  }
  mu = s  + (h^alpha)*r*N[j]*(1-(N[j]/K))
  return(mu)
  
}

# log_likelihood for M6
log_L_M6 = function(theta){
  s1 = dnorm(N[1], mean = theta[5], sd = sigma_sim, log = TRUE)
  s2 = dnorm(N[2],mean = N[1]+ h*r_sim*N[1]*(1-(N[1]/K_sim)),
             sd = sigma_sim, log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M6(N, t,i, h, theta[1],
                                           theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE)
  }
  z = s1+ s2 + s3
  return(z)
}


cat("model{
  N[1] ~ dnorm(mu_0,1/sigma_sim^2)
  N[2] ~ dnorm(N[1] + h*r_sim*N[1]*(1-N[1]/K_sim), 1/sigma_sim^2)
  for(i in 2:(n-1)){
    for(k in 1:i){
      temp[i,k] = ((-1)^(k+1))*(prod(rep(alpha,k)-c(0:(k-1)))/exp(logfact(k)))*N[i-k+1]
      }
    mu[i+1] = sum(temp[i,1:i]) + (h^alpha)*r*N[i]*(1-(N[i]/K))
  }
  for (z in 3:n) {
    N[z] ~ dnorm(mu[z], tau)
    }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dnorm(30000, 1/10000)
  tau ~ dgamma(0.01, 0.01)
  mu_0 ~ dnorm(5000, 1/1000)
 }", file = "M6_Grunwald.txt")

## Posterior distribution using JAGS
jagsfit_M6 = jags(data = c('N', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim'), inits, params,
                  n.iter = n_iter, model.file = "M6_Grunwald.txt", n.chains = 3, n.thin = 1)
fit_M6 = as.mcmc(jagsfit_M6)
class(fit_M6)

# checking posterior distribution convergence
gelman.diag(fit_M6)
gelman.plot(fit_M6)
traceplot(fit_M6)

# storing the posterior distribution of parameters of model M6
C6 = as.matrix(fit_M6)
write.table(C6, "M6_output.txt", row.names = FALSE)

## function for sampling from posterior distribution
draw6 = function(){
  C6[sample(dim(C6)[1], 1, replace = TRUE), -which(colnames(C6) == "deviance")][c(1,4,2,5,3)]
}

######################################################################################
######################################################################################
#### 3. final step: Multimodel inference using RJMCMC


library(rjmcmc)
post_model = rjmcmcpost(post.draw = list(draw1, draw2, draw3, draw4, draw5, draw6), 
                        g = list(g, g, g, g, g,g), ginv = list(ginv, ginv, ginv, ginv, ginv, ginv),
                        likelihood = list(log_L_M1,log_L_M2, log_L_M3,log_L_M4, log_L_M5, log_L_M6),
                        param.prior = list(log_prior, log_prior, log_prior, log_prior, log_prior,log_prior),
                        model.prior = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6),chainlength = 2e4, TM.thin = 10 , 
                        save.all = TRUE)
write.table(post_model$result$`Posterior Model Probabilities`, "post_mod_prob_M1_1234.txt", )
write.table(post_model$progress$prb, "post_mod_progress_M1_1234.txt", row.names = FALSE)
capture.output(post_model$result, file = "post_rjmcmc_result_M1_1234.txt")

####################################################################################
####################################################################################
###### 4. visualization of outputs by plotting them

#setwd("specify path of stored file here")

####### barplot of posterior model probabilities
post_prob = read.table("post_mod_prob_M1_1234.txt", header = TRUE)
post_prob_temp = c(post_prob[1:6,1])

post_progress = as.matrix(read.csv("post_mod_progress_M1_1234.txt", header = TRUE, sep = " "))
matplot(post_progress, type = "l", lwd = 2, col = 6:1, ylab = "posterior model probability")


# barplot of posterior model probabilities

library(ggplot2)
Mod_labels = c("M1", "M2", "M3", "M4", "M5", "M6")
df = data.frame(Model = Mod_labels, len = Re(post_prob_temp))

p = ggplot(data = df, aes(x = Model, y = len)) +
  labs(x = "Model", y = "Posterior Model Probability") +
  geom_bar(stat = "identity", fill = "aquamarine4") +
  geom_text(aes(label = round(len, 4)), vjust = -0.3, size = 3.5) +
  theme_minimal() +
  scale_x_discrete(labels = Mod_labels) +
  scale_y_continuous(expand = c(0, 0)) +  # Ensure the y = 0 line aligns with x-axis
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    panel.grid.major = element_blank(),  # Remove all grid lines
    panel.grid.minor.y = element_line(colour = "white"),
    axis.line.x = element_line(color = "black")  # Add the x-axis line in black
  ) + coord_cartesian(ylim = c(0, 1.1))
# to save the plot
filename = paste0("POSTERIOR_MODEL_PROBABILITIES_BARPLOT_simulation study.jpeg")
ggsave(filename, plot = p, width = 6, height = 4, dpi = 400) # to save the plot



## posterior histogram of model parameters along with posterior density curve for all the models

for (m in 1:6) {
  #setwd("specify path where files has been stored")
  # output reading 
  filename = paste0("M",m, "_output", ".txt")
  C_post = read.table(filename, header = TRUE)
  df = data.frame(C_post[,c(1,5, 3, 6, 4)])
  
  library(ggplot2)
  
  # Define a light background color
  light_background_color = "white"
  
  # posterior histogram along with posterior density plot for parameter alpha
  h1 = ggplot(df, aes(x = alpha, after_stat(density))) + 
    geom_histogram(aes(y = after_stat(density)),
                   colour = 1, fill = "#F3DF6C", alpha = 0.6) +
    geom_density(lwd = 1.5, colour = "#E3B710", alpha = 0.25)+
    labs(x = expression(alpha), y = "density") +
    theme_classic() +  
    ggtitle(bquote(M[.(m)]))+
    theme(panel.background = element_rect(fill = light_background_color),  # Set background color
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  # to save the plot
  
  filename = paste0("POST_hist_alpha_M", m,"_", "simulation study", ".jpeg")
  ggsave(filename, plot = h1, width = 6, height = 4, dpi = 400) # to save the plot
  
  
  # posterior histogram along with posterior density plot for parameter r
  h2 =ggplot(df, aes(x = r, ..density..)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "#ABDDDE", alpha = 0.2) +
    geom_density(lwd = 1.5, colour = "#46ACC8", alpha = 0.25)+
    labs(x = expression(r), y = "density") +
    theme_classic() + 
    ggtitle(bquote(M[.(m)]))+
    theme(panel.background = element_rect(fill = light_background_color),  # Set background color
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  # to save the plot
  filename = paste0("POST_hist_r_M", m,"_", "simulation study", ".jpeg")
  ggsave(filename, plot = h2, width = 6, height = 4, dpi = 400) # to save the plot
  
  # posterior histogram along with posterior density plot for parameter K
  h3 = ggplot(df, aes(x = K, ..density..)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "#7294D4", alpha = 0.9) +
    geom_density(lwd = 1.5, colour = "blue2", alpha = 0.25)+
    labs(x = expression(K), y = "density") +
    theme_classic() +  
    ggtitle(bquote(M[.(m)]))+
    theme(panel.background = element_rect(fill = light_background_color),  # Set background color
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  # to save the plot
  filename = paste0("POST_hist_K_M", m,"_", "simulation study", ".jpeg")
  ggsave(filename, plot = h3, width = 6, height = 4, dpi = 400) # to save the plot
  
  # posterior histogram along with posterior density plot for parameter tau
  h4 = ggplot(df, aes(x = tau, ..density..)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "brown1", alpha = 0.2) +
    geom_density(lwd = 1.5, colour = "brown3", alpha = 0.25)+
    labs(x = expression(tau == 1/sigma^2), y = "density") +
    theme_classic() +  
    ggtitle(bquote(M[.(m)]))+
    theme(panel.background = element_rect(fill = light_background_color),  # Set background color
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  # to save the plot
  filename = paste0("POST_hist_tau_M", m,"_", "simulation study", ".jpeg")
  ggsave(filename, plot = h4, width = 6, height = 4, dpi = 400) # to save the plot
  
  # posterior histogram along with posterior density plot for parameter tau
  h5 = ggplot(df, aes(x = mu_0, ..density..)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "lightgreen", alpha = 0.2) +
    geom_density(lwd = 1.5, colour = "green4", alpha = 0.25)+
    labs(x = expression(mu[0]), y = "density") +
    theme_classic() +  
    ggtitle(bquote(M[.(m)]))+
    theme(panel.background = element_rect(fill = light_background_color),  # Set background color
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  # to save the plot
  filename = paste0("POST_hist_mu_0_M", m,"_", "simulation study", ".jpeg")
  ggsave(filename, plot = h5, width = 6, height = 4, dpi = 400) # to save the plot
  
}

