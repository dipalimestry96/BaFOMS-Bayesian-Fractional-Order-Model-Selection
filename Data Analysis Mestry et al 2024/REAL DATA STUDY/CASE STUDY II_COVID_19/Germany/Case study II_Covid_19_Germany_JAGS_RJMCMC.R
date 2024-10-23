## R code for case study II: Covid-19 data for country Germany
## we have considered the total cases under natural logarithmic scale for modelling using fractional logistic model
## Data file name: "Covid_Germany_data.csv"
## keep this file in the working directionary
## data is from 31-12-2019 to 08-07-2020
## till the data 08-03-2020, there is no significant relative change in the total cases
## so we have considered the data of total cases from 09-03-2020 to 08-07-2020 of length 122 for the study
## We have total six fractional logistic models under six different fractional derivative (FD) definitions
##   M1: Under Caputo FD definition
##   M2: Under Caputo-Fabrizio FD definition
##   M3: Under ABC FD definition
##   M4: Under Riemann-Liouville FD definition
##   M5: Under CF-RL FD definition
##   M6: Under Grunwal-Letnikov definition

## The sequence of code is as follows:
## 1. Getting posterior distribution of model parameters and posterior model probabilities
## 2. Barplot of posterior model probabilities
## 3. Total cases data (log_e scale) along with fitted best model growth curve
## 4. posterior histogram of model parameters along with posterior density curve for all the models
## 5. fitting Covid-19 new cases using posterior mode
## 6. Parameter estimate (mode) and 95% posterior credible interval

###########################################################################################
###########################################################################################

#### 1. Getting posterior distribution of model parameters and posterior model probabilities

set.seed(12345)   # for reproducibility of all results

## Importing all needed packages
library(rjags)
library(R2jags)
library(coda)
library(madness)
library(pracma)

#setwd("specify the path of folder containing data file 'Covid_Germany_data.csv' ")
setwd("specify path here")

##data reading part
covid_data = read.csv("Covid_Germany_data.csv")  # covid-19 data of the country Germany
covid_data = as.data.frame(covid_data)

par(mfrow = c(1,1))
data = covid_data$total_cases
plot(data, col = "blue", pch = 19)

## here we have observed that there is no significant change in the total cases at the initial phase so
## we have discarded starting some values for the analysis.
## analysis are carried using data on natural logarithmic scale as the eventually total cases becomes too large
## so making computationally feasible we have taken the log scale

N = log(data[70:nrow(covid_data)]) 
which(is.na(N))   # to check the NA values
length(N)
plot( N, main = "Total covid-19 cases for Germany", type="p", col = "blue", lwd = 2, xlab = "", pch = 19)

n = length(N)   # length of data
M_alpha = 1     # normalizing function   
ABC_alpha = 1   # normalizing function


var_samp = (1/(n-1))*sum((N-mean(N))^2)   #sample variance
sigma_samp = sqrt(var_samp)    # SD of sample variance


t = numeric(length = length(N))   # time vector
t[1] = 0
h = 1    # time step
for (i in 1:(length(N)-1)) {
  t[i+1] = t[i] + h
  
}

N_1 = N[1]
max_N = max(N)

## theta = c(alpha, r, K, tau) where tau = 1/sigma^2
## prior specification of all model parameters
# alpha ~ unif(0,1)
# r ~ gamma(2,1/2)
# K ~ unif(2, max(N) + 10*sigma_samp)
# tau ~ gamma(0.01, 0.01)


# log-prior for all model parameters
log_prior = function(theta){
  sum(dunif(theta[1], min = 0, max = 1, log = TRUE),
      dgamma(theta[2], 2, 1/2, log = TRUE) ,
      dunif(theta[3], min = 2, max = max_N + 10*sigma_samp, log = TRUE),
      dgamma(theta[4], 0.01, 0.01, log = TRUE ))
}

##initial values
inits = function(){
  list(alpha = 0.5, r = 3, K = max_N, tau = 1/var_samp)
}

params = c("alpha", "r", "K", "tau")  # parameters to be estimated

## bijection and its inverse
g = function(psi){ psi }
ginv = function(theta){ theta }


n_iter = 1e5  # no of iteration for jags
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
  s1 = dnorm(N[1], mean = N_1, sd = 1/sqrt(theta[4]), log = TRUE)
  s2 = dnorm(N[2],mean = N_1+ h*theta[2]*N_1*(1-(N_1/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[1] ~ dnorm(N_1, tau)
  N[2] ~ dnorm(N_1 + h*r*N_1*(1-N_1/K), tau)
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
  K ~ dunif(2, max_N + 10*sigma_samp)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M1_Caputo.txt")

## Posterior distribution using JAGS
jagsfit_M1 = jags(data = c('N', 't', 'n','sigma_samp', 'h','max_N', 'N_1'), inits, params,
                  n.iter = n_iter, model.file = "M1_Caputo.txt", n.chains = 3, n.thin = 10)
fit_M1 = as.mcmc(jagsfit_M1) # mcmc after burnin period=iter/2

# checking posterior distribution convergence
gelman.diag(fit_M1)
gelman.plot(fit_M1)
par(mfrow = c(2,3))
traceplot(fit_M1)

C1 = as.matrix(fit_M1)

par(mfrow = c(2,2))
# to check autocorrelation
acf(C1[,1])   
acf(C1[,3])
acf(C1[,4])
acf(C1[,5])
write.table(C1, "M1_output.txt", row.names = FALSE) # storing posterior values

## function for sampling from posterior
draw1 = function(){
  C1[sample(dim(C1)[1], 1, replace = TRUE), -which(colnames(C1) == "deviance")][c(1,3,2,4)]
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
  s1 = dnorm(N[1], mean = N_1, sd = 1/sqrt(theta[4]), log = TRUE)
  s2 = dnorm(N[2],mean = N_1+ h*theta[2]*N_1*(1-(N_1/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M2(N[i], N[i-1],M_alpha,h,theta[1], theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE) 
  }
  z = s1+ s2 + s3
  return(z)
}

# model txt file of M2
cat("model{
  N[1] ~ dnorm(N_1, tau)
  N[2] ~ dnorm(N_1 + h*r*N_1*(1-N_1/K), tau)
  for(i in 3:n){
    mu[i] = N[i-1] + ((1-alpha)/M_alpha + (3*alpha*h)/(2*M_alpha))*r*N[i-1]*(1-(N[i-1]/K)) 
               - ((1-alpha)/M_alpha + (alpha*h)/(2*M_alpha))*r*N[i-2]*(1-(N[i-2]/K))
    }
   for (j in 3:n){    
      N[j] ~ dnorm(mu[j], tau)
  }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dunif(2, max_N + 10*sigma_samp)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M2_Caputo_Fabrizio.txt")


## Posterior distribution using JAGS
jagsfit_M2 = jags(data = c('N', 'n', 'h','M_alpha','N_1', 'max_N', 'sigma_samp'), inits, params,
                  n.iter = n_iter, model.file = "M2_Caputo_Fabrizio.txt", n.chains = 3, n.thin = 20)
fit_M2 = as.mcmc(jagsfit_M2)

# checking convergence of posterior distribution
gelman.diag(fit_M2)
gelman.plot(fit_M2)
par(mfrow = c(2,3))
traceplot(fit_M2)

C2 = as.matrix(fit_M2)
# to check autocorrelation
par(mfrow = c(2,2))
acf(C2[,1])
acf(C2[,3])
acf(C2[,4])
acf(C2[,5])

write.table(C2, "M2_output.txt", row.names = FALSE) # storing posterior values

## function for sampling from posterior
draw2 = function(){
  C2[sample(dim(C2)[1], 1, replace = TRUE), -which(colnames(C2) == "deviance")][c(1,3,2,4)]
}

###################################################################################################################
## for model M3: ABC logistic model

##mean function
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
  s1 = dnorm(N[1], mean = N_1, sd = 1/sqrt(theta[4]), log = TRUE)
  s2 = dnorm(N[2],mean = N_1+ h*theta[2]*N_1*(1-(N_1/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M3(N[i], N[i-1],t[i+1],t[i],t[i-1],h,ABC_alpha,theta[1], theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE) 
  }
  z = s1+ s2 + s3
  return(z)
}


# model txt file of M3
cat("model{
  N[1] ~ dnorm(N_1, tau)
  N[2] ~ dnorm(N_1 + h*r*N_1*(1-N_1/K), tau)
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
  K ~ dunif(2, max_N + 10*sigma_samp)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M3_ABC.txt")

## Posterior distribution using JAGS
jagsfit_M3 = jags(data = c('N', 't', 'n', 'h','ABC_alpha','N_1', 'max_N', 'sigma_samp'), inits, params,
                  n.iter = n_iter, model.file = "M3_ABC.txt", n.chains = 3, n.thin = 20)
fit_M3 = as.mcmc(jagsfit_M3)
class(fit_M3)

# checking convergence of posterior distribution
gelman.diag(fit_M3)
par(mfrow = c(2,3))
gelman.plot(fit_M3)
traceplot(fit_M3)

C3 = as.matrix(fit_M3)

# to check autocorrelation
par(mfrow = c(2,2))
acf(C3[,1])
acf(C3[,3])
acf(C3[,4])
acf(C3[,5])

write.table(C3, "M3_output.txt", row.names = FALSE) # storing posterior values

##sampling function from posterior
draw3 = function(){
  C3[sample(dim(C3)[1], 1, replace = TRUE), -which(colnames(C3) == "deviance")][c(1,3,2,4)]
}

############################################################################################################
## for model M4: RL logistic model

##Mean function
mean_M4 = function(N, t,i, h, alpha, r, K ){
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
  s1 = dnorm(N[1], mean = N_1, sd = 1/sqrt(theta[4]), log = TRUE)
  s2 = dnorm(N[2],mean = N_1+ h*theta[2]*N_1*(1-(N_1/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[1] ~ dnorm(N_1, tau)
  N[2] ~ dnorm(N_1 + h*r*N_1*(1-N_1/K), tau)
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
  K ~ dunif(2, max_N + 10*sigma_samp)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M4_RL.txt")

##initial values
inits_4 = function(){
  list(alpha = 0.5, r = 3, K = max_N, tau = 1/var_samp)
}

## Posterior distribution using JAGS
jagsfit_M4 = jags(data = c('N', 't', 'n', 'h','N_1',"max_N","sigma_samp"), inits_4, params,
                  n.iter = n_iter, model.file = "M4_RL.txt", n.chains = 3, n.thin = 40)
fit_M4 = as.mcmc(jagsfit_M4)
class(fit_M4)

#checking convergence of posterior distribution
gelman.diag(fit_M4)
par(mfrow = c(2,3))
gelman.plot(fit_M4)
traceplot(fit_M4)

C4 = as.matrix(fit_M4)

# to check autocorrelation
par(mfrow = c(2,2))
acf(C4[,1])
acf(C4[,3])
acf(C4[,4])
acf(C4[,5])

write.table(C4, "M4_output.txt", row.names = FALSE) # storing posterior values

## function for sampling from posterior distribution
draw4 = function(){
  C4[sample(dim(C4)[1], 1, replace = TRUE), -which(colnames(C4) == "deviance")][c(1,3,2,4)]
}

##################################################################################################################################

##for model M5: RL-CF logistic model

##mean function
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
  s1 = dnorm(N[1], mean = N_1, sd = 1/sqrt(theta[4]), log = TRUE)
  s2 = dnorm(N[2],mean = N_1+ h*theta[2]*N_1*(1-(N_1/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[1] ~ dnorm(N_1, tau)
  N[2] ~ dnorm(N_1 + h*r*N_1*(1-N_1/K), tau)
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
  K ~ dunif(2, max_N + 10*sigma_samp)
  tau ~ dgamma(0.01, 0.01)
  }", file = "M5_RL_CF.txt")

##initial values
inits_5 = function(){
  list(alpha = 0.5, r = 0.2, K = max_N, tau = 1/var_samp)
}

## Posterior distribution using JAGS
jagsfit_M5 = jags(data = c('N', 't', 'n', 'h', 'M_alpha','N_1',"max_N","sigma_samp"), inits_5, params,
                  n.iter = n_iter, model.file = "M5_RL_CF.txt", n.chains = 3,
                  n.thin = 50)
fit_M5 = as.mcmc(jagsfit_M5)
class(fit_M5)

# checking convergence of posterior distribution
gelman.diag(fit_M5)
par(mfrow = c(2,3))
gelman.plot(fit_M5)
gelman.diag(fit_M5)
traceplot(fit_M5)

C5 = as.matrix(fit_M5)

# to check autocorrelation
par(mfrow = c(2,2))
acf(C5[,1])
acf(C5[,3])
acf(C5[,4])
acf(C5[,5])

write.table(C5, "M5_output.txt", row.names = FALSE) # storing posterior values

## function for sampling from the posterior distribution
draw5 = function(){
  C5[sample(dim(C5)[1], 1, replace = TRUE), -which(colnames(C5) == "deviance")][c(1,3,2,4)]
}


#######################################################################################################
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
  s1 = dnorm(N[1], mean = N_1, sd = 1/sqrt(theta[4]), log = TRUE)
  s2 = dnorm(N[2],mean = N_1+ h*theta[2]*N_1*(1-(N_1/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M6(N, t,i, h, theta[1],
                                           theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE)
  }
  z = s1+ s2 + s3
  return(z)
}

# Model file for M6
cat("model{
  N[1] ~ dnorm(N_1, tau)
  N[2] ~ dnorm(N_1 + h*r*N_1*(1-N_1/K), tau)
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
  K ~ dunif(2, max_N + 10*sigma_samp)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M6_Grunwald.txt")

##initial values
inits_6 = function(){
  list(alpha = 0.7, r = 2, K = max_N, tau = 1/var_samp)
}

## Posterior distribution using JAGS
jagsfit_M6 = jags(data = c('N', 'n', 'h', 'N_1',"max_N", "sigma_samp"), inits_6, params,
                  n.iter = n_iter, model.file = "M6_Grunwald.txt", n.chains = 3, n.thin = 40)
fit_M6 = as.mcmc(jagsfit_M6)
class(fit_M6)

# checking convergence of posterior distribution
gelman.diag(fit_M6)
par(mfrow = c(2,3))
gelman.plot(fit_M6)
traceplot(fit_M6)

C6 = as.matrix(fit_M6)
# to check autocorrelation
par(mfrow = c(2,2))
acf(C6[,1])
acf(C6[,3])
acf(C6[,4])
acf(C6[,5])

write.table(C6, "M6_output.txt", row.names = FALSE) # for storing the posterior values

## function for sampling from posterior values
draw6 = function(){
  C6[sample(dim(C6)[1], 1, replace = TRUE), -which(colnames(C6) == "deviance")][c(1,3,2,4)]
}
###############################################################################################################

####final step: RJMCMC

library(rjmcmc)
post_model = rjmcmcpost(post.draw = list(draw1, draw2, draw3, draw4, draw5, draw6), 
                        g = list(g, g, g, g, g,g), ginv = list(ginv, ginv, ginv, ginv, ginv, ginv),
                        likelihood = list(log_L_M1,log_L_M2, log_L_M3,log_L_M4, log_L_M5, log_L_M6),
                        param.prior = list(log_prior, log_prior, log_prior, log_prior, log_prior,log_prior),
                        model.prior = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6),chainlength = 5e4, TM.thin = 10 , 
                        save.all = TRUE)

write.table(post_model$result$`Posterior Model Probabilities`, "post_mod_prob_real_data_Germany.txt" )
write.table(post_model$progress$prb, "post_mod_progress_real_data_Germany.txt", row.names = FALSE)
capture.output(post_model$result, file = "post_rjmcmc_result_real_data_Germany.txt")

########################################################################################################
########################################################################################################

#### 2. Barplot of posterior model probabilities

set.seed(12345)   # for reproducibility of all results
par(mfrow = c(1,1))
post_prob = read.table("post_mod_prob_real_data_Germany.txt", header = TRUE)
post_prob_temp = c(Re(post_prob[1:6,1]))

Mod_labels = c("M1", "M2", "M3", "M4", "M5", "M6") # model labels
df = data.frame(Model = Mod_labels,len=post_prob_temp) # data frame for plotting
Mod_labels = c(expression(M[1], M[2], M[3], M[4], M[5], M[6])) # labels in mathematical expression

# posterior model probabilities barplot
library(ggplot2)
p = ggplot(data = df, aes(x = Model, y = len)) +
  labs(x = "Model", y = "Posterior Model Probability") +
  geom_bar(stat = "identity", fill = "aquamarine4") +
  geom_text(aes(label = round(len, 4)), vjust = -0.3, size = 4.5) +
  theme_minimal() +
  scale_x_discrete(labels = Mod_labels) +
  scale_y_continuous(expand = c(0, 0)) +  # Ensure the y = 0 line aligns with x-axis
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    panel.grid.major = element_blank(),  # Remove all grid lines
    panel.grid.minor.y = element_line(colour = "white"),
    axis.line.x = element_line(color = "black")  # Add the x-axis line in black
  ) + coord_cartesian(ylim = c(0, 1.1))
# to save the plot
filename = paste0("posterior model prob_barplot_Germany covid_19_data(total cases).jpeg")
ggsave(filename, plot = p, width = 6, height = 4, dpi = 600) # to save the plot


########################################################################################################
########################################################################################################

#### 3. Total cases data (log_e scale) along with fitted best model growth curve

C6 = read.table("M6_output.txt", header = TRUE) # posterior values of model M6 parameters
C6 = as.matrix(C6)

# posterior mode as a parameter estimate
density_est_alpha = density(C6[,1])
alpha = density_est_alpha$x[which.max(density_est_alpha$y)]

density_est_r = density(C6[,4])
r = density_est_r$x[which.max(density_est_r$y)]

density_est_K = density(C6[,3])
K = density_est_K$x[which.max(density_est_K$y)]

density_est_tau = density(C6[,5])
tau = density_est_tau$x[which.max(density_est_tau$y)]
theta = c(alpha, r, K, tau)

data = numeric(length(N)) # array for fittted values
data[1] = N[1]
data[2] = data[1]+ h*theta[2]*data[1]*(1-(data[1]/theta[3]))

for (i in 2:(n-1)) {
  data[i+1] = mean_M6(data, t,i, h, theta[1], theta[2], theta[3])
}

df = data.frame(t = t, y1 = N, y2 = data) # dataframe for plotting

# total cases data (log_e scale) along with the best fitted G-L fractional growth model curve
ppp = ggplot(df, aes(x = t)) +
  theme_classic() +
  geom_point(aes(y = y1, color = "Germany: Covid-19 Total cases"), pch = 19, size = 3) +
  geom_line(aes(y = y2, color = "G-L Fractional Growth Model"), linewidth = 1.2) +
  labs(x = "t (days): 9 March 2020 - 8 July 2020", y = expression(log[e]("Covid-19 Total Cases"))) +
  scale_color_manual(
    values = c("Germany: Covid-19 Total cases" = "blue", "G-L Fractional Growth Model" = "red")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = c(NA, 19), linetype = c("solid", "blank"),size = c(1.2, 3)))  # Use NA to exclude shapes for lines
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = c(0.65, 0.2),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_blank()
  )
# to save the plot
filename = paste0("Total cases with fitted best model curve for_Germany covid_19_data(total cases).jpeg")
ggsave(filename, plot = ppp, width = 6, height = 4, dpi = 600) # to save the plot

###########################################################################################################
###########################################################################################################

## 4. posterior histogram of model parameters along with posterior density curve for all the models

for (m in 1:6) {
  # output reading
  filename = paste0("M",m, "_output", ".txt")
  C_post = read.table(filename, header = TRUE)
  df = data.frame(C_post[,c(1,4,3,5)])
  
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
  filename = paste0("POST_hist_alpha_M", m,"_", "Germany_covid_19_data(total cases)", ".jpeg")
  ggsave(filename, plot = h1, width = 6, height = 4, dpi = 600) # to save the plot
  
  
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
  filename = paste0("POST_hist_r_M", m,"_", "Germany_covid_19_data(total cases)", ".jpeg")
  ggsave(filename, plot = h2, width = 6, height = 4, dpi = 600) # to save the plot
  
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
  filename = paste0("POST_hist_K_M", m,"_", "Germany_covid_19_data(total cases)", ".jpeg")
  ggsave(filename, plot = h3, width = 6, height = 4, dpi = 600) # to save the plot
  
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
  filename = paste0("POST_hist_tau_M", m,"_", "Germany_covid_19_data(total cases)", ".jpeg")
  ggsave(filename, plot = h4, width = 6, height = 4, dpi = 600) # to save the plot
  
}

#####################################################################################################
#####################################################################################################

#### 5. fitting Covid-19 new cases using posterior mode

## Posterior mode as a parameter estimate of parameters 
filename = paste0("M",6, "_output", ".txt")
C6 = read.table(filename, header = TRUE)[,c(1,4,3,5)]

## Posterior mode as a parameter estimate of parameters 
alpha = mode(C6[,1])
r = mode(C6[,2])
K = mode(C6[,3])
tau = mode(C6[,4])

theta = c(alpha, r, K, tau)

data = numeric(length(N))
data[1] = N[1]
data[2] = data[1]+ h*theta[2]*data[1]*(1-(data[1]/theta[3]))

for (i in 2:(n-1)) {
  data[i+1] = mean_M6(data, t,i, h, theta[1], theta[2], theta[3])
}

library(ggplot2)

# Calculate differences
diff = (exp(data[2:length(data)]) - exp(data[1:(length(data)-1)]))/h
# Create a new data frame for the differences
diff_data = data.frame(t = t[-(length(t)-1)], data_point = covid_data$new_cases[70:(nrow(covid_data)-1)], diff = diff)


pp = ggplot() +
  geom_point(data = diff_data, aes(x = t, y = data_point, color = "Points"), pch = 19, size = 3) +
  geom_line(data = diff_data, aes(x = t, y = diff, color = "Line"), lwd = 1.2) +
  labs(y = "Germany: Covid-19 New Cases", x = "t(days): 9 March 2020 - 8 July 2020") +
  scale_color_manual(values = c("Points" = "blue", "Line" = "red"),name = "",
                     labels = c("N'(t) using G-L FDE","Germany: new COVID-19 cases" )) +
  guides(
    color = guide_legend(override.aes = list(shape = c(NA, 19), linetype = c("solid", "blank"),size = c(1.2, 3)))  # Use NA to exclude shapes for lines
  )+
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = c(0.7, 0.9),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 13))

# to save the plot
filename = paste0("New cases fitting_Germany covid_19_data(total cases).jpeg")
ggsave(filename, plot = pp, width = 6, height = 4, dpi = 600) # to save the plot

############################################################################################################
############################################################################################################
#### 6. Parameter estimate (mode) and 95% posterior credible interval


GPDD_ID = 9330
# function for calculating mode
mode = function(w){
  density_est = density(w)
  return(density_est$x[which.max(density_est$y)])
}

post_mod_prob = as.numeric(read.table("post_mod_prob_real_data_Germany.txt", header = TRUE)[1:6,1]) #posterior model probabilities
m = which.max(post_mod_prob)
filename = paste0("M",6, "_output", ".txt")
C = read.table(filename, header = TRUE)[,c(1,4,3,5)]   # posterior values of parameters obtained using the JAGS
head(C)

param_est_mode = numeric(4)   # array to store the parameter estimate
param_post_CI = matrix(NA, ncol = 2, nrow = 4)  # array to store 95% posterior credible interval (row-wise, for alpha, r, K and tau) 
for (i in 1:4) {
  param_est_mode[i] = mode(C[,i])
  param_post_CI[i,1] = as.numeric(quantile(C[,i],probs = 0.025))
  param_post_CI[i,2] = as.numeric(quantile(C[,i], probs = 0.975))
}
round(param_est_mode,4)
param_post_CI


