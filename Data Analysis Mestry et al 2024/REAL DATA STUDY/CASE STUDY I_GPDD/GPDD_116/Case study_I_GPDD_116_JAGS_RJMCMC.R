## R code for case study I: Global population dynamics database (GPDD)
## One data set are extracted from GPDD and analysis are carried out for this dataset
## GPDD ID = 116
## GPDD ID = 116 : # a mammalian species, namely Ursus Americanus 
# Commonly known as the American black bear, GPDD Main ID: 116 
# spanning a duration of 12 years (1970-1981)

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
## 3. posterior histogram of model parameters along with posterior density curve for all the models
## 4. Parameter estimate (mode) and 95% posterior credible interval



##########################################################################################################
#### 1. Getting posterior distribution of model parameters and posterior model probabilities

set.seed(1234) # for reproducibility of results

# installation of rgpdd package
#install.packages("remotes",dependencies = TRUE)
#remotes::install_github("ropensci/rgpdd")

# loading required libraries
library(rgpdd) # loading package rgpdd for extracting data from GPDD
library(R2jags)
library(coda)
library(madness)


 
setwd("specify path of local folder on system for saving the output")

GPDD_ID = 116

## data extraction to excel file
data = subset(gpdd_data, MainID == GPDD_ID)
GPDD_sample_year_and_population_data = data.frame(data$SampleYear, data$Population)
filename = paste0("GPDD Data_MainID_", GPDD_ID,'.csv')
write.csv(GPDD_sample_year_and_population_data,filename)

data = subset(gpdd_data, MainID == GPDD_ID)$Population

N = data      # population density data
n = length(N)  # length of data
M_alpha = 1  # normalising function value
ABC_alpha = 1 # normalising function value

var_samp = (1/(n-1))*sum((N-mean(N))^2)   #sample variance
sigma_samp = sqrt(var_samp)    # SD of sample variance

t = numeric(length = length(N))   # time grid
t[1] = 0
h = 0.01                # time step
for (j in 1:(length(N)-1)) {
  t[j+1] = t[j] + h
}

N_1 = N[1]      # mean for initial population size
max_N = max(N)  # maximum of population densities
min_N = min(N)  # minimum of population densities

Max = 10000

params = c("alpha", "r", "K", "tau")   # parameters to be estimated

log_prior = function(theta){   #log prior density function
  sum(dunif(theta[1], min = 0, max = 1, log = TRUE),
      dgamma(theta[2], 2, 1/2, log = TRUE) ,
      dunif(theta[3], min = min_N, max = Max, log = TRUE),
      dgamma(theta[4], 0.01, 0.01, log = TRUE ))
}

## initial value of parameters
inits = function(){
  list(alpha = 0.7, r = 0.7, K = max_N, tau = 1/var_samp)
}



## bijection and its inverse for RJMCMC to transition between model space
g = function(psi){ psi }
ginv = function(theta){ theta }

n_iter = 2e5  # no of iteration for jags


##########################################################################################
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
  s2 = dnorm(N[2],mean = N[1]+ h*theta[2]*N[1]*(1-(N[1]/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[2] ~ dnorm(N[1] + h*r*N[1]*(1-N[1]/K), tau)
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
  K ~ dunif(min_N, Max)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M1_Caputo.txt")

## Posterior distribution using JAGS
jagsfit_M1 = jags(data = c('N', 't', 'n', 'h','min_N', 'N_1','Max'), inits, params,
                  n.iter = n_iter, model.file = "M1_Caputo.txt", n.chains = 2, n.thin = 10)
fit_M1 = as.mcmc(jagsfit_M1) # mcmc after burnin period=iter/2

C1 = as.matrix(fit_M1)
filename = paste0("post_M1_GPDD_ID_", GPDD_ID, ".txt")
write.table(C1, filename, row.names = FALSE)

##function for sampling from posterior distribution
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
  s2 = dnorm(N[2],mean = N[1]+ h*theta[2]*N[1]*(1-(N[1]/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[2] ~ dnorm(N[1] + h*r*N[1]*(1-N[1]/K), tau)
  for(i in 3:n){
    mu[i] = N[i-1] + ((1-alpha)/M_alpha + (3*alpha*h)/(2*M_alpha))*r*N[i-1]*(1-(N[i-1]/K)) 
               - ((1-alpha)/M_alpha + (alpha*h)/(2*M_alpha))*r*N[i-2]*(1-(N[i-2]/K))
    }
   for (j in 3:n){    
      N[j] ~ dnorm(mu[j], tau)
  }
  alpha ~ dunif(0,1)
  r ~ dgamma(2,1/2)
  K ~ dunif(min_N, Max)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M2_Caputo_Fabrizio.txt")

## Posterior distribution using JAGS
jagsfit_M2 = jags(data = c('N', 'n', 'h','M_alpha','N_1', 'min_N','Max'), inits, params,
                  n.iter = n_iter, model.file = "M2_Caputo_Fabrizio.txt", n.chains = 2, n.thin = 10)
fit_M2 = as.mcmc(jagsfit_M2)
class(fit_M2)

C2 = as.matrix(fit_M2)
filename = paste0("post_M2_GPDD_ID_", GPDD_ID, ".txt")
write.table(C2, filename, row.names = FALSE)

## function for sampling from posterior distribution
draw2 = function(){
  C2[sample(dim(C2)[1], 1, replace = TRUE), -which(colnames(C2) == "deviance")][c(1,3,2,4)]
}

#############################################################################################################################

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
  s2 = dnorm(N[2],mean = N[1]+ h*theta[2]*N[1]*(1-(N[1]/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[2] ~ dnorm(N[1] + h*r*N[1]*(1-N[1]/K), tau)
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
  K ~ dunif(min_N, Max)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M3_ABC.txt")

## Posterior distribution using JAGS
jagsfit_M3 = jags(data = c('N', 't', 'n', 'h','ABC_alpha','N_1', 'min_N', 'Max'), inits, params,
                  n.iter = n_iter, model.file = "M3_ABC.txt", n.chains = 2, n.thin = 10)
fit_M3 = as.mcmc(jagsfit_M3)
class(fit_M3)

C3 = as.matrix(fit_M3)
filename = paste0("post_M3_GPDD_ID_", GPDD_ID, ".txt")
write.table(C3, filename, row.names = FALSE)

##function for sampling from posterior distribution
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
  s2 = dnorm(N[2],mean = N[1]+ h*theta[2]*N[1]*(1-(N[1]/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[2] ~ dnorm(N[1] + h*r*N[1]*(1-N[1]/K), tau)
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
  K ~ dunif(min_N, Max)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M4_RL.txt")

## Posterior distribution using JAGS
jagsfit_M4 = jags(data = c('N', 't', 'n', 'h','N_1','min_N','Max'), inits, params,
                  n.iter = n_iter, model.file = "M4_RL.txt", n.chains = 2, n.thin = 10)
fit_M4 = as.mcmc(jagsfit_M4)
class(fit_M4)

C4 = as.matrix(fit_M4)
filename = paste0("post_M4_GPDD_ID_", GPDD_ID, ".txt")
write.table(C4, filename, row.names = FALSE)

##function for sampling from posterior distribution
draw4 = function(){
  C4[sample(dim(C4)[1], 1, replace = TRUE), -which(colnames(C4) == "deviance")][c(1,3,2,4)]
}

##################################################################################################################################

##for model M5: CF-RL logistic model

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
  s2 = dnorm(N[2],mean = N[1]+ h*theta[2]*N[1]*(1-(N[1]/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
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
  N[2] ~ dnorm(N[1] + h*r*N[1]*(1-N[1]/K), tau)
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
  K ~ dunif(min_N, Max)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M5_RL_CF.txt")

## Posterior distribution using JAGS
jagsfit_M5 = jags(data = c('N', 't', 'n', 'h', 'M_alpha','N_1','min_N','Max'), inits, params,
                  n.iter = n_iter, model.file = "M5_RL_CF.txt", n.chains = 2,
                  n.thin = 10)
fit_M5 = as.mcmc(jagsfit_M5)
class(fit_M5)

C5 = as.matrix(fit_M5)
filename = paste0("post_M5_GPDD_ID_", GPDD_ID, ".txt")
write.table(C5, filename, row.names = FALSE)

##function for sampling from posterior distribution
draw5 = function(){
  C5[sample(dim(C5)[1], 1, replace = TRUE), -which(colnames(C5) == "deviance")][c(1,3,2,4)]
}

##################################################################################################################
## for model M6: G-L logistic model

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
  s2 = dnorm(N[2],mean = N[1]+ h*theta[2]*N[1]*(1-(N[1]/theta[3])),sd = 1/sqrt(theta[4]), log = TRUE)
  s3 = 0
  for (i in 2:(n-1)) {
    s3 = s3 + dnorm(N[i+1], mean = mean_M6(N, t,i, h, theta[1],
                                           theta[2], theta[3]), sd = 1/sqrt(theta[4]), log = TRUE)
  }
  z = s1+ s2 + s3
  return(z)
}

# model file  for M6
cat("model{
  N[1] ~ dnorm(N_1, tau)
  N[2] ~ dnorm(N[1] + h*r*N[1]*(1-N[1]/K), tau)
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
  K ~ dunif(min_N, Max)
  tau ~ dgamma(0.01, 0.01)
 }", file = "M6_Grunwald.txt")

## Posterior distribution using JAGS
jagsfit_M6 = jags(data = c('N', 'n', 'h', 'N_1', 'min_N','Max'), inits, params,
                  n.iter = n_iter, model.file = "M6_Grunwald.txt", n.chains = 2, n.thin = 10)
fit_M6 = as.mcmc(jagsfit_M6)
class(fit_M6)

C6 = as.matrix(fit_M6)
filename = paste0("post_M6_GPDD_ID_", GPDD_ID, ".txt")
write.table(C6, filename, row.names = FALSE)

##function for sampling from posterior distribution
draw6 = function(){
  C6[sample(dim(C6)[1], 1, replace = TRUE), -which(colnames(C6) == "deviance")][c(1,3,2,4)]
}


###############################################################################################################

####final step: RJMCMC: posterior model probabilities calculation

library(rjmcmc)
post_model = rjmcmcpost(post.draw = list(draw1, draw2, draw3, draw4, draw5, draw6), 
                        g = list(g, g, g, g, g,g), ginv = list(ginv, ginv, ginv, ginv, ginv, ginv),
                        likelihood = list(log_L_M1,log_L_M2, log_L_M3,log_L_M4, log_L_M5, log_L_M6),
                        param.prior = list(log_prior, log_prior, log_prior, log_prior, log_prior,log_prior),
                        model.prior = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6),chainlength = 2e4, TM.thin = 10 , 
                        save.all = TRUE)


# output storing of RJMCMC in the working directory
filename = paste0("post_mod_prob_real_data_GPDD_ID_", GPDD_ID, ".txt")
write.table(post_model$result$`Posterior Model Probabilities`, filename,row.names = FALSE)
filename = paste0("post_mod_progress_real_data_GPDD_ID_", GPDD_ID, ".txt")
write.table(post_model$progress$prb, filename, row.names = FALSE)
filename = paste0("post_rjmcmc_result_real_data_GPDD_ID_", GPDD_ID, ".txt")
capture.output(post_model$result, file = filename)


###################################################################################################################
#### 2. Barplot of posterior model probabilities ##############

library(ggplot2)
post_mod_prob_116 = as.numeric(Re(read.table("post_mod_prob_real_data_GPDD_ID_116.txt", header = TRUE)[1:6,1]))

Mod_lebels = c("M1", "M2", "M3", "M4", "M5", "M6")
df_116 = data.frame(Model = Mod_lebels,post_mod_prob_116)

Mod_labels = c(expression(M[1], M[2], M[3], M[4], M[5], M[6])) # model names in mathematical symbols

# barplot of posterior model probabilities for GPDD ID 116
p1 = ggplot(data = df_116, aes(x = Model, y = post_mod_prob_116)) +
  labs(x = "Model", y = "Posterior Model Probability", title = "GPDD ID: 116") +
  geom_bar(stat = "identity", fill = "aquamarine4") +
  geom_text(aes(label = round(post_mod_prob_116, 4)), vjust = -0.3, size = 3.5) +
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
ggsave("POST_MOD_PROB_GPDD_ID_116.jpeg", plot = p1, width = 8, height = 6, dpi = 400)



#########################################################################################
#### 3. posterior histogram of model parameters along with posterior density curve for all the models

GPDD_ID = 116
for (m in 1:6) {
  # output reading
  filename = paste0("post_M",m, "_GPDD_ID_", GPDD_ID,".txt")
  C_post = read.table(filename, header = TRUE)
  df = data.frame(C_post[,c(1,4,3,5)])
  
  library(ggplot2)
  
  # Define a light background color
  light_background_color = "white"
  
  # posterior histogram along with posterior density plot for parameter alpha
  h1_116 = ggplot(df, aes(x = alpha, ..density..)) + 
    geom_histogram(aes(y = ..density..),
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
  filename = paste0("POST_hist_alpha_GPDD_ID_", GPDD_ID,"_","M",m, ".jpeg")
  ggsave(filename, plot = h1_116, width = 6, height = 4, dpi = 400) # to save the plot
  
  
  # posterior histogram along with posterior density plot for parameter r
  h2_116 =ggplot(df, aes(x = r, ..density..)) + 
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
  filename = paste0("POST_hist_r_GPDD_ID_", GPDD_ID, "_","M",m,".jpeg")
  ggsave(filename, plot = h2_116, width = 6, height = 4, dpi = 400) # to save the plot
  
  # posterior histogram along with posterior density plot for parameter K
  h3_116 = ggplot(df, aes(x = K, ..density..)) + 
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
  filename = paste0("POST_hist_K_GPDD_ID_", GPDD_ID,"_","M",m, ".jpeg")
  ggsave(filename, plot = h3_116, width = 6, height = 4, dpi = 400) # to save the plot
  
  # posterior histogram along with posterior density plot for parameter tau
  h4_116 = ggplot(df, aes(x = tau, ..density..)) + 
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
  filename = paste0("POST_hist_tau_GPDD_ID_", GPDD_ID,"_","M",m, ".jpeg")
  ggsave(filename, plot = h4_116, width = 6, height = 4, dpi = 400) # to save the plot
  
}


##################################################################################################
#### 4. Parameter estimate (mode) and 95% posterior credible interval

GPDD_ID = 116
# function for calculating mode
mode = function(w){
  density_est = density(w)
  return(density_est$x[which.max(density_est$y)])
}

post_mod_prob_116 = as.numeric(Re(read.table("post_mod_prob_real_data_GPDD_ID_116.txt", header = TRUE)[1:6,1])) #posterior model probabilities
m = which.max(post_mod_prob_116)
filename = paste0("post_M",m, "_GPDD_ID_", GPDD_ID,".txt")
C = read.table(filename, header = TRUE)[,c(1,4,3,5)]
head(C)

param_est_mode = numeric(4)   # array to store the parameter estimate
param_post_CI = matrix(NA, ncol = 2, nrow = 4)  # array to store 95% posterior credible interval (row-wise, for alpha, r, K and tau) 
for (i in 1:4) {
  param_est_mode[i] = mode(C[,i])
  param_post_CI[i,1] = as.numeric(quantile(C[,i],probs = 0.025))
  param_post_CI[i,2] = as.numeric(quantile(C[,i], probs = 0.975))
}
param_est_mode
param_post_CI