#### simulation study: data generated using the M4 model with 5 different sets of parameters

setwd("set here path for saving outputs")

## Needed packages
library(R2jags)
library(coda)
library(madness)
library(rjmcmc)

seeds = 141:145 # used seeds for reproducibility of results
set = matrix(NA, ncol = 4, nrow = length(seeds)) # matrix for saving the different sets values
set[1,] = c(0.2, 0.8, 30000, 300) # set 1
set[2,] = c(0.3, 0.6, 20000, 200) # set 2
set[3,] = c(0.25, 0.7, 30000, 350) # set 3
set[4,] = c(0.35, 0.65, 30000, 200) # set 4
set[5,] = c(0.4, 0.45, 20000, 250) # set 5

post_prob_M4 = matrix(NA, length(seeds), 6) # matrix for saving the posterior model probabilities set wise

## data simulation step from the model M4
for (w in 1:length(seeds)) {
  # simulated data generation for analysis using different sets
  set.seed(seeds[w])
  S = set[w,]
  alpha_sim = S[1]
  r_sim = S[2]
  K_sim = S[3]
  sigma_sim = S[4]
  tau_sim = 1/sigma_sim^2
  n = 40
  t0 = 0
  tn = 20
  
  t = seq(t0, tn, length.out = n) # time grid
  h = t[2]-t[1]
  
  M_alpha = 1
  ABC_alpha = 1
  
  ##Mean function
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
  
  
  N = numeric(n)
  mu_0 = 1000 # mean population size for N[1]
  N[1] = rnorm(1,mean = mu_0, sd = sigma_sim)
  mu_star = N[1] + h*r_sim*N[1]*(1-N[1]/K_sim)
  N[2] = rnorm(1, mean = mu_star, sd = sigma_sim )
  for (i in 2:(n-1)) {
    mean = mean_M4(N, t,i,h,alpha_sim, r_sim, K_sim)
    N[i+1] = rnorm(1, mean = mean, sd = sigma_sim )
    
  }
  
  par(mfrow = c(1,1))
  plot(t, N, type = "b", col= "red", pch = 16, main = "Simulated data from M4")
  
  
  ##################################################################################################
  ##################################################################################################
  
  
  # log-prior for all model parameters
  log_prior = function(theta){
    sum(dunif(theta[1], min = 0, max = 1, log = TRUE),
        dgamma(theta[2], 2, 1/2, log = TRUE) ,
        dnorm(theta[3], mean = 30000, sd = 100, log = TRUE),
        dgamma(theta[4], 0.01, 0.01, log = TRUE ),
        dnorm(theta[5], 5000, sqrt(1000), log = TRUE))
  }
  
  ##initial values
  inits = function(){
    list(alpha = 0.7, r = 0.7, K = 1000, tau = 2, mu_0 = 1000)
  }
  
  params = c("alpha", "r", "K", "tau", "mu_0")
  
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
  
  jagsfit_M1 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim'), inits, params,
                    n.iter = n_iter, model.file = "M1_Caputo.txt", n.chains = 3, n.thin = 1)
  fit_M1 = as.mcmc(jagsfit_M1) # mcmc after burnin period=iter/2
  
  
  C1 = as.matrix(fit_M1)
  write.table(C1, "M1_output.txt", row.names = FALSE)
  
  ##sampling function from posterior
  
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
  
  jagsfit_M2 = jags(data = c('N', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim','M_alpha'), inits, params,
                    n.iter = n_iter, model.file = "M2_Caputo_Fabrizio.txt", n.chains = 3, n.thin = 1)
  fit_M2 = as.mcmc(jagsfit_M2)
  class(fit_M2)
  
  C2 = as.matrix(fit_M2)
  write.table(C2, "M2_output.txt", row.names = FALSE)
  
  ##sampling function from posterior
  
  draw2 = function(){
    C2[sample(dim(C2)[1], 1, replace = TRUE), -which(colnames(C2) == "deviance")][c(1,4,2,5,3)]
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
  
  
  jagsfit_M3 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim','ABC_alpha'), inits, params,
                    n.iter = n_iter, model.file = "M3_ABC.txt", n.chains = 3, n.thin = 1)
  fit_M3 = as.mcmc(jagsfit_M3)
  class(fit_M3)
  
  C3 = as.matrix(fit_M3)
  write.table(C3, "M3_output.txt", row.names = FALSE)
  
  ##sampling function from posterior
  
  draw3 = function(){
    C3[sample(dim(C3)[1], 1, replace = TRUE), -which(colnames(C3) == "deviance")][c(1,4,2,5,3)]
  }
  
  ############################################################################################################
  
  ## for model M4: RL logistic model
  
  ##Mean function
  
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
  
  jagsfit_M4 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim'), inits, params,
                    n.iter = n_iter, model.file = "M4_RL.txt", n.chains = 3, n.thin = 1)
  fit_M4 = as.mcmc(jagsfit_M4)
  class(fit_M4)
  
  C4 = as.matrix(fit_M4)
  write.table(C4, "M4_output.txt", row.names = FALSE)
  
  ##sampling function from posterior
  
  draw4 = function(){
    C4[sample(dim(C4)[1], 1, replace = TRUE), -which(colnames(C4) == "deviance")][c(1,4,2,5,3)]
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
  
  
  
  jagsfit_M5 = jags(data = c('N', 't', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim', 'M_alpha'), inits, params,
                    n.iter = n_iter, model.file = "M5_RL_CF.txt", n.chains = 3, n.thin = 1)
  fit_M5 = as.mcmc(jagsfit_M5)
  class(fit_M5)
  C5 = as.matrix(fit_M5)
  write.table(C5, "M5_output.txt", row.names = FALSE)
  
  ##sampling function from posterior
  
  draw5 = function(){
    C5[sample(dim(C5)[1], 1, replace = TRUE), -which(colnames(C5) == "deviance")][c(1,4,2,5,3)]
  }
  
  
  ###################################################################################################################
  
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
  
  
  jagsfit_M6 = jags(data = c('N', 'n', 'sigma_sim', 'h', 'r_sim', 'K_sim'), inits, params,
                    n.iter = n_iter, model.file = "M6_Grunwald.txt", n.chains = 3, n.thin = 1)
  fit_M6 = as.mcmc(jagsfit_M6)
  class(fit_M6)
  
  C6 = as.matrix(fit_M6)
  write.table(C6, "M6_output.txt", row.names = FALSE)
  
  ##sampling function from posterior
  
  draw6 = function(){
    C6[sample(dim(C6)[1], 1, replace = TRUE), -which(colnames(C6) == "deviance")][c(1,4,2,5,3)]
  }
  
  ###############################################################################################################
  
  #### final step for multimodel inference: RJMCMC
  
  post_model = rjmcmcpost(post.draw = list(draw1, draw2, draw3, draw4, draw5, draw6), 
                          g = list(g, g, g, g, g,g), ginv = list(ginv, ginv, ginv, ginv, ginv, ginv),
                          likelihood = list(log_L_M1,log_L_M2, log_L_M3,log_L_M4, log_L_M5, log_L_M6),
                          param.prior = list(log_prior, log_prior, log_prior, log_prior, log_prior,log_prior),
                          model.prior = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6),chainlength = 5e4, TM.thin = 10 , 
                          save.all = TRUE)
  post_prob_M4[w,] = post_model$result$`Posterior Model Probabilities`
  print(post_prob_M4[w,])
  filename = paste0("post_mod_prob_M4_set_", w, ".txt")
  write.table(post_prob_M4[w,], filename, row.names = FALSE)
  
}
write.table(post_prob_M4, "posterior_prob_mat_M4.txt",row.names = FALSE) # writing output file in given local folder on system
View(round(Re(post_prob_M4), digits = 5))


post_mod_prob = read.table("posterior_prob_mat_M4.txt", header = TRUE)
post_mod_prob = as.matrix(post_mod_prob)
round(Re(post_mod_prob),digits = 4)
