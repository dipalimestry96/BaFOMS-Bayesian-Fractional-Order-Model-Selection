## Here we have calculated the 95% Prediction interval by including 20 years future time points
## for that we have generated the parameter values from their joint posterior
## then using that generated the population density trajectory
## we repeated this process for 1000 times and using quantile we have
## calculated the 95% Prediction interval for the population density
## we have taken mode of these generated density values at each time point
## as a fitted value.


#################################################################
######################## Mean functions for the models

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


## for model M2: Caputo-Fabrizio logistic model
# mean function
mean_M2 = function(N_p, N_pp, M_alpha, h, alpha, r, K){
  A = ((1-alpha)/M_alpha + (3*alpha*h)/(2*M_alpha))*r*N_p*(1-(N_p/K))
  B = ((1-alpha)/M_alpha + (alpha*h)/(2*M_alpha))*r*N_pp*(1-(N_pp/K))
  return(N_p + A - B)
}


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



#####################################################################
############################# 95% prediction interval 

library(rgpdd)
data = subset(gpdd_data, MainID == 200)$Population    # population data for GPDD ID 200
N = data
n_previous = length(N)

n = n_previous + 20   # future 20 years
M_alpha = 1
ABC_alpha = 1

t = numeric(length = n)                       # time grid
t[1] = 0                                              #initial time point
h = 0.01                                              # step size
for (j in 1:(n-1)) {
  t[j+1] = t[j] + h
}

# function for calculating mode
mode = function(w){
  density_est = density(w)
  return(density_est$x[which.max(density_est$y)])
}


rep =  1000                                            # number of repetitions for 95% prediction interval calculation
set.seed(123) # for reproducibility of results

setwd("Specify path of local folder for reading and storing outputs")

###########################################################################
######### For M1: prediction interval ##############################
############## posterior distribution of parameters for model M1
filename = paste0("post_M1_GPDD_ID_", "200", ".txt")
C1 = as.matrix(read.table(filename, header = TRUE))[,c(1,4,3,5)]   # posterior values of model (M1) parameters obtained using JAGS

index = sample(1:nrow(C1), rep, replace = TRUE)       # index used for sampling from joint posterior distribution of model parameters
sim_M1 = matrix(data = NA, nrow = rep, ncol = n )      # matrix to store generated population values

######### Generation of population traectories
for (i in 1:rep) {
  alpha = C1[index[i],1]
  r = C1[index[i],2]
  K = C1[index[i],3]
  tau = C1[index[i],4]
  
  sim_M1[i,1] = rnorm(1, mean = N[1], sd = 1/sqrt(tau))
  sim_M1[i,2] = rnorm(1, mean =  N[1] + h*r*N[1]*(1-(N[1]/K)),sd = 1/sqrt(tau) )
  for (s in 3:n) {
    mu = mean_M1(sim_M1[i,s-1],sim_M1[i,s-2], t[s], t[s-1], t[s-2], h, alpha, r,K)
    sim_M1[i,s] =rnorm(1, mu, 1/sqrt(tau))
  }
}

sim_M1[sim_M1<0]=NA
zz_M1 = na.omit(sim_M1)  # omision of NA values
PI_matrix_M1 = matrix(data = NA, n, 3) # matrix with 3 columns for storing the 95% CI bounds and fitted mean profile

## lower and upper bound of 95% CI calculation 
for(y in 1:n){
  PI_matrix_M1[y,1] = quantile(zz_M1[,y], 0.025) #lower bound
  PI_matrix_M1[y,2] = mode(zz_M1[,y])   
  PI_matrix_M1[y,3] = quantile(zz_M1[,y], 0.975) # upper bound
}

filename = paste0("95_per_Prediction_Interval_MATRIX_M1_", "200", ".txt")
write.table(PI_matrix_M1, filename, row.names = FALSE)


filename = paste0("95_per_Prediction_Interval_MATRIX_M", "1", "_", "200", ".txt")  
PI_mat = as.matrix(read.table(filename, header = TRUE)) # 95% CI values and mean profile matrix reading
df = data.frame(t = t*100, PI_mat[,2], PI_mat[,1], PI_mat[,3]) # data frame for plotting
dff = data.frame(t_previous = t[1:n_previous]*100, N) 
dfff = data.frame(t_future = t[(n_previous+1):n]*100, N_future = PI_mat[(n_previous+1):n, 2], 
                  PI_mat[(n_previous+1):n,1], PI_mat[(n_previous+1):n, 3]) 

library(ggplot2)

xlabel = "t (years): 1931-2001"


p = ggplot(df, aes(x = t, y = PI_mat[,2])) +
  geom_ribbon(data = dfff, aes(x = t_future, y = N_future, ymin = PI_mat[(n_previous+1):n,1], ymax = PI_mat[(n_previous+1):n,3], fill = "95% Prediction Interval"), alpha = 0.45) +
  geom_point(data = dff, aes(x = t_previous, y = N, color = "Observed Population Density"), size = 3) +  
  geom_point(data = dfff, aes(x = t_future, y = N_future, color = "Predicted Population Density"), size = 3) +  
  scale_fill_manual(values = "lightblue") +
  scale_color_manual(values = c("Observed Population Density" = "red3", "Predicted Population Density" = "blue4")) +
  labs(x = xlabel, y = "Population Density",
       color = NULL, fill = NULL, title = NULL) +
  theme_classic() +
  ggtitle(bquote(M[.(1)])) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.3, 0.9),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13)
  )


filename = paste0("95_Prediction_Interval_","M", "1", "_", "200(mode)", ".jpeg" )
ggsave(filename, plot = p, width = 8, height = 6, dpi = 600) # to save plot


###########################################################################
######### For M2: prediction interval ##############################

############## posterior distribution of parameters for model M2
filename = paste0("post_M2_GPDD_ID_", "200", ".txt")
C2 = as.matrix(read.table(filename, header = TRUE))[,c(1,4,3,5)]   # posterior values of model (M2) parameters obtained using JAGS

index = sample(1:nrow(C2), rep, replace = TRUE)       # index used for sampling from joint posterior distribution of model parameters
sim_M2 = matrix(data = NA, nrow = rep, ncol = n )      # matrix to store generated population values

######### Generation of population traectories
for (i in 1:rep) {
  alpha = C2[index[i],1]
  r = C2[index[i],2]
  K = C2[index[i],3]
  tau = C2[index[i],4]
  
  sim_M2[i,1] = rnorm(1, mean = N[1], sd = 1/sqrt(tau))
  sim_M2[i,2] = rnorm(1, mean =  N[1] + h*r*N[1]*(1-(N[1]/K)),sd = 1/sqrt(tau) )
  for (s in 3:n) {
    mu = mean_M2(sim_M2[i,s-1],sim_M2[i,s-2],M_alpha, h, alpha, r, K)
    sim_M2[i,s] =rnorm(1, mu, 1/sqrt(tau))
  }
}

sim_M2[sim_M2<0]=NA
zz_M2 = na.omit(sim_M2)  # omision of NA values
PI_matrix_M2 = matrix(data = NA, n, 3) # matrix with 3 columns for storing the 95% CI bounds and fitted mean profile

## lower and upper bound of 95% CI calculation 
for(y in 1:n){
  PI_matrix_M2[y,1] = quantile(zz_M2[,y], 0.025) #lower bound
  PI_matrix_M2[y,2] = mode(zz_M2[,y])   # mean profile
  PI_matrix_M2[y,3] = quantile(zz_M2[,y], 0.975) # upper bound
}

filename = paste0("95_per_Prediction_Interval_MATRIX_M2_", "200", ".txt")
write.table(PI_matrix_M2, filename, row.names = FALSE)


filename = paste0("95_per_Prediction_Interval_MATRIX_M", "2", "_", "200", ".txt")  
PI_mat = as.matrix(read.table(filename, header = TRUE)) # 95% CI values and mean profile matrix reading
df = data.frame(t = t*100, PI_mat[,2], PI_mat[,1], PI_mat[,3]) # data frame for plotting
dff = data.frame(t_previous = t[1:n_previous]*100, N) 
dfff = data.frame(t_future = t[(n_previous+1):n]*100, N_future = PI_mat[(n_previous+1):n, 2], 
                  PI_mat[(n_previous+1):n,1], PI_mat[(n_previous+1):n, 3]) 

library(ggplot2)

xlabel = "t (years): 1931-2001"


p = ggplot(df, aes(x = t, y = PI_mat[,2])) +
  geom_ribbon(data = dfff, aes(x = t_future, y = N_future, ymin = PI_mat[(n_previous+1):n,1], ymax = PI_mat[(n_previous+1):n,3], fill = "95% Prediction Interval"), alpha = 0.45) +
  geom_point(data = dff, aes(x = t_previous, y = N, color = "Observed Population Density"), size = 3) +  
  geom_point(data = dfff, aes(x = t_future, y = N_future, color = "Predicted Population Density"), size = 3) +  
  scale_fill_manual(values = "lightblue") +
  scale_color_manual(values = c("Observed Population Density" = "red3", "Predicted Population Density" = "blue4")) +
  labs(x = xlabel, y = "Population Density",
       color = NULL, fill = NULL, title = NULL) +
  theme_classic() +
  ggtitle(bquote(M[.(2)])) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.3, 0.9),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13)
  )

filename = paste0("95_Prediction_Interval_","M", "2", "_", "200", ".jpeg")
ggsave(filename, plot = p, width = 8, height = 6, dpi = 600) # to save plot

###########################################################################
######### For M3: prediction interval ##############################

############## posterior distribution of parameters for model M3
filename = paste0("post_M3_GPDD_ID_", "200", ".txt")
C3 = as.matrix(read.table(filename, header = TRUE))[,c(1,4,3,5)]   # posterior values of model (M3) parameters obtained using JAGS

index = sample(1:nrow(C3), rep, replace = TRUE)       # index used for sampling from joint posterior distribution of model parameters
sim_M3 = matrix(data = NA, nrow = rep, ncol = n )      # matrix to store generated population values

######### Generation of population traectories
for (i in 1:rep) {
  alpha = C3[index[i],1]
  r = C3[index[i],2]
  K = C3[index[i],3]
  tau = C3[index[i],4]
  
  sim_M3[i,1] = rnorm(1, mean = N[1], sd = 1/sqrt(tau))
  sim_M3[i,2] = rnorm(1, mean =  N[1] + h*r*N[1]*(1-(N[1]/K)),sd = 1/sqrt(tau) )
  for (s in 3:n) {
    mu = mean_M3(sim_M3[i,s-1],sim_M3[i,s-2],t[s], t[s-1], t[s-2]
                 , h, ABC_alpha, alpha, r,K)
    sim_M3[i,s] = rnorm(1, mu, 1/sqrt(tau))
  }
}

sim_M3[sim_M3<0]=NA
zz_M3 = na.omit(sim_M3)  # omision of NA values
PI_matrix_M3 = matrix(data = NA, n, 3) # matrix with 3 columns for storing the 95% CI bounds and fitted mean profile

## lower and upper bound of 95% CI calculation 
for(y in 1:n){
  PI_matrix_M3[y,1] = quantile(zz_M3[,y], 0.025) #lower bound
  PI_matrix_M3[y,2] = mode(zz_M3[,y])   # mean profile
  PI_matrix_M3[y,3] = quantile(zz_M3[,y], 0.975) # upper bound
}

filename = paste0("95_per_Prediction_Interval_MATRIX_M3_", "200", ".txt")
write.table(PI_matrix_M3, filename, row.names = FALSE)


filename = paste0("95_per_Prediction_Interval_MATRIX_M", "3", "_", "200", ".txt")  
PI_mat = as.matrix(read.table(filename, header = TRUE)) # 95% CI values and mean profile matrix reading
df = data.frame(t = t*100, PI_mat[,2], PI_mat[,1], PI_mat[,3]) # data frame for plotting
dff = data.frame(t_previous = t[1:n_previous]*100, N) 
dfff = data.frame(t_future = t[(n_previous+1):n]*100, N_future = PI_mat[(n_previous+1):n, 2], 
                  PI_mat[(n_previous+1):n,1], PI_mat[(n_previous+1):n, 3]) 

library(ggplot2)

xlabel = "t (years): 1931-2001"


p = ggplot(df, aes(x = t, y = PI_mat[,2])) +
  geom_ribbon(data = dfff, aes(x = t_future, y = N_future, ymin = PI_mat[(n_previous+1):n,1], ymax = PI_mat[(n_previous+1):n,3], fill = "95% Prediction Interval"), alpha = 0.45) +
  geom_point(data = dff, aes(x = t_previous, y = N, color = "Observed Population Density"), size = 3) +  
  geom_point(data = dfff, aes(x = t_future, y = N_future, color = "Predicted Population Density"), size = 3) +  
  scale_fill_manual(values = "lightblue") +
  scale_color_manual(values = c("Observed Population Density" = "red3", "Predicted Population Density" = "blue4")) +
  labs(x = xlabel, y = "Population Density",
       color = NULL, fill = NULL, title = NULL) +
  theme_classic() +
  ggtitle(bquote(M[.(3)])) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.3, 0.9),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13)
  )

filename = paste0("95_Prediction_Interval_","M", "3", "_", "200", ".jpeg")
ggsave(filename, plot = p, width = 8, height = 6, dpi = 600) # to save plot

###########################################################################
######### For M4: prediction interval ##############################

############## posterior distribution of parameters for model M4
filename = paste0("post_M4_GPDD_ID_", "200", ".txt")
C4 = as.matrix(read.table(filename, header = TRUE))[,c(1,4,3,5)]   # posterior values of model (M4) parameters obtained using JAGS

index = sample(1:nrow(C4), rep, replace = TRUE)       # index used for sampling from joint posterior distribution of model parameters
sim_M4 = matrix(data = NA, nrow = rep, ncol = n )      # matrix to store generated population values

######### Generation of population traectories
for (i in 1:rep) {
  alpha = C4[index[i],1]
  r = C4[index[i],2]
  K = C4[index[i],3]
  tau = C4[index[i],4]
  
  sim_M4[i,1] = rnorm(1, mean = N[1], sd = 1/sqrt(tau))
  sim_M4[i,2] = rnorm(1, mean =  N[1] + h*r*N[1]*(1-(N[1]/K)),sd = 1/sqrt(tau) )
  for (s in 3:n) {
    mu = mean_M4(sim_M4[i,],t, s-1, h,alpha, r,K)
    sim_M4[i,s] =rnorm(1, mu, 1/sqrt(tau))
  }
}

sim_M4[sim_M4<0]=NA
zz_M4 = na.omit(sim_M4)  # omision of NA values
PI_matrix_M4 = matrix(data = NA, n, 3) # matrix with 3 columns for storing the 95% CI bounds and fitted mean profile

## lower and upper bound of 95% CI calculation 
for(y in 1:n){
  PI_matrix_M4[y,1] = quantile(zz_M4[,y], 0.025) #lower bound
  PI_matrix_M4[y,2] = mode(zz_M4[,y])   # mean profile
  PI_matrix_M4[y,3] = quantile(zz_M4[,y], 0.975) # upper bound
}
filename = paste0("95_per_Prediction_Interval_MATRIX_M4_", "200", ".txt")
write.table(PI_matrix_M4, filename, row.names = FALSE)


filename = paste0("95_per_Prediction_Interval_MATRIX_M", "4", "_", "200", ".txt")  
PI_mat = as.matrix(read.table(filename, header = TRUE)) # 95% CI values and mean profile matrix reading
df = data.frame(t = t*100, PI_mat[,2], PI_mat[,1], PI_mat[,3]) # data frame for plotting
dff = data.frame(t_previous = t[1:n_previous]*100, N) 
dfff = data.frame(t_future = t[(n_previous+1):n]*100, N_future = PI_mat[(n_previous+1):n, 2], 
                  PI_mat[(n_previous+1):n,1], PI_mat[(n_previous+1):n, 3]) 

library(ggplot2)

xlabel = "t (years): 1931-2001"


p = ggplot(df, aes(x = t, y = PI_mat[,2])) +
  geom_ribbon(data = dfff, aes(x = t_future, y = N_future, ymin = PI_mat[(n_previous+1):n,1], ymax = PI_mat[(n_previous+1):n,3], fill = "95% Prediction Interval"), alpha = 0.45) +
  geom_point(data = dff, aes(x = t_previous, y = N, color = "Observed Population Density"), size = 3) +  
  geom_point(data = dfff, aes(x = t_future, y = N_future, color = "Predicted Population Density"), size = 3) +  
  scale_fill_manual(values = "lightblue") +
  scale_color_manual(values = c("Observed Population Density" = "red3", "Predicted Population Density" = "blue4")) +
  labs(x = xlabel, y = "Population Density",
       color = NULL, fill = NULL, title = NULL) +
  theme_classic() +
  ggtitle(bquote(M[.(4)])) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.3, 0.9),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13)
  )
filename = paste0("95_Prediction_Interval_","M", "4", "_", "200", ".jpeg")
ggsave(filename, plot = p, width = 8, height = 6, dpi = 600) # to save plot


###########################################################################
######### For M5: prediction interval ##############################

############## posterior distribution of parameters for model M5
filename = paste0("post_M5_GPDD_ID_", "200", ".txt")
C5 = as.matrix(read.table(filename, header = TRUE))[,c(1,4,3,5)]   # posterior values of model (M5) parameters obtained using JAGS

index = sample(1:nrow(C5), rep, replace = TRUE)       # index used for sampling from joint posterior distribution of model parameters
sim_M5 = matrix(data = NA, nrow = rep, ncol = n )      # matrix to store generated population values

######### Generation of population traectories
for (i in 1:rep) {
  alpha = C5[index[i],1]
  r = C5[index[i],2]
  K = C5[index[i],3]
  tau = C5[index[i],4]
  
  sim_M5[i,1] = rnorm(1, mean = N[1], sd = 1/sqrt(tau))
  sim_M5[i,2] = rnorm(1, mean =  N[1] + h*r*N[1]*(1-(N[1]/K)),sd = 1/sqrt(tau) )
  for (s in 3:n) {
    mu = mean_M5(sim_M5[i,],t,s-1,M_alpha,h,alpha, r, K) 
    sim_M5[i,s] =rnorm(1, mu, 1/sqrt(tau))
  }
}

sim_M5[sim_M5<0]=NA
zz_M5 = na.omit(sim_M5)  # omision of NA values
PI_matrix_M5 = matrix(data = NA, n, 3) # matrix with 3 columns for storing the 95% CI bounds and fitted mean profile

## lower and upper bound of 95% CI calculation 
for(y in 1:n){
  PI_matrix_M5[y,1] = quantile(zz_M5[,y], 0.025) #lower bound
  PI_matrix_M5[y,2] = mode(zz_M5[,y])   # mean profile
  PI_matrix_M5[y,3] = quantile(zz_M5[,y], 0.975) # upper bound
}
filename = paste0("95_per_Prediction_Interval_MATRIX_M5_", "200", ".txt")
write.table(PI_matrix_M5, filename, row.names = FALSE)


filename = paste0("95_per_Prediction_Interval_MATRIX_M", "5", "_", "200", ".txt")  
PI_mat = as.matrix(read.table(filename, header = TRUE)) # 95% CI values and mean profile matrix reading
df = data.frame(t = t*100, PI_mat[,2], PI_mat[,1], PI_mat[,3]) # data frame for plotting
dff = data.frame(t_previous = t[1:n_previous]*100, N) 
dfff = data.frame(t_future = t[(n_previous+1):n]*100, N_future = PI_mat[(n_previous+1):n, 2], 
                  PI_mat[(n_previous+1):n,1], PI_mat[(n_previous+1):n, 3]) 

library(ggplot2)

xlabel = "t (years): 1931-2001"


p = ggplot(df, aes(x = t, y = PI_mat[,2])) +
  geom_ribbon(data = dfff, aes(x = t_future, y = N_future, ymin = PI_mat[(n_previous+1):n,1], ymax = PI_mat[(n_previous+1):n,3], fill = "95% Prediction Interval"), alpha = 0.45) +
  geom_point(data = dff, aes(x = t_previous, y = N, color = "Observed Population Density"), size = 3) +  
  geom_point(data = dfff, aes(x = t_future, y = N_future, color = "Predicted Population Density"), size = 3) +  
  scale_fill_manual(values = "lightblue") +
  scale_color_manual(values = c("Observed Population Density" = "red3", "Predicted Population Density" = "blue4")) +
  labs(x = xlabel, y = "Population Density",
       color = NULL, fill = NULL, title = NULL) +
  theme_classic() +
  ggtitle(bquote(M[.(5)])) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.3, 0.9),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13)
  )

filename = paste0("95_Prediction_Interval_","M", "5", "_", "200", ".jpeg")
ggsave(filename, plot = p, width = 8, height = 6, dpi = 600) # to save plot


########################################################################
######### For M6: prediction interval ###########################

############## posterior distribution of parameters for model M6
filename = paste0("post_M6_GPDD_ID_", "200", ".txt")
C6 = as.matrix(read.table(filename, header = TRUE))[,c(1,4,3,5)]   # posterior values of model (M6) parameters obtained using JAGS

index = sample(1:nrow(C6), rep, replace = TRUE)       # index used for sampling from joint posterior distribution of model parameters
sim_M6 = matrix(data = NA, nrow = rep, ncol = n )      # matrix to store generated population values

######### Generation of population traectories 
for (i in 1:rep) {
  alpha = C6[index[i],1]
  r = C6[index[i],2]
  K = C6[index[i],3]
  tau = C6[index[i],4]
  
  sim_M6[i,1] = rnorm(1, mean = N[1], sd = 1/sqrt(tau))
  sim_M6[i,2] = rnorm(1, mean =  N[1] + h*r*N[1]*(1-(N[1]/K)),sd = 1/sqrt(tau) )
  for (s in 3:n) {
    mu = mean_M6(sim_M6[i,],t,s-1,h,alpha, r,K)
    sim_M6[i,s] = rnorm(1, mu, 1/sqrt(tau))
  }
}

sim_M6[sim_M6<0]=NA
zz_M6 = na.omit(sim_M6)  # omision of NA values
PI_matrix_M6 = matrix(data = NA, n, 3) # matrix with 3 columns for storing the 95% CI bounds and fitted mean profile

## lower and upper bound of 95% CI calculation 
for(y in 1:n){
  PI_matrix_M6[y,1] = quantile(zz_M6[,y], 0.025) #lower bound
  PI_matrix_M6[y,2] = mode(zz_M6[,y])   # mean profile
  PI_matrix_M6[y,3] = quantile(zz_M6[,y], 0.975) # upper bound
}
filename = paste0("95_per_Prediction_Interval_MATRIX_M6_", "200", ".txt")
write.table(PI_matrix_M6, filename, row.names = FALSE)


filename = paste0("95_per_Prediction_Interval_MATRIX_M", "6", "_", "200", ".txt")  
PI_mat = as.matrix(read.table(filename, header = TRUE)) # 95% CI values and mean profile matrix reading
df = data.frame(t = t*100, PI_mat[,2], PI_mat[,1], PI_mat[,3]) # data frame for plotting
dff = data.frame(t_previous = t[1:n_previous]*100, N) 
dfff = data.frame(t_future = t[(n_previous+1):n]*100, N_future = PI_mat[(n_previous+1):n, 2], 
                  PI_mat[(n_previous+1):n,1], PI_mat[(n_previous+1):n, 3]) 

library(ggplot2)

xlabel = "t (years): 1931-2001"


p = ggplot(df, aes(x = t, y = PI_mat[,2])) +
  geom_ribbon(data = dfff, aes(x = t_future, y = N_future, ymin = PI_mat[(n_previous+1):n,1], ymax = PI_mat[(n_previous+1):n,3], fill = "95% Prediction Interval"), alpha = 0.45) +
  geom_point(data = dff, aes(x = t_previous, y = N, color = "Observed Population Density"), size = 3) +  
  geom_point(data = dfff, aes(x = t_future, y = N_future, color = "Predicted Population Density"), size = 3) +  
  scale_fill_manual(values = "lightblue") +
  scale_color_manual(values = c("Observed Population Density" = "red3", "Predicted Population Density" = "blue4")) +
  labs(x = xlabel, y = "Population Density",
       color = NULL, fill = NULL, title = NULL) +
  theme_classic() +
  ggtitle(bquote(M[.(6)])) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.3, 0.9),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13)
  )

filename = paste0("95_Prediction_Interval_","M", "6", "_", "200", ".jpeg" )
ggsave(filename, plot = p, width = 8, height = 6, dpi = 600) # to save plot
