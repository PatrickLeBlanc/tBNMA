library(R2jags)

#set parameters
I = 150
K = 5
T = 50

#make some data structures
true_model = NULL
ran_model = NULL
d_rmse = NULL
dkt_rmse = NULL
dlast_rmse = NULL
time_vec = NULL
I_vec = NULL
K_vec = NULL
T_vec = NULL

##############
# Generate BNMA #
##############


set.seed(1234)

#list studies and treatments
studies = 1:I
treatments = 1:K

#simulate years
years = sample(1:T,I,replace = TRUE)

#simulate treatment 1 and treatment 2
t1 = sample(1:K,I,replace = TRUE)
t2 = rep(0,I)
for(i in 1:I){
  #treatment 2 can't be treatment 1
  pos_treat = 1:K
  pos_treat = pos_treat[-t1[i]]
  t2[i] = sample(pos_treat,1)
  
  #want t1 < t2 - if not the case, reorder
  if(t1[i] >= t2[i]){
    temp = t2[i]
    t2[i] = t1[i]
    t1[i] = temp
  }
}

#simulate sd from unif(0,1)
sd = runif(1,0,0.5)
tau = 1/sd^2

#simulate d according to standard normal
d = rep(0,K)
for(k in 2:K){
  d[k] = rnorm(1,0,1)
}

#simulate d_kt?

#find how many datapoints for each treatment
tsize = rep(0,K)
for(k in 1:K){
  for(i in 1:I){
    if(t1[i] == k | t2[i] == k){
      tsize[k] = tsize[k] + 1
    } 
  }
}

TS = max(tsize)

#which years each study appears, indexed by which study they appear in
long_year_ki = matrix(NA,nrow = K, ncol = I)
for(k in 1:K){
  for(i in 1:I){
    if(t1[i] == k | t2[i] == k){
      long_year_ki[k,i] = years[i]
    } 
  }
}

short_year_ki = matrix(NA,nrow = K, ncol = max(tsize))
for(k in 1:K){
  temp = long_year_ki[k,]
  temp = temp[!is.na(temp)]
  for(t in 1:tsize[k]){
    short_year_ki[k,t] = temp[t]
  }
}

#need a map leading from study to short index for each treatment
lts_ind = matrix(NA,nrow = K, ncol = I)
for(k in 1:K){
  sum = 0
  for(i in 1:I){
    if(is.na(long_year_ki[k,i]) == FALSE){
      sum = sum + 1
      lts_ind[k,i] = sum
    }
  }
}

d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    d_kt[k,t] = d[k]
  }
}

#simulate delta for each study
#implicitly assume that delta for t1 is zero and that there are only 2 treatments
delta = matrix(0,nrow = I, ncol = 2)
for(i in 1:I){
  mu = d_kt[t2[i],lts_ind[t2[i],i]] - d_kt[t1[i],lts_ind[t1[i],i]]
  delta[i,2] = rnorm(1,mu,sd = sd)
}

#simulate trial specific baseline effects
mu_i = rnorm(I,0,sd = 1)

#simulate number of pateitns
lambdas = c(50,250,1000)
n_ik = matrix(0,nrow = I, ncol = 2)
for(i in 1:I){
  draw = rpois(1,sample(lambdas,1,prob = c(7/16,7/16,2/16)))
  for(k in 1:2){
    n_ik[i,k] = draw
  }
}

#simulate y
y_ik = matrix(0,nrow = I, ncol = 2)
for(i in 1:I){
  
  for(k in 1:2){
    p = exp(mu_i[i] + delta[i,k])/(1+exp(mu_i[i] + delta[i,k]))
    y_ik[i,k] = rbinom(1,n_ik[i,k],p)
  }
  
}

d_last = d

########
# BNMA #
########
m_mu = 0
prec_mu = 1

#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K")

#note which params to save
jags_params <- c("d","sd",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm.bug"

#fit model

begin = Sys.time()
jags_fit <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")

post_d = jags_fit$BUGSoutput$mean$d
post_d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    post_d_kt[k,t] = post_d[k]
  }
}
post_d_last = post_d

true_model = c(true_model,"BNMA")
ran_model = c(ran_model,"BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)

##################
# BNMA Sigmoidal #
##################

pi = 0.5
m_mu = 0
prec_mu = 1

jags_data <- list("I","n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind", "T",
                  "pi")

#note which params to save
jags_params <- c("z_k","d", "d_kt",
                 "height_k",
                 "scale_k", "center_k",
                 "psi")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_Sigmoidal_v2.bug"

#fit model
begin = Sys.time()
jags_fit_sig <- R2jags::jags(data = jags_data, 
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")


post_d = jags_fit_sig$BUGSoutput$mean$d
post_d_kt = jags_fit_sig$BUGSoutput$mean$d_kt

post_z = jags_fit_sig$BUGSoutput$mean$z_k
post_height = jags_fit_sig$BUGSoutput$mean$height_k
post_scale = jags_fit_sig$BUGSoutput$mean$scale_k
post_center = jags_fit_sig$BUGSoutput$mean$center_k

post_psi = jags_fit_sig$BUGSoutput$mean$psi


pred_x = T
post_d_last = post_z*((post_d - post_height) +  2*(post_d + post_height)/(1+exp(-post_scale*(pred_x - post_center)))) + (1-post_z)*rep(post_d,length(pred_x))

true_model = c(true_model,"BNMA")
ran_model = c(ran_model,"Sig-BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)

###################
# Meta-regression #
###################

m_mu = 0
prec_mu = 1

#impelment a time model
#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K",
                  "years")

#note which params to save
jags_params <- c("b","d","sd")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time.bug"

#fit model
set.seed(123)

begin = Sys.time()
jags_fit_meta = R2jags::jags(data = jags_data, 
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")


post_b = jags_fit_meta$BUGSoutput$mean$b
post_d = jags_fit_meta$BUGSoutput$mean$d

post_d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    post_d_kt[k,t] = post_d[k] + short_year_ki[k,t]*post_b[k]
  }
}

post_d_last = post_d + T*post_b

true_model = c(true_model,"BNMA")
ran_model = c(ran_model,"Meta-BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)

###########
# GP BNMA #
###########

m_mu = 0
prec_mu = 1

#JAGS Now

#save data for jags
jags_data <- list("I", "n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind")

#note which params to save
jags_params <- c("d","d_kt","sd",
                 "phi","rho","psi",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_GP.bug"

#fit model
begin = Sys.time()
jags_fit_gp <- R2jags::jags(data = jags_data,
                              inits = jags_inits,
                              parameters.to.save = jags_params,
                              n.chains = 3, n.iter = 9000,n.burnin = 1000,
                              model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")


post_d = jags_fit_gp$BUGSoutput$mean$d
post_d_kt = jags_fit_gp$BUGSoutput$mean$d_kt



post_rho = jags_fit_gp$BUGSoutput$mean$rho
post_phi = jags_fit_gp$BUGSoutput$mean$phi
post_psi = jags_fit_gp$BUGSoutput$mean$psi
post_sd = jags_fit_gp$BUGSoutput$mean$sd

post_d_last = rep(0,K)
for(k in 2:K){
  obs_x =  short_year_ki[k,1:tsize[k]]
  y = post_d_kt[k,1:tsize[k]]
  
  pred_x = T
  
  mu1 = rep(post_d[k],length(pred_x))
  mu2 = rep(post_d[k],length(obs_x))
  
  
  #find covariance
  tot_years = c(pred_x,obs_x)
  Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
  for(i in 1:nrow(Sigma)){
    
    if(i <= length(pred_x)){
      Sigma[i,i] = post_phi[k]^2
    } else {
      Sigma[i,i] = post_psi^2 + post_phi[k]^2
    }
    
    if(i < nrow(Sigma)){
      for(j in (i+1):ncol(Sigma)){
        Sigma[i,j] = post_phi[k]^2*exp(-post_rho[k]*(tot_years[i] - tot_years[j])^2 )
        Sigma[j,i] = Sigma[i,j]
      }
    }
    
  }
  
  one_ind = 1:length(pred_x)
  two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))
  Sigma11 = Sigma[one_ind,one_ind]
  Sigma12 = Sigma[one_ind,two_ind]
  Sigma22 = Sigma[two_ind,two_ind]
  Sigma22_inv = solve(Sigma22)
  Sigma21 = Sigma[two_ind,one_ind]
  
  post_d_last[k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
}

true_model = c(true_model,"BNMA")
ran_model = c(ran_model,"GP-BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)

######################
# Generate Sigmoidal #
######################

#for k = 3, put sigmoidal effects with d[3] as mean
height = rgamma(1,1,1)
scale = rgamma(1,1,1)
center = runif(1,0,T)

k = 3
for(t in 1:tsize[k]){
  time = short_year_ki[k,t]
  d_kt[k,t] = (d[k] - height) +  2*(height)/(1+exp(-scale*(time - center)))
}


########
# BNMA #
########
m_mu = 0
prec_mu = 1

#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K")

#note which params to save
jags_params <- c("d","sd",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm.bug"

#fit model

begin = Sys.time()
jags_fit <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")

post_d = jags_fit$BUGSoutput$mean$d
post_d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    post_d_kt[k,t] = post_d[k]
  }
}
post_d_last = post_d

true_model = c(true_model,"Sig-BNMA")
ran_model = c(ran_model,"BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)

##################
# BNMA Sigmoidal #
##################

pi = 0.5
m_mu = 0
prec_mu = 1

jags_data <- list("I","n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind", "T",
                  "pi")

#note which params to save
jags_params <- c("z_k","d", "d_kt",
                 "height_k",
                 "scale_k", "center_k",
                 "psi")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_Sigmoidal_v2.bug"

#fit model
begin = Sys.time()
jags_fit_sig <- R2jags::jags(data = jags_data, 
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")


post_d = jags_fit_sig$BUGSoutput$mean$d
post_d_kt = jags_fit_sig$BUGSoutput$mean$d_kt

post_z = jags_fit_sig$BUGSoutput$mean$z_k
post_height = jags_fit_sig$BUGSoutput$mean$height_k
post_scale = jags_fit_sig$BUGSoutput$mean$scale_k
post_center = jags_fit_sig$BUGSoutput$mean$center_k

post_psi = jags_fit_sig$BUGSoutput$mean$psi


pred_x = T
post_d_last = post_z*((post_d - post_height) +  2*(post_d + post_height)/(1+exp(-post_scale*(pred_x - post_center)))) + (1-post_z)*rep(post_d,length(pred_x))

true_model = c(true_model,"Sig-BNMA")
ran_model = c(ran_model,"Sig-BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)

###################
# Meta-regression #
###################

m_mu = 0
prec_mu = 1

#impelment a time model
#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K",
                  "years")

#note which params to save
jags_params <- c("b","d","sd")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time.bug"

#fit model
set.seed(123)

begin = Sys.time()
jags_fit_meta = R2jags::jags(data = jags_data, 
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")


post_b = jags_fit_meta$BUGSoutput$mean$b
post_d = jags_fit_meta$BUGSoutput$mean$d

post_d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    post_d_kt[k,t] = post_d[k] + short_year_ki[k,t]*post_b[k]
  }
}

post_d_last = post_d + T*post_b

true_model = c(true_model,"Sig-BNMA")
ran_model = c(ran_model,"Meta-BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)

###########
# GP BNMA #
###########

m_mu = 0
prec_mu = 1

#JAGS Now

#save data for jags
jags_data <- list("I", "n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind")

#note which params to save
jags_params <- c("d","d_kt","sd",
                 "phi","rho","psi",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_GP.bug"

#fit model
begin = Sys.time()
jags_fit_gp <- R2jags::jags(data = jags_data,
                            inits = jags_inits,
                            parameters.to.save = jags_params,
                            n.chains = 3, n.iter = 9000,n.burnin = 1000,
                            model.file =  jags_file
)
end = Sys.time()
time = as.numeric(end - begin,units = "secs")


post_d = jags_fit_gp$BUGSoutput$mean$d
post_d_kt = jags_fit_gp$BUGSoutput$mean$d_kt



post_rho = jags_fit_gp$BUGSoutput$mean$rho
post_phi = jags_fit_gp$BUGSoutput$mean$phi
post_psi = jags_fit_gp$BUGSoutput$mean$psi
post_sd = jags_fit_gp$BUGSoutput$mean$sd

post_d_last = rep(0,K)
for(k in 2:K){
  obs_x =  short_year_ki[k,1:tsize[k]]
  y = post_d_kt[k,1:tsize[k]]
  
  pred_x = T
  
  mu1 = rep(post_d[k],length(pred_x))
  mu2 = rep(post_d[k],length(obs_x))
  
  
  #find covariance
  tot_years = c(pred_x,obs_x)
  Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
  for(i in 1:nrow(Sigma)){
    
    if(i <= length(pred_x)){
      Sigma[i,i] = post_phi[k]^2
    } else {
      Sigma[i,i] = post_psi^2 + post_phi[k]^2
    }
    
    if(i < nrow(Sigma)){
      for(j in (i+1):ncol(Sigma)){
        Sigma[i,j] = post_phi[k]^2*exp(-post_rho[k]*(tot_years[i] - tot_years[j])^2 )
        Sigma[j,i] = Sigma[i,j]
      }
    }
    
  }
  
  one_ind = 1:length(pred_x)
  two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))
  Sigma11 = Sigma[one_ind,one_ind]
  Sigma12 = Sigma[one_ind,two_ind]
  Sigma22 = Sigma[two_ind,two_ind]
  Sigma22_inv = solve(Sigma22)
  Sigma21 = Sigma[two_ind,one_ind]
  
  post_d_last[k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
}

true_model = c(true_model,"Sig-BNMA")
ran_model = c(ran_model,"GP-BNMA")
d_rmse = c(d_rmse,sqrt(mean((post_d - d)^2)))
dkt_rmse = c(dkt_rmse,sqrt(mean((post_d_kt - d_kt)^2)))
dlast_rmse = c(dlast_rmse,sqrt(mean((post_d_last - d_last)^2)))
time_vec = c(time_vec,time)
I_vec = c(I_vec,I)
K_vec = c(K_vec,K)
T_vec = c(T_vec,T)