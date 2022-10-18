library(ggplot2)
library(R2jags)

set.seed(1234)

#set parameters
I = 140
K = 5
T = 14

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

#for k = 3, implement sigmoidal function with mean 0
k = 2
# height = rgamma(1,1,1)
# scale = rgamma(1,1,1)
# center = runif(1,0,T)
height = 2
scale = 1/3
center = T/2
for(t in 1:tsize[k]){
  d_kt[k,t] = d[k] - height +  2*height/(1+exp(-scale*(short_year_ki[k,t] - center)))
}
plot(short_year_ki[k,1:tsize[k]],d_kt[k,1:tsize[k]])

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
# 
# ##################
# # BNMA Sigmoidal #
# ##################
# # given d_kt
# 
# pi = 0.5
# m_mu = 0
# prec_mu = 1
# 
# jags_data <- list("I","n_ik","y_ik",
#                   "m_mu","prec_mu",
#                   "t2","t1","K",
#                   "short_year_ki",
#                   "tsize","TS",
#                   "lts_ind", "T",
#                   "pi", "d_kt")
# 
# #note which params to save
# jags_params <- c("z_k","d",
#                  "min_k", "height_k",
#                  "scale_k", "center_k",
#                  "psi")
# 
# #define inititailization values
# jags_inits <- function(){
#   list("sd" = runif(1,0,5)
#   )
# }
# 
# jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_Sigmoidal_given_dkt.bug"
# 
# #fit model
# 
# jags_fit_sig <- R2jags::jags(data = jags_data, 
#                              inits = jags_inits,
#                              parameters.to.save = jags_params,
#                              n.chains = 3, n.iter = 9000,n.burnin = 1000,
#                              model.file =  jags_file
# )
# 
# print(jags_fit_sig)
# 
# post_d = jags_fit_sig$BUGSoutput$sims.list$d
# apply(post_d,2,mean)
# 
# post_z = jags_fit_sig$BUGSoutput$sims.list$z_k
# apply(post_z,2,mean)
# 
# post_height = jags_fit_sig$BUGSoutput$sims.list$height_k
# post_scale = jags_fit_sig$BUGSoutput$sims.list$scale_k
# post_center = jags_fit_sig$BUGSoutput$sims.list$center_k
# 
# post_psi = jags_fit_sig$BUGSoutput$sims.list$psi
# 
# N = length(post_psi)
# 
# plot_years = NULL
# plot_pred_mu = NULL
# plot_pred_low = NULL
# plot_pred_high = NULL
# plot_k = NULL
# 
# find_quant = function(x){
#   return(quantile(x,c(0.025,0.975)))
# }
# 
# for(k in 2:K){
#   pred_x = 0:(10*T)/10
#   
#   
#   pred_mu_mat = matrix(0,nrow = N, ncol = length(pred_x))
#   for(n in 1:N){
#     pred_mu_mat[n,] = post_z[n,k]*((post_d[n,k] - post_height[n,k]) +  2*( post_height[n,k])/(1+exp(-post_scale[n,k]*(pred_x - post_center[n,k])))) + (1-post_z[n,k])*rep(post_d[n,k],length(pred_x))
#   }
#   # pred_mu =   post_min[k] +  post_height[k]/(1+exp(-post_scale[k]*(pred_x - post_center[k])))
#   pred_mu =   apply(pred_mu_mat,2,mean)
#   quant_mat = apply(pred_mu_mat,2,find_quant)
#   
#   
#   # 
#   #   pred_low = pred_mu - 1.95 * sqrt(diag(pred_Sigma))
#   #   pred_high = pred_mu + 1.95 * sqrt(diag(pred_Sigma))
#   
#   plot_years = c(plot_years,pred_x)
#   plot_pred_mu = c(plot_pred_mu,pred_mu)
#   plot_pred_low = c(plot_pred_low,quant_mat[1,])
#   plot_pred_high = c(plot_pred_high,quant_mat[2,])
#   plot_k = c(plot_k,rep(treatments[k],length(pred_x)))
# }
# 
# df = data.frame("Years" = plot_years,
#                 "Mean" = plot_pred_mu,
#                 "Low" = plot_pred_low,
#                 "High" = plot_pred_high,
#                 "K" = as.factor(plot_k))
# 
# ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
#   geom_line()
# 

##################
# BNMA Sigmoidal #
##################
# given delta

pi = 0.5
m_mu = 0
prec_mu = 1

jags_data <- list("I","n_ik","y_ik",
                  "m_mu","prec_mu",
                  "t2","t1","K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind", "T",
                  "pi", "delta")

#note which params to save
jags_params <- c("z_k","d",
                 "min_k", "height_k",
                 "scale_k", "center_k",
                 "psi")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/tBNMA/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_Sigmoidal_given_delta.bug"

#fit model

jags_fit_sig <- R2jags::jags(data = jags_data, 
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)

print(jags_fit_sig)

post_d = jags_fit_sig$BUGSoutput$sims.list$d
apply(post_d,2,mean)

post_z = jags_fit_sig$BUGSoutput$sims.list$z_k
apply(post_z,2,mean)

post_height = jags_fit_sig$BUGSoutput$sims.list$height_k
post_scale = jags_fit_sig$BUGSoutput$sims.list$scale_k
post_center = jags_fit_sig$BUGSoutput$sims.list$center_k

post_psi = jags_fit_sig$BUGSoutput$sims.list$psi

N = length(post_psi)

plot_years = NULL
plot_pred_mu = NULL
plot_pred_low = NULL
plot_pred_high = NULL
plot_k = NULL

find_quant = function(x){
  return(quantile(x,c(0.025,0.975)))
}

for(k in 2:K){
  pred_x = 0:(10*T)/10
  
  
  pred_mu_mat = matrix(0,nrow = N, ncol = length(pred_x))
  for(n in 1:N){
    pred_mu_mat[n,] = post_z[n,k]*((post_d[n,k] - post_height[n,k]) +  2*( post_height[n,k])/(1+exp(-post_scale[n,k]*(pred_x - post_center[n,k])))) + (1-post_z[n,k])*rep(post_d[n,k],length(pred_x))
  }
  # pred_mu =   post_min[k] +  post_height[k]/(1+exp(-post_scale[k]*(pred_x - post_center[k])))
  pred_mu =   apply(pred_mu_mat,2,mean)
  quant_mat = apply(pred_mu_mat,2,find_quant)
  
  
  # 
  #   pred_low = pred_mu - 1.95 * sqrt(diag(pred_Sigma))
  #   pred_high = pred_mu + 1.95 * sqrt(diag(pred_Sigma))
  
  plot_years = c(plot_years,pred_x)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,quant_mat[1,])
  plot_pred_high = c(plot_pred_high,quant_mat[2,])
  plot_k = c(plot_k,rep(treatments[k],length(pred_x)))
}

df = data.frame("Years" = plot_years,
                "Mean" = plot_pred_mu,
                "Low" = plot_pred_low,
                "High" = plot_pred_high,
                "K" = as.factor(plot_k))

ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
  geom_line()

