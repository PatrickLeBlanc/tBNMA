library(ggplot2)
library(R2jags)

set.seed(1234)

#set parameters
I = 150
K = 5
T = 50

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

# #for k = 2, implement sigmoidal function with mean 0
# k = 2
# # height = rgamma(1,1,1)
# # scale = rgamma(1,1,1)
# # center = runif(1,0,T)
# height = 2
# scale = 1/3
# center = T/2
# for(t in 1:tsize[k]){
#   d_kt[k,t] = d[k] - height +  2*height/(1+exp(-scale*(short_year_ki[k,t] - center)))
# }
#for k = 3, implement sigmoidal function with mean 0
k = 3
amp = d[k]
for(t in 1:tsize[k]){
  time = short_year_ki[k,t]
  d_kt[k,t] = -abs(amp) +  2*abs(amp)/(1+exp(-1/3*(time - T/2)))
}
d[k] = mean(d_kt[k,1:tsize[k]])
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

jags_file = "/jags/BNMA_Like_Bin_Trial_Two_Arm.bug"

#fit model

jags_fit <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)

print(jags_fit)

post_d = jags_fit$BUGSoutput$mean$d
plot(d,post_d, xlab = "True d_{1k}", ylab = "Posterior Mean d_{1k}", main = "BNMA")
lines(-5:5,-5:5)
#seems to give accurate inference to mean
#bad at estimating end-time effect

post_delta = jags_fit$BUGSoutput$mean$delta
plot(delta,post_delta)
lines(-5:5,-5:5)
#again, accurate inference to mean

#calculate sucra

#extract mcmc chains
d_mat = jags_fit$BUGSoutput$sims.matrix[,1:5]
#calc mc rank prob for each treatment
bnma_rank_prob = matrix(0,nrow = K, ncol = K)
for(n in 1:nrow(d_mat)){
  temp = order(d_mat[n,],decreasing = TRUE)
  for(i in 1:length(temp)){
    bnma_rank_prob[temp[i],i] = bnma_rank_prob[temp[i],i] + 1
  }
}
for(k in 1:K){
  bnma_rank_prob[k,] = bnma_rank_prob[k,]/sum(bnma_rank_prob[k,])
}
#calc sucra for each treatment
bnma_sucra = rep(0,K)
for(k in 1:K){
  for(i in 1:(K-1)){
    bnma_sucra[k] = bnma_sucra[k] + sum(bnma_rank_prob[k,1:i])
  }
}
bnma_sucra = bnma_sucra/(K-1)

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
jags_params <- c("z_k","d",
                 "min_k", "height_k",
                 "scale_k", "center_k",
                 "psi")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_Sigmoidal_v2.bug"

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

ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
  geom_line() +
  facet_grid(~K) +
  geom_ribbon(aes(ymin = Low, ymax = High),
              alpha = 0.2, fill = "deepskyblue4")


# post_delta = jags_fit_sig$BUGSoutput$mean$delta
# plot(delta,post_delta)
# lines(-5:5,-5:5)
#again, accurate inference to mean


#calculate sucra

#extract mcmc chains
d_mat = matrix(0,nrow = nrow(jags_fit_sig$BUGSoutput$sims.matrix), ncol = K)
#find each d_k^T (might be easier to do this in JAGS in future?)
for(k in 1:K){
  
  d_k = jags_fit_sig$BUGSoutput$sims.matrix[,K + k]
  a_k = jags_fit_sig$BUGSoutput$sims.matrix[,2*K + 1 + k]
  b_k = jags_fit_sig$BUGSoutput$sims.matrix[,3*K + 2 + k]
  c_k = jags_fit_sig$BUGSoutput$sims.matrix[,k]
  psi = jags_fit_sig$BUGSoutput$sims.matrix[,17]
  
  z_k = jags_fit_sig$BUGSoutput$sims.matrix[,22 + k]
  
  for(it in 1:nrow(d_mat)){
    sig_draw = rnorm(1,(d_k[it] - a_k[it]) + 2*(d_k[it] + a_k[it])/(1+exp(-b_k[it]*(T - c_k[it]))),psi^2)
    d_mat[it,k] = z_k[it]*sig_draw + (1-z_k[it])*d_k[it]
  }
  
}

#calc mc rank prob for each treatment
bnma_sig_rank_prob = matrix(0,nrow = K, ncol = K)
for(n in 1:nrow(d_mat)){
  temp = order(d_mat[n,],decreasing = TRUE)
  for(i in 1:length(temp)){
    bnma_sig_rank_prob[temp[i],i] = bnma_sig_rank_prob[temp[i],i] + 1
  }
}
for(k in 1:K){
  bnma_sig_rank_prob[k,] = bnma_sig_rank_prob[k,]/sum(bnma_sig_rank_prob[k,])
}
#calc sucra for each treatment
bnma_sig_sucra = rep(0,K)
for(k in 1:K){
  for(i in 1:(K-1)){
    bnma_sig_sucra[k] = bnma_sig_sucra[k] + sum(bnma_sig_rank_prob[k,1:i])
  }
}
bnma_sig_sucra = bnma_sig_sucra/(K-1)
# 
# ###################
# # Meta-regression #
# ###################
# 
# m_mu = 0
# prec_mu = 1
# 
# #impelment a time model
# #save data for jags
# jags_data <- list("I","n_ik","y_ik",
#                   "m_mu","prec_mu",
#                   "t2","t1","K",
#                   "years")
# 
# #note which params to save
# jags_params <- c("b","d","sd")
# 
# #define inititailization values
# jags_inits <- function(){
#   list("sd" = runif(1,0,5)
#   )
# }
# 
# jags_file = "/jags/BNMA_Like_Bin_Trial_Two_Arm_Time.bug"
# 
# #fit model
# set.seed(123)
# 
# jags_fit_meta = R2jags::jags(data = jags_data, 
#                              inits = jags_inits,
#                              parameters.to.save = jags_params,
#                              n.chains = 3, n.iter = 20000,n.burnin = 10000,
#                              model.file =  jags_file
# )
# 
# print(jags_fit_meta)
# 
# b3 = jags_fit_meta$BUGSoutput$sims.matrix[,3]
# d3 = jags_fit_meta$BUGSoutput$sims.matrix[,8]
# 
# mean(d3 + 50*b3)
# quantile(d3 + 50*b3,c(0.025,0.975))
# 
# 
# #calculate sucra
# 
# #extract mcmc chains
# d_mat = matrix(0,nrow = nrow(jags_fit$BUGSoutput$sims.matrix), ncol = K)
# #find each d_k^T (might be easier to do this in JAGS in future?)
# for(k in 1:K){
#   d_mat[,k] = jags_fit_meta$BUGSoutput$sims.matrix[,K + k] + T*jags_fit_meta$BUGSoutput$sims.matrix[,k]
# }
# 
# #calc mc rank prob for each treatment
# bnma_meta_rank_prob = matrix(0,nrow = K, ncol = K)
# for(n in 1:nrow(d_mat)){
#   temp = order(d_mat[n,],decreasing = TRUE)
#   for(i in 1:length(temp)){
#     bnma_meta_rank_prob[temp[i],i] = bnma_meta_rank_prob[temp[i],i] + 1
#   }
# }
# for(k in 1:K){
#   bnma_meta_rank_prob[k,] = bnma_meta_rank_prob[k,]/sum(bnma_meta_rank_prob[k,])
# }
# #calc sucra for each treatment
# bnma_meta_sucra = rep(0,K)
# for(k in 1:K){
#   for(i in 1:(K-1)){
#     bnma_meta_sucra[k] = bnma_meta_sucra[k] + sum(bnma_meta_rank_prob[k,1:i])
#   }
# }
# bnma_meta_sucra = bnma_meta_sucra/(K-1)
# 
# ###########
# # GP BNMA #
# ###########
# 
# m_mu = 0
# prec_mu = 1
# 
# #JAGS Now
# 
# #save data for jags
# jags_data <- list("I", "n_ik","y_ik",
#                   "m_mu","prec_mu",
#                   "t2","t1","K",
#                   "short_year_ki",
#                   "tsize","TS",
#                   "lts_ind")
# 
# #note which params to save
# jags_params <- c("d","d_kt","sd",
#                  "phi","rho","psi",
#                  "delta")
# 
# #define inititailization values
# jags_inits <- function(){
#   list("sd" = runif(1,0,5)
#   )
# }
# 
# jags_file = "/jags/BNMA_Like_Bin_Trial_Two_Arm_Time_GP.bug"
# 
# #fit model
# jags_fit_time <- R2jags::jags(data = jags_data,
#                               inits = jags_inits,
#                               parameters.to.save = jags_params,
#                               n.chains = 3, n.iter = 9000,n.burnin = 1000,
#                               model.file =  jags_file
# )
# 
# print(jags_fit_time)
# 
# post_d_time = jags_fit_time$BUGSoutput$mean$d
# plot(d,post_d_time, xlab = "True d_{1k}", ylab = "Posterior Mean d_{1k}", main = "GP-BNMA")
# lines(-5:5,-5:5)
# 
# post_delta_time = jags_fit_time$BUGSoutput$mean$delta
# plot(delta,post_delta_time)
# lines(-5:5,-5:5)
# 
# d_kt_post = jags_fit_time$BUGSoutput$mean$d_kt
# 
# d_kt_list = NULL
# years_kt = NULL
# for(k in 1:K){
#   print(mean(d_kt_post[k,1:tsize[k]]))
#   d_kt_list[[k]] = d_kt_post[k,1:tsize[k]]
#   years_kt[[k]] = short_year_ki[k,1:tsize[k]]
# }
# #very slightly different results - found telehealth more effective?
# 
# post_rho = jags_fit_time$BUGSoutput$mean$rho
# post_phi = jags_fit_time$BUGSoutput$mean$phi
# post_psi = jags_fit_time$BUGSoutput$mean$psi
# post_d = jags_fit_time$BUGSoutput$mean$d
# 
# jags_fit_time$BUGSoutput$summary[1:5,c(1,3,7)]
# 
# 
# post_sd = jags_fit_time$BUGSoutput$mean$sd
# 
# plot_years = NULL
# plot_pred_mu = NULL
# plot_pred_low = NULL
# plot_pred_high = NULL
# plot_k = NULL
# 
# for(k in 2:K){
#   obs_x = years_kt[[k]]
#   y = d_kt_list[[k]]
#   pred_x = 0:(10*T)/10
#   
#   mu1 = rep(post_d[k],length(pred_x))
#   mu2 = rep(post_d[k],length(obs_x))
#   
#   
#   #find covariance
#   tot_years = c(pred_x,obs_x)
#   Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
#   for(i in 1:nrow(Sigma)){
#     
#     if(i <= length(pred_x)){
#       Sigma[i,i] = post_phi[k]^2
#     } else {
#       Sigma[i,i] = post_psi^2 + post_phi[k]^2
#     }
#     
#     if(i < nrow(Sigma)){
#       for(j in (i+1):ncol(Sigma)){
#         Sigma[i,j] = post_phi[k]^2*exp(-post_rho[k]*(tot_years[i] - tot_years[j])^2 )
#         Sigma[j,i] = Sigma[i,j]
#       }
#     }
#     
#   }
#   
#   one_ind = 1:length(pred_x)
#   two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))
#   Sigma11 = Sigma[one_ind,one_ind]
#   Sigma12 = Sigma[one_ind,two_ind]
#   Sigma22 = Sigma[two_ind,two_ind]
#   Sigma22_inv = solve(Sigma22)
#   Sigma21 = Sigma[two_ind,one_ind]
#   
#   pred_mu = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
#   pred_Sigma = Sigma11 - Sigma12 %*% Sigma22_inv %*% Sigma21
#   
#   pred_low = pred_mu - 1.95 * sqrt(diag(pred_Sigma))
#   pred_high = pred_mu + 1.95 * sqrt(diag(pred_Sigma))
#   
#   plot_years = c(plot_years,pred_x)
#   plot_pred_mu = c(plot_pred_mu,pred_mu)
#   plot_pred_low = c(plot_pred_low,pred_low)
#   plot_pred_high = c(plot_pred_high,pred_high)
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
# ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
#   geom_line() +
#   facet_grid(~K) +
#   geom_ribbon(aes(ymin = Low, ymax = High),
#               alpha = 0.2, fill = "deepskyblue4")
# 
# #extract mcmc chains
# d_mat = matrix(0,nrow = nrow(jags_fit$BUGSoutput$sims.matrix), ncol = K)
# #find each d_k^T (might be easier to do this in JAGS in future?)
# for(k in 2:K){
#   #find
#   obs_x = years_kt[[k]]
#   y = d_kt_list[[k]]
#   pred_x =T
#   mu1 = rep(post_d[k],length(pred_x))
#   mu2 = rep(post_d[k],length(obs_x))
#   one_ind = 1:length(pred_x)
#   two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))
#   
#   tot_years = c(pred_x,obs_x)
#   
#   phi_vec = jags_fit_time$BUGSoutput$sims.list$phi[,k]
#   rho_vec = jags_fit_time$BUGSoutput$sims.list$rho[,k]
#   psi_vec = jags_fit_time$BUGSoutput$sims.list$psi
#   
#   
#   for(n in 1:nrow(d_mat)){
#     
#     #find covariance
#     Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
#     for(i in 1:nrow(Sigma)){
#       
#       if(i <= length(pred_x)){
#         Sigma[i,i] = phi_vec[n]^2
#       } else {
#         Sigma[i,i] = psi_vec[n]^2 + phi_vec[n]^2
#       }
#       
#       if(i < nrow(Sigma)){
#         for(j in (i+1):ncol(Sigma)){
#           Sigma[i,j] = phi_vec[n]^2*exp(-rho_vec[n]*(tot_years[i] - tot_years[j])^2 )
#           Sigma[j,i] = Sigma[i,j]
#         }
#       }
#       
#     }
#     
#     Sigma11 = Sigma[one_ind,one_ind]
#     Sigma12 = Sigma[one_ind,two_ind]
#     Sigma22 = Sigma[two_ind,two_ind]
#     Sigma22_inv = solve(Sigma22)
#     Sigma21 = Sigma[two_ind,one_ind]
#     
#     d_mat[n,k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
#     
#     if(n %% 50 ==0){
#       print(n)
#     }
#   }
# }
# #calc mc rank prob for each treatment
# gpbnma_rank_prob = matrix(0,nrow = K, ncol = K)
# for(n in 1:nrow(d_mat)){
#   temp = order(d_mat[n,],decreasing = TRUE)
#   for(i in 1:length(temp)){
#     gpbnma_rank_prob[temp[i],i] = gpbnma_rank_prob[temp[i],i] + 1
#   }
# }
# for(k in 1:K){
#   gpbnma_rank_prob[k,] = gpbnma_rank_prob[k,]/sum(gpbnma_rank_prob[k,])
# }
# #calc sucra for each treatment
# gpbnma_sucra = rep(0,K)
# for(k in 1:K){
#   for(i in 1:(K-1)){
#     gpbnma_sucra[k] = gpbnma_sucra[k] + sum(gpbnma_rank_prob[k,1:i])
#   }
# }
# bnma_gp_sucra = gpbnma_sucra/(K-1)
