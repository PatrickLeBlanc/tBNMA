library(ggplot2)
library(R2jags)
library(gridExtra)

find_quant = function(x){
  return(quantile(x,c(0.025,0.975)))
}

#Make data structures
plot_years = NULL
plot_pred_mu = NULL
plot_pred_low = NULL
plot_pred_high = NULL
plot_mod = NULL
plot_truth = NULL

dat = read.csv("Data/data_date_tbnma.csv")

########
# make data structures

#find number of studies/treatments
uniq_studies = unique(dat$Study)
I = length(uniq_studies)
K = length(table(dat$Treatment))

#find dates in which the studies took place
#normalize so that the first year is 0
days = NULL
years = NULL
for(i in 1:length(uniq_studies)){
  temp_dat = dat[dat$Study == uniq_studies[i],]
  temp = temp_dat$Date[1]
  days = c(days,strftime(as.POSIXct(temp), format = "%j"))
  years = c(years,strftime(as.POSIXct(temp), format = "%Y"))
}
days = as.numeric(days)
years = as.numeric(years)
dates = years + days/365 #normalize assuming 365 days/year - close enough?

#reconfigure names to match up with old syntax
years = dates

first_year = min(years)
years = years - first_year
T = max(years)

#list studies and treatments
studies = 1:I
treatments = 1:K

#make map between study number and study
study_map = data.frame("Number" = studies,
                       "Study" = uniq_studies)

#make map between treatment number and treatmetns
treatment_map = data.frame("Number" = treatments,
                           "Treatment" = unique(dat$Treatment))

#make lin 1 and van 2
treatment_map$Treatment[1] = "LIN"
treatment_map$Treatment[2] = "VAN"

#find num_treat
num_treat = rep(0,I)
for(i in 1:I){
  temp_dat = dat[which(dat$Study == uniq_studies[i]),]
  num_treat[i] = nrow(temp_dat)
}

#find treatment matrix
#t_mat is three columns, one for each treamtent arm (blanks are 0)
#indicator for which treatment is in whcih arm
t_mat = matrix(0,nrow = I, ncol = max(num_treat))
for(i in 1:I){
  
  temp_dat = dat[which(dat$Study == study_map$Study[i]),]
  
  for(j in 1:num_treat[i]){
    t_mat[i,j] = which(treatment_map$Treatment == temp_dat$Treatment[j])
  }
  #want lower numbered treatments first
  t_mat[i,1:num_treat[i]] = sort(t_mat[i,1:num_treat[i]])
}

#############
# First Sim #
#############

set.seed(2)

#simulate sd from unif(0,1)
sd = runif(1,0,0.5)
tau = 1/sd^2

#simulate d according to standard normal
d = rep(0,K)
for(k in 2:K){
  d[k] = rnorm(1,0,1)
}
# d = post_d

#simulate d_kt?

#find how many datapoints for each treatment
tsize = rep(0,K)
for(k in 1:K){
  for(i in 1:I){
    if(sum(t_mat[i,] ==k) > 0 ){
      tsize[k] = tsize[k] + 1
    } 
  }
}

TS = max(tsize)


#which years each study appears, indexed by which study they appear in
#treatments are rows, columns are studies
long_year_ki = matrix(NA,nrow = K, ncol = I)
for(k in 1:K){
  for(i in 1:I){
    if(sum(t_mat[i,] ==k) > 0 ){
      long_year_ki[k,i] = years[i]
    } 
  }
}

#short_year_ki[k,1:tsize[k]] returns a list indicating the indices i
#of all studies where treatment k is from
short_year_ki = matrix(NA,nrow = K, ncol = max(tsize))
for(k in 1:K){
  temp = long_year_ki[k,]
  temp = temp[!is.na(temp)]
  for(t in 1:tsize[k]){
    short_year_ki[k,t] = temp[t]
  }
}

#years_kt[[k]] contains the years in which treatment k appears
#years_kt[[k]] = short_year_ki[k,]
years_kt = NULL
for(k in 1:K){
  years_kt[[k]] = short_year_ki[k,1:tsize[k]]
}

#need a map leading from study to short index for each treatment
#lts_ind[k,i], for k in 1:K and i in 1:I, returns an index
# which will map to the short_year_ki format
#that is, we take an index on the space of all studies, and map to an index
# on the space only of studies involving treatment k
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

#####################
# Flat time effects #
#####################

d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    d_kt[k,t] = d[k]
  }
}


#simulate delta for each study
#implicitly assume that delta for t1 is zero and that there are only 2 treatments
delta = matrix(0,nrow = I, ncol = max(num_treat))
w = matrix(0,nrow = nrow(delta),ncol = ncol(delta))
sw = matrix(0,nrow = nrow(delta),ncol = ncol(delta))
for(i in 1:I){
  #first arm
  delta[i,1] = 0
  w[i,1] = 0
  
  
  #other arms
  for(j in 2:num_treat[i]){
    #do multi-arm adjusemtns
    taud = tau*2*(j-1)/j
    sdd = taud^(-1/2)
    w[i,j] = delta[i,j] - d_kt[t_mat[i,j],lts_ind[t_mat[i,j],i]] + d_kt[t_mat[i,1],lts_ind[t_mat[i,1],i]]
    sw[i,j] = sum(w[i,1:(j-1)])/(j-1)
    
    #find means and sim
    md = d_kt[t_mat[i,j],lts_ind[t_mat[i,j],i]] - d_kt[t_mat[i,1],lts_ind[t_mat[i,1],i]] + sw[i,j]
    delta[i,j] = rnorm(1,md,sdd)
    # delta[i,j] = md
  }
  
}

#simulate trial specific baseline effects
mu_i = rnorm(I,0,sd = 1)
# mu_i = rep(0,I)
# mu_i[12] = abs(mu_i[12])

#simulate number of pateitns
n_ik = matrix(0,nrow = I, ncol = max(num_treat))
for(i in 1:I){
  study = study_map$Study[i]
  temp_dat = dat[dat$Study == study,]
  for(j in 1:num_treat[i]){
    #find the treatment at hand
    treat = temp_dat$Treatment[j]
    
    #find numeric label for treat
    treat_label = which(treatment_map$Treatment == treat)
    
    #find out which treatment arm this treatment should belong to
    treat_ind = which(t_mat[i,] == treat_label)
    
    #update trials/successes
    n_ik[i,treat_ind] = temp_dat$Trials[j]
  }
}


#simulate y
y_ik = matrix(0,nrow = I, ncol = max(num_treat))
for(i in 1:I){
  
  for(k in 1:num_treat[i]){
    p = exp(mu_i[i] + delta[i,k])/(1+exp(mu_i[i] + delta[i,k]))
    y_ik[i,k] = rbinom(1,n_ik[i,k],p)
  }
  
}

########
# BNMA #
########

#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "num_treat", "t_mat",
                  "K")

#note which params to save
jags_params <- c("d","sd_d", "m_d",
                 "sd",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}


dir = getwd()
jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm.bug")

#fit model

jags_fit <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)

#Extract posterior credible/posterior predictive distributions

k = 2
obs_x = years_kt[[k]]
pred_x = 0:(10*T)/10

#find credible intervals for d_kt
#just replicate d across all time points, lol
d_vec = jags_fit$BUGSoutput$sims.list$d[,k]
pred_mu = rep(mean(d_vec),length(pred_x))
quant_mat = find_quant(d_vec)
pred_low = rep(quant_mat[1],length(pred_x))
pred_high = rep(quant_mat[2],length(pred_x))

#put into dataframe
plot_years = c(plot_years,first_year + pred_x)
plot_pred_mu = c(plot_pred_mu,pred_mu)
plot_pred_low = c(plot_pred_low,pred_low)
plot_pred_high = c(plot_pred_high,pred_high)
plot_mod = c(plot_mod,rep("BNMA",length(pred_x)))
plot_truth = c(plot_truth,rep("Constant Effect",length(pred_x)))

#get true values
plot_years = c(plot_years,first_year + pred_x)
plot_pred_mu = c(plot_pred_mu,rep(d[k],length(pred_x)))
plot_pred_low = c(plot_pred_low,rep(d[k],length(pred_x)))
plot_pred_high = c(plot_pred_high,rep(d[k],length(pred_x)))
plot_mod = c(plot_mod,rep("Truth",length(pred_x)))
plot_truth = c(plot_truth,rep("Constant Effect",length(pred_x)))

###################
# Meta-regression #
###################

#impelment a time model
#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "num_treat","t_mat",
                  "K",
                  "years")

#note which params to save
jags_params <- c("b","d","sd","z_k")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

dir = getwd()
jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm_Time_Z.bug")

jags_fit_meta = R2jags::jags(data = jags_data,
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)


k = 2
obs_x = years_kt[[k]]
pred_x = 0:(10*T)/10

#find credible intervals for d_kt
#just replicate d across all time points, lol
d_vec = jags_fit_meta$BUGSoutput$sims.list$d[,k]
b_vec = jags_fit_meta$BUGSoutput$sims.list$b[,k]

#find posterior credible intervals for the d_kt
for(t in pred_x){
  pred_vec = d_vec + b_vec*t
  
  pred_mu = mean(pred_vec)
  
  quant_mat = find_quant(pred_vec)
  pred_low = quant_mat[1]
  pred_high = quant_mat[2]
  
  plot_years = c(plot_years,first_year + t)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_mod = c(plot_mod,"Meta-BNMA")
  plot_truth = c(plot_truth,"Constant Effect")
}

###########
# GP BNMA #
###########
#JAGS Now

#save data for jags
jags_data <- list("I", "n_ik","y_ik",
                  "num_treat","t_mat",
                  "K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind")

#note which params to save
jags_params <- c("d","d_kt","sd",
                 "phi","rho","psi",
                 "sb","sl",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

dir = getwd()
jags_file = paste0(dir,"/jags/tBNMA_GP_LIN_AR1.bug")

#fit model
jags_fit_time <- R2jags::jags(data = jags_data,
                              inits = jags_inits,
                              parameters.to.save = jags_params,
                              n.chains = 3, n.iter = 9000,n.burnin = 1000,
                              model.file =  jags_file
)


# for(k in c(2,8,9,10,15)){
for(k in c(2)){

  # obs_x = years_kt[[k]]
  # pred_x = obs_x
  # d_kt_mat = jags_fit_time$BUGSoutput$sims.list$d_kt[,k,]
  #
  # pred_mu = apply(d_kt_mat,2,mean)
  # quant_mat = apply(d_kt_mat,2,find_quant)
  # pred_low = quant_mat[1,]
  # pred_high = quant_mat[2,]

  #find
  obs_x = years_kt[[k]]
  pred_x = 0:(10*T)/10
  # pred_x = 0:(10*(T+1))/10
  # mu1 = rep(post_d[k],length(pred_x))
  # mu2 = rep(post_d[k],length(obs_x))
  one_ind = 1:length(pred_x)
  two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))

  tot_years = c(pred_x,obs_x)

  phi_vec = jags_fit_time$BUGSoutput$sims.list$phi[,k]
  rho_vec = jags_fit_time$BUGSoutput$sims.list$rho[,k]
  psi_vec = jags_fit_time$BUGSoutput$sims.list$psi
  sb_vec = jags_fit_time$BUGSoutput$sims.list$sb[,k]
  sl_vec = jags_fit_time$BUGSoutput$sims.list$sl[,k]
  d_vec =  jags_fit_time$BUGSoutput$sims.list$d[,k]
  d_kt_mat = jags_fit_time$BUGSoutput$sims.list$d_kt[,k,]


  out_ntk = array(0,dim = c(length(psi_vec),length(pred_x),K))
  for(n in 1:nrow(d_kt_mat)){
    #find y
    y = d_kt_mat[n,1:tsize[k]]

    #find mean
    mu1 = rep(d_vec[n],length(pred_x))
    mu2 = rep(d_vec[n],length(obs_x))

    #find covariance
    Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
    for(i in 1:nrow(Sigma)){

      if(i <= length(pred_x)){
        Sigma[i,i] = phi_vec[n]^2 + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[i]
      } else {
        Sigma[i,i] = psi_vec[n]^2 + phi_vec[n]^2 + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[i]
      }

      if(i < nrow(Sigma)){
        for(j in (i+1):ncol(Sigma)){
          Sigma[i,j] = phi_vec[n]^2*exp(-rho_vec[n]*abs(tot_years[i] - tot_years[j])) + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[j]
          Sigma[j,i] = Sigma[i,j]
        }
      }

    }

    Sigma11 = Sigma[one_ind,one_ind]
    Sigma12 = Sigma[one_ind,two_ind]
    Sigma22 = Sigma[two_ind,two_ind]
    Sigma22_inv = solve(Sigma22)
    Sigma21 = Sigma[two_ind,one_ind]

    # out_ntk[n,,k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)

    cond_mu = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
    cond_Sigma = Sigma11 - Sigma12 %*% Sigma22_inv %*% Sigma21

    # A = chol(cond_Sigma)
    A = svd(cond_Sigma)
    A1 = A$u %*% diag(sqrt(A$d))

    out_ntk[n,,k] = (A1) %*% rnorm(length(pred_x)) + cond_mu

    if(n %% 500 ==0){
      print(n)
    }
  }

  pred_mu = apply(out_ntk[,,k],2,mean)
  quant_mat = apply(out_ntk[,,k],2,find_quant)
  pred_low = quant_mat[1,]
  pred_high = quant_mat[2,]

  plot_years = c(plot_years,first_year + pred_x)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_mod = c(plot_mod,rep("tBNMA",length(pred_x)))
  plot_truth = c(plot_truth,rep("Constant Effect",length(pred_x)))

}


##############
# Second Sim #
##############

set.seed(2)

#simulate sd from unif(0,1)
sd = runif(1,0,0.5)
tau = 1/sd^2

#simulate d according to standard normal
d = rep(0,K)
for(k in 2:K){
  d[k] = rnorm(1,0,1)
}
# d = post_d

#simulate d_kt?

#find how many datapoints for each treatment
tsize = rep(0,K)
for(k in 1:K){
  for(i in 1:I){
    if(sum(t_mat[i,] ==k) > 0 ){
      tsize[k] = tsize[k] + 1
    } 
  }
}

TS = max(tsize)


#which years each study appears, indexed by which study they appear in
#treatments are rows, columns are studies
long_year_ki = matrix(NA,nrow = K, ncol = I)
for(k in 1:K){
  for(i in 1:I){
    if(sum(t_mat[i,] ==k) > 0 ){
      long_year_ki[k,i] = years[i]
    } 
  }
}

#short_year_ki[k,1:tsize[k]] returns a list indicating the indices i
#of all studies where treatment k is from
short_year_ki = matrix(NA,nrow = K, ncol = max(tsize))
for(k in 1:K){
  temp = long_year_ki[k,]
  temp = temp[!is.na(temp)]
  for(t in 1:tsize[k]){
    short_year_ki[k,t] = temp[t]
  }
}

#need a map leading from study to short index for each treatment
#lts_ind[k,i], for k in 1:K and i in 1:I, returns an index
# which will map to the short_year_ki format
#that is, we take an index on the space of all studies, and map to an index
# on the space only of studies involving treatment k
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

##########################
# Quadratic time effects #
##########################

d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    d_kt[k,t] = d[k]
  }
}

k = 2
# height = rgamma(1,1,1)
# scale = rgamma(1,1,1)
# center = runif(1,0,T)
beta0 = -1
beta1 = -0.1
beta2 = 0.04
for(t in 1:tsize[k]){
  d_kt[k,t] = d[k] + +beta0 + beta1*(short_year_ki[k,t] - T/2) + beta2*(short_year_ki[k,t] - T/2)^2
}
# plot(short_year_ki[k,1:tsize[k]],d_kt[k,1:tsize[k]])


#simulate delta for each study
#implicitly assume that delta for t1 is zero and that there are only 2 treatments
delta = matrix(0,nrow = I, ncol = max(num_treat))
w = matrix(0,nrow = nrow(delta),ncol = ncol(delta))
sw = matrix(0,nrow = nrow(delta),ncol = ncol(delta))
for(i in 1:I){
  #first arm
  delta[i,1] = 0
  w[i,1] = 0
  
  
  #other arms
  for(j in 2:num_treat[i]){
    #do multi-arm adjusemtns
    taud = tau*2*(j-1)/j
    sdd = taud^(-1/2)
    w[i,j] = delta[i,j] - d_kt[t_mat[i,j],lts_ind[t_mat[i,j],i]] + d_kt[t_mat[i,1],lts_ind[t_mat[i,1],i]]
    sw[i,j] = sum(w[i,1:(j-1)])/(j-1)
    
    #find means and sim
    md = d_kt[t_mat[i,j],lts_ind[t_mat[i,j],i]] - d_kt[t_mat[i,1],lts_ind[t_mat[i,1],i]] + sw[i,j]
    delta[i,j] = rnorm(1,md,sdd)
    # delta[i,j] = md
  }
  
}

#simulate trial specific baseline effects
mu_i = rnorm(I,0,sd = 1)
# mu_i = rep(0,I)
# mu_i[12] = abs(mu_i[12])

#simulate number of pateitns
n_ik = matrix(0,nrow = I, ncol = max(num_treat))
for(i in 1:I){
  study = study_map$Study[i]
  temp_dat = dat[dat$Study == study,]
  for(j in 1:num_treat[i]){
    #find the treatment at hand
    treat = temp_dat$Treatment[j]
    
    #find numeric label for treat
    treat_label = which(treatment_map$Treatment == treat)
    
    #find out which treatment arm this treatment should belong to
    treat_ind = which(t_mat[i,] == treat_label)
    
    #update trials/successes
    n_ik[i,treat_ind] = temp_dat$Trials[j]
  }
}


#simulate y
y_ik = matrix(0,nrow = I, ncol = max(num_treat))
for(i in 1:I){
  
  for(k in 1:num_treat[i]){
    p = exp(mu_i[i] + delta[i,k])/(1+exp(mu_i[i] + delta[i,k]))
    y_ik[i,k] = rbinom(1,n_ik[i,k],p)
  }
  
}

########
# BNMA #
########

#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "num_treat", "t_mat",
                  "K")

#note which params to save
jags_params <- c("d","sd_d", "m_d",
                 "sd",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}


dir = getwd()
jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm.bug")

#fit model

jags_fit <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)

#Extract posterior credible/posterior predictive distributions

#Make data structures
k = 2
obs_x = years_kt[[k]]
pred_x = 0:(10*T)/10

#find credible intervals for d_kt
#just replicate d across all time points, lol
d_vec = jags_fit$BUGSoutput$sims.list$d[,k]
pred_mu = rep(mean(d_vec),length(pred_x))
quant_mat = find_quant(d_vec)
pred_low = rep(quant_mat[1],length(pred_x))
pred_high = rep(quant_mat[2],length(pred_x))

#put into dataframe
plot_years = c(plot_years,first_year + pred_x)
plot_pred_mu = c(plot_pred_mu,pred_mu)
plot_pred_low = c(plot_pred_low,pred_low)
plot_pred_high = c(plot_pred_high,pred_high)
plot_mod = c(plot_mod,rep("BNMA",length(pred_x)))
plot_truth = c(plot_truth,rep("Quadratic Effect",length(pred_x)))

#get true values 
beta0 = -1
beta1 = -0.1
beta2 = 0.04
for(t in 1:tsize[k]){
  d_kt[k,t] = d[k] + +beta0 + beta1*(short_year_ki[k,t] - T/2) + beta2*(short_year_ki[k,t] - T/2)^2
}

pred_mu = beta0 + beta1*(pred_x - T/2) + beta2*(pred_x-T/2)^2
pred_low = pred_mu
pred_high = pred_mu

#put true values into dataframe
plot_years = c(plot_years,first_year + pred_x)
plot_pred_mu = c(plot_pred_mu,pred_mu)
plot_pred_low = c(plot_pred_low,pred_low)
plot_pred_high = c(plot_pred_high,pred_high)
plot_mod = c(plot_mod,rep("Truth",length(pred_x)))
plot_truth = c(plot_truth,rep("Quadratic Effect",length(pred_x)))

###################
# Meta-regression #
###################

#impelment a time model
#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "num_treat","t_mat",
                  "K",
                  "years")

#note which params to save
jags_params <- c("b","d","sd","z_k")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

dir = getwd()
jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm_Time_Z.bug")

jags_fit_meta = R2jags::jags(data = jags_data,
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)


k = 2
obs_x = years_kt[[k]]
pred_x = 0:(10*T)/10

#find credible intervals for d_kt
#just replicate d across all time points, lol
d_vec = jags_fit_meta$BUGSoutput$sims.list$d[,k]
b_vec = jags_fit_meta$BUGSoutput$sims.list$b[,k]

#find posterior credible intervals for the d_kt
for(t in pred_x){
  pred_vec = d_vec + b_vec*t
  
  pred_mu = mean(pred_vec)
  
  quant_mat = find_quant(pred_vec)
  pred_low = quant_mat[1]
  pred_high = quant_mat[2]
  
  plot_years = c(plot_years,first_year + t)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_mod = c(plot_mod,"Meta-BNMA")
  plot_truth = c(plot_truth,"Quadratic Effect")
}

###########
# GP BNMA #
###########
#JAGS Now

#save data for jags
jags_data <- list("I", "n_ik","y_ik",
                  "num_treat","t_mat",
                  "K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind")

#note which params to save
jags_params <- c("d","d_kt","sd",
                 "phi","rho","psi",
                 "sb","sl",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

dir = getwd()
jags_file = paste0(dir,"/jags/tBNMA_GP_LIN_AR1.bug")

#fit model
jags_fit_time <- R2jags::jags(data = jags_data,
                              inits = jags_inits,
                              parameters.to.save = jags_params,
                              n.chains = 3, n.iter = 9000,n.burnin = 1000,
                              model.file =  jags_file
)


# for(k in c(2,8,9,10,15)){
for(k in c(2)){
  
  # obs_x = years_kt[[k]]
  # pred_x = obs_x
  # d_kt_mat = jags_fit_time$BUGSoutput$sims.list$d_kt[,k,]
  #
  # pred_mu = apply(d_kt_mat,2,mean)
  # quant_mat = apply(d_kt_mat,2,find_quant)
  # pred_low = quant_mat[1,]
  # pred_high = quant_mat[2,]
  
  #find
  obs_x = years_kt[[k]]
  pred_x = 0:(10*T)/10
  # pred_x = 0:(10*(T+1))/10
  # mu1 = rep(post_d[k],length(pred_x))
  # mu2 = rep(post_d[k],length(obs_x))
  one_ind = 1:length(pred_x)
  two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))
  
  tot_years = c(pred_x,obs_x)
  
  phi_vec = jags_fit_time$BUGSoutput$sims.list$phi[,k]
  rho_vec = jags_fit_time$BUGSoutput$sims.list$rho[,k]
  psi_vec = jags_fit_time$BUGSoutput$sims.list$psi
  sb_vec = jags_fit_time$BUGSoutput$sims.list$sb[,k]
  sl_vec = jags_fit_time$BUGSoutput$sims.list$sl[,k]
  d_vec =  jags_fit_time$BUGSoutput$sims.list$d[,k]
  d_kt_mat = jags_fit_time$BUGSoutput$sims.list$d_kt[,k,]
  
  
  out_ntk = array(0,dim = c(length(psi_vec),length(pred_x),K))
  for(n in 1:nrow(d_kt_mat)){
    #find y
    y = d_kt_mat[n,1:tsize[k]]
    
    #find mean
    mu1 = rep(d_vec[n],length(pred_x))
    mu2 = rep(d_vec[n],length(obs_x))
    
    #find covariance
    Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
    for(i in 1:nrow(Sigma)){
      
      if(i <= length(pred_x)){
        Sigma[i,i] = phi_vec[n]^2 + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[i]
      } else {
        Sigma[i,i] = psi_vec[n]^2 + phi_vec[n]^2 + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[i]
      }
      
      if(i < nrow(Sigma)){
        for(j in (i+1):ncol(Sigma)){
          Sigma[i,j] = phi_vec[n]^2*exp(-rho_vec[n]*abs(tot_years[i] - tot_years[j])) + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[j]
          Sigma[j,i] = Sigma[i,j]
        }
      }
      
    }
    
    Sigma11 = Sigma[one_ind,one_ind]
    Sigma12 = Sigma[one_ind,two_ind]
    Sigma22 = Sigma[two_ind,two_ind]
    Sigma22_inv = solve(Sigma22)
    Sigma21 = Sigma[two_ind,one_ind]
    
    # out_ntk[n,,k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
    
    cond_mu = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
    cond_Sigma = Sigma11 - Sigma12 %*% Sigma22_inv %*% Sigma21
    
    # A = chol(cond_Sigma)
    A = svd(cond_Sigma)
    A1 = A$u %*% diag(sqrt(A$d))
    
    out_ntk[n,,k] = (A1) %*% rnorm(length(pred_x)) + cond_mu
    
    if(n %% 500 ==0){
      print(n)
    }
  }
  
  pred_mu = apply(out_ntk[,,k],2,mean)
  quant_mat = apply(out_ntk[,,k],2,find_quant)
  pred_low = quant_mat[1,]
  pred_high = quant_mat[2,]
  
  plot_years = c(plot_years,first_year + pred_x)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_mod = c(plot_mod,rep("tBNMA",length(pred_x)))
  plot_truth = c(plot_truth,rep("Quadratic Effect",length(pred_x)))
  
}


#############
# Third Sim #
#############

set.seed(2)

#simulate sd from unif(0,1)
sd = runif(1,0,0.5)
tau = 1/sd^2

#simulate d according to standard normal
d = rep(0,K)
for(k in 2:K){
  d[k] = rnorm(1,0,1)
}
# d = post_d

#simulate d_kt?

#find how many datapoints for each treatment
tsize = rep(0,K)
for(k in 1:K){
  for(i in 1:I){
    if(sum(t_mat[i,] ==k) > 0 ){
      tsize[k] = tsize[k] + 1
    } 
  }
}

TS = max(tsize)


#which years each study appears, indexed by which study they appear in
#treatments are rows, columns are studies
long_year_ki = matrix(NA,nrow = K, ncol = I)
for(k in 1:K){
  for(i in 1:I){
    if(sum(t_mat[i,] ==k) > 0 ){
      long_year_ki[k,i] = years[i]
    } 
  }
}

#short_year_ki[k,1:tsize[k]] returns a list indicating the indices i
#of all studies where treatment k is from
short_year_ki = matrix(NA,nrow = K, ncol = max(tsize))
for(k in 1:K){
  temp = long_year_ki[k,]
  temp = temp[!is.na(temp)]
  for(t in 1:tsize[k]){
    short_year_ki[k,t] = temp[t]
  }
}

#need a map leading from study to short index for each treatment
#lts_ind[k,i], for k in 1:K and i in 1:I, returns an index
# which will map to the short_year_ki format
#that is, we take an index on the space of all studies, and map to an index
# on the space only of studies involving treatment k
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

##########################
# Sigmoidal time effects #
##########################

d_kt = matrix(0,nrow = K, ncol = TS)
for(k in 2:K){
  #set everything equal to d
  for(t in 1:tsize[k]){
    d_kt[k,t] = d[k]
  }
}

k = 2
# height = rgamma(1,1,1)
# scale = rgamma(1,1,1)
# center = runif(1,0,T)
height = 3
scale = 1
center = T/2
for(t in 1:tsize[k]){
  d_kt[k,t] = d[k] - height +  2*height/(1+exp(-scale*(short_year_ki[k,t] - center)))
}
# plot(short_year_ki[k,1:tsize[k]],d_kt[k,1:tsize[k]])


#simulate delta for each study
#implicitly assume that delta for t1 is zero and that there are only 2 treatments
delta = matrix(0,nrow = I, ncol = max(num_treat))
w = matrix(0,nrow = nrow(delta),ncol = ncol(delta))
sw = matrix(0,nrow = nrow(delta),ncol = ncol(delta))
for(i in 1:I){
  #first arm
  delta[i,1] = 0
  w[i,1] = 0
  
  
  #other arms
  for(j in 2:num_treat[i]){
    #do multi-arm adjusemtns
    taud = tau*2*(j-1)/j
    sdd = taud^(-1/2)
    w[i,j] = delta[i,j] - d_kt[t_mat[i,j],lts_ind[t_mat[i,j],i]] + d_kt[t_mat[i,1],lts_ind[t_mat[i,1],i]]
    sw[i,j] = sum(w[i,1:(j-1)])/(j-1)
    
    #find means and sim
    md = d_kt[t_mat[i,j],lts_ind[t_mat[i,j],i]] - d_kt[t_mat[i,1],lts_ind[t_mat[i,1],i]] + sw[i,j]
    delta[i,j] = rnorm(1,md,sdd)
    # delta[i,j] = md
  }
  
}

#simulate trial specific baseline effects
mu_i = rnorm(I,0,sd = 1)
# mu_i = rep(0,I)
# mu_i[12] = abs(mu_i[12])

#simulate number of pateitns
n_ik = matrix(0,nrow = I, ncol = max(num_treat))
for(i in 1:I){
  study = study_map$Study[i]
  temp_dat = dat[dat$Study == study,]
  for(j in 1:num_treat[i]){
    #find the treatment at hand
    treat = temp_dat$Treatment[j]
    
    #find numeric label for treat
    treat_label = which(treatment_map$Treatment == treat)
    
    #find out which treatment arm this treatment should belong to
    treat_ind = which(t_mat[i,] == treat_label)
    
    #update trials/successes
    n_ik[i,treat_ind] = temp_dat$Trials[j]
  }
}


#simulate y
y_ik = matrix(0,nrow = I, ncol = max(num_treat))
for(i in 1:I){
  
  for(k in 1:num_treat[i]){
    p = exp(mu_i[i] + delta[i,k])/(1+exp(mu_i[i] + delta[i,k]))
    y_ik[i,k] = rbinom(1,n_ik[i,k],p)
  }
  
}

########
# BNMA #
########

#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "num_treat", "t_mat",
                  "K")

#note which params to save
jags_params <- c("d","sd_d", "m_d",
                 "sd",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}


dir = getwd()
jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm.bug")

#fit model

jags_fit <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)

#Extract posterior credible/posterior predictive distributions

#Make data structures
k = 2
obs_x = years_kt[[k]]
pred_x = 0:(10*T)/10

#find credible intervals for d_kt
#just replicate d across all time points, lol
d_vec = jags_fit$BUGSoutput$sims.list$d[,k]
pred_mu = rep(mean(d_vec),length(pred_x))
quant_mat = find_quant(d_vec)
pred_low = rep(quant_mat[1],length(pred_x))
pred_high = rep(quant_mat[2],length(pred_x))

#put into dataframe
plot_years = c(plot_years,first_year + pred_x)
plot_pred_mu = c(plot_pred_mu,pred_mu)
plot_pred_low = c(plot_pred_low,pred_low)
plot_pred_high = c(plot_pred_high,pred_high)
plot_mod = c(plot_mod,rep("BNMA",length(pred_x)))
plot_truth = c(plot_truth,rep("Sigmoidal Effect",length(pred_x)))

#get true values
pred_mu = d[k] - height + 2*height/(1+exp(-scale*(pred_x - center)))
pred_low = pred_mu
pred_high = pred_mu

#put true values into dataframe
plot_years = c(plot_years,first_year + pred_x)
plot_pred_mu = c(plot_pred_mu,pred_mu)
plot_pred_low = c(plot_pred_low,pred_low)
plot_pred_high = c(plot_pred_high,pred_high)
plot_mod = c(plot_mod,rep("Truth",length(pred_x)))
plot_truth = c(plot_truth,rep("Sigmoidal Effect",length(pred_x)))

###################
# Meta-regression #
###################

#impelment a time model
#save data for jags
jags_data <- list("I","n_ik","y_ik",
                  "num_treat","t_mat",
                  "K",
                  "years")

#note which params to save
jags_params <- c("b","d","sd","z_k")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

dir = getwd()
jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm_Time_Z.bug")

jags_fit_meta = R2jags::jags(data = jags_data,
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 9000,n.burnin = 1000,
                             model.file =  jags_file
)


k = 2
obs_x = years_kt[[k]]
pred_x = 0:(10*T)/10

#find credible intervals for d_kt
#just replicate d across all time points, lol
d_vec = jags_fit_meta$BUGSoutput$sims.list$d[,k]
b_vec = jags_fit_meta$BUGSoutput$sims.list$b[,k]

#find posterior credible intervals for the d_kt
for(t in pred_x){
  pred_vec = d_vec + b_vec*t
  
  pred_mu = mean(pred_vec)
  
  quant_mat = find_quant(pred_vec)
  pred_low = quant_mat[1]
  pred_high = quant_mat[2]
  
  plot_years = c(plot_years,first_year + t)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_mod = c(plot_mod,"Meta-BNMA")
  plot_truth = c(plot_truth,"Sigmoidal Effect")
}

###########
# GP BNMA #
###########
#JAGS Now

#save data for jags
jags_data <- list("I", "n_ik","y_ik",
                  "num_treat","t_mat",
                  "K",
                  "short_year_ki",
                  "tsize","TS",
                  "lts_ind")

#note which params to save
jags_params <- c("d","d_kt","sd",
                 "phi","rho","psi",
                 "sb","sl",
                 "delta")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

dir = getwd()
jags_file = paste0(dir,"/jags/tBNMA_GP_LIN_AR1.bug")

#fit model
jags_fit_time <- R2jags::jags(data = jags_data,
                              inits = jags_inits,
                              parameters.to.save = jags_params,
                              n.chains = 3, n.iter = 9000,n.burnin = 1000,
                              model.file =  jags_file
)


# for(k in c(2,8,9,10,15)){
for(k in c(2)){
  
  # obs_x = years_kt[[k]]
  # pred_x = obs_x
  # d_kt_mat = jags_fit_time$BUGSoutput$sims.list$d_kt[,k,]
  #
  # pred_mu = apply(d_kt_mat,2,mean)
  # quant_mat = apply(d_kt_mat,2,find_quant)
  # pred_low = quant_mat[1,]
  # pred_high = quant_mat[2,]
  
  #find
  obs_x = years_kt[[k]]
  pred_x = 0:(10*T)/10
  # pred_x = 0:(10*(T+1))/10
  # mu1 = rep(post_d[k],length(pred_x))
  # mu2 = rep(post_d[k],length(obs_x))
  one_ind = 1:length(pred_x)
  two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))
  
  tot_years = c(pred_x,obs_x)
  
  phi_vec = jags_fit_time$BUGSoutput$sims.list$phi[,k]
  rho_vec = jags_fit_time$BUGSoutput$sims.list$rho[,k]
  psi_vec = jags_fit_time$BUGSoutput$sims.list$psi
  sb_vec = jags_fit_time$BUGSoutput$sims.list$sb[,k]
  sl_vec = jags_fit_time$BUGSoutput$sims.list$sl[,k]
  d_vec =  jags_fit_time$BUGSoutput$sims.list$d[,k]
  d_kt_mat = jags_fit_time$BUGSoutput$sims.list$d_kt[,k,]
  
  
  out_ntk = array(0,dim = c(length(psi_vec),length(pred_x),K))
  for(n in 1:nrow(d_kt_mat)){
    #find y
    y = d_kt_mat[n,1:tsize[k]]
    
    #find mean
    mu1 = rep(d_vec[n],length(pred_x))
    mu2 = rep(d_vec[n],length(obs_x))
    
    #find covariance
    Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
    for(i in 1:nrow(Sigma)){
      
      if(i <= length(pred_x)){
        Sigma[i,i] = phi_vec[n]^2 + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[i]
      } else {
        Sigma[i,i] = psi_vec[n]^2 + phi_vec[n]^2 + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[i]
      }
      
      if(i < nrow(Sigma)){
        for(j in (i+1):ncol(Sigma)){
          Sigma[i,j] = phi_vec[n]^2*exp(-rho_vec[n]*abs(tot_years[i] - tot_years[j])) + sb_vec[n]^2 + sl_vec[n]^2*tot_years[i]*tot_years[j]
          Sigma[j,i] = Sigma[i,j]
        }
      }
      
    }
    
    Sigma11 = Sigma[one_ind,one_ind]
    Sigma12 = Sigma[one_ind,two_ind]
    Sigma22 = Sigma[two_ind,two_ind]
    Sigma22_inv = solve(Sigma22)
    Sigma21 = Sigma[two_ind,one_ind]
    
    # out_ntk[n,,k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
    
    cond_mu = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
    cond_Sigma = Sigma11 - Sigma12 %*% Sigma22_inv %*% Sigma21
    
    # A = chol(cond_Sigma)
    A = svd(cond_Sigma)
    A1 = A$u %*% diag(sqrt(A$d))
    
    out_ntk[n,,k] = (A1) %*% rnorm(length(pred_x)) + cond_mu
    
    if(n %% 500 ==0){
      print(n)
    }
  }
  
  pred_mu = apply(out_ntk[,,k],2,mean)
  quant_mat = apply(out_ntk[,,k],2,find_quant)
  pred_low = quant_mat[1,]
  pred_high = quant_mat[2,]
  
  plot_years = c(plot_years,first_year + pred_x)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_mod = c(plot_mod,rep("tBNMA",length(pred_x)))
  plot_truth = c(plot_truth,rep("Sigmoidal Effect",length(pred_x)))
  
}

##############
# Plot Stuff #
##############

df = data.frame("Years" = plot_years,
                "Effect" = plot_pred_mu,
                "Low" = plot_pred_low,
                "High" = plot_pred_high,
                "Model" = plot_mod,
                "Truth" = plot_truth)

ggplot(data = df, aes(x = Years, y = Effect, group = Model, color = Model)) +
  geom_line() +
  facet_grid(Model ~ Truth,
             ) +
  theme(strip.placement = "outside") + 
  geom_ribbon(aes(ymin = Low, ymax = High),
              alpha = 0.2, fill = "deepskyblue4") +
  ggtitle("Credible Intervals")
