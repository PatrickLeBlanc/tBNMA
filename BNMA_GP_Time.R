library(ggplot2)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

setwd("C:/Users/psmit/Desktop/Research/NMA")

# # for 
# 
# #read in Walsh et al 2019
# # https://pubmed.ncbi.nlm.nih.gov/30829399/
# data_file = "Data/CD007868StatsDataOnly.rm5"
# data = read.rm5(data_file)
# data$year[5] =2008
# data$year[185] = 2003
# 
# #just do a BNMA first - no time effects
# 
# #do normally distributed responses first
# norm_rows = !is.na(data$mean.e)
# 
# studlabs = data$studlab[norm_rows]
# TE = data$TE[norm_rows]
# complab = data$complab[norm_rows]
# year = data$year
# 
# #find number of unique studies
# unique_studlabs = unique(studlabs)
# I = length(unique_studlabs)
# 
# #find number of arms in each study
# for(i in 1:I){
#   sum(studlabs == unique_studlabs[i])
# }


library(dmetar)

data("TherapyFormats")

#find if anything appears more than once (don't want to deal with this yet)
temp = table(TherapyFormats$author)
temp[temp > 1]

rows = TherapyFormats$author != "Breiman, 2001"
TherapyFormats = TherapyFormats[rows,]

studies = unique(TherapyFormats$author)
studies_map = data.frame("index" = 1:length(studies),
                         "study" = studies)
I = length(studies)

treatments = rev(unique(c(TherapyFormats$treat1,TherapyFormats$treat2)))
treatment_map = data.frame("index" = 1:7,
                           "treatment" = treatments)
K = length(treatments)

#get treatments in vector format
replace_treat = function(x){
  which(treatment_map$treatment == x)
}

t2 = unlist(lapply(TherapyFormats$treat1,replace_treat))
t1 = unlist(lapply(TherapyFormats$treat2,replace_treat))

y = TherapyFormats$TE

se = TherapyFormats$seTE

#add year effects?
extract_year = function(y){
  names = strsplit(y, split = " ")
  return(names[[1]][length(names[[1]])])
}
years = unlist(lapply(studies,extract_year))
years = as.numeric(years)
first_year = min(years)
years = years - min(years)
uniq_years = unique(years)
T = max(years)


#find how many datapoints for each treatment
tsize = rep(0,K)
for(k in 1:K){
  for(i in 1:I){
    if(t1[i] == k | t2[i] == k){
      tsize[k] = tsize[k] + 1
    } 
  }
}

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
# 
# #save data for jags
# jags_data <- list("I","se","y",
#                   "t2","t1","K")
# 
# #note which params to save
# jags_params <- c("d","sd")
# 
# #define inititailization values
# jags_inits <- function(){
#   list("sd" = runif(1,0,5)
#   )
# }
# 
# jags_file = "C:/Users/psmit/Desktop/Research/NMA/jags/BNMA_Like_Norm_Trial_Two_Arm.bug"
# 
# #fit model
# library(R2jags)
# set.seed(123)
# 
# jags_fit <- R2jags::jags(data = jags_data, 
#                          inits = jags_inits,
#                          parameters.to.save = jags_params,
#                          n.chains = 3, n.iter = 9000,n.burnin = 1000,
#                          model.file =  jags_file
# )
# 
# print(jags_fit)
# 
# d_post = data.frame("index" = 1:7,
#                     "treatment" = treatments,
#                     "d" = jags_fit$BUGSoutput$mean$d)

########
# Time #
########

TS = max(tsize)

#save data for jags
jags_data <- list("I","se","y",
                  "t2","t1","K",
                  "short_year_ki", 
                  "tsize","TS",
                  "lts_ind")

#note which params to save
jags_params <- c("d","d_kt","sd",
                 "phi","rho","psi")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/NMA/jags/BNMA_Like_Norm_Trial_Two_Arm_Time_GP.bug"

#fit model
library(R2jags)
set.seed(123)

jags_fit_time <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)

print(jags_fit_time)

d_post_time = data.frame("index" = 1:7,
                    "treatment" = treatments,
                    "d" = jags_fit_time$BUGSoutput$mean$d)

d_kt_post = jags_fit_time$BUGSoutput$mean$d_kt
matplot(t(d_kt_post),type="l")

d_kt_list = NULL
years_kt = NULL
for(k in 1:K){
  print(mean(d_kt_post[k,1:tsize[k]]))
  d_kt_list[[k]] = d_kt_post[k,1:tsize[k]]
  years_kt[[k]] = short_year_ki[k,1:tsize[k]]
}
#very slightly different results - found telehealth more effective?

plot(years_kt[[1]],d_kt_list[[1]],type="l",ylim = c(-1.25,1))
for(k in 2:K){
  new_ind = order(years_kt[[k]])
  lines(years_kt[[k]][new_ind],d_kt_list[[k]][new_ind],col = k)
}

post_rho = jags_fit_time$BUGSoutput$mean$rho
post_phi = jags_fit_time$BUGSoutput$mean$phi
post_psi = jags_fit_time$BUGSoutput$mean$psi
post_d = jags_fit_time$BUGSoutput$mean$d

post_sd = jags_fit_time$BUGSoutput$mean$sd

plot_years = NULL
plot_pred_mu = NULL
plot_pred_low = NULL
plot_pred_high = NULL
plot_k = NULL

for(k in 2:K){
  obs_x = years_kt[[k]]
  y = d_kt_list[[k]]
  pred_x = 0:(10*T)/10
  
  mu1 = rep(post_d[k],length(pred_x))
  mu2 = rep(post_d[k],length(obs_x))
  
  
  #find covariance
  tot_years = c(pred_x,obs_x)
  Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
  for(i in 1:nrow(Sigma)){
    
    if(i <= length(pred_x)){
      Sigma[i,i] = post_phi^2
    } else {
      Sigma[i,i] = post_psi^2 + post_phi^2
    }

    if(i < nrow(Sigma)){
      for(j in (i+1):ncol(Sigma)){
        Sigma[i,j] = post_phi^2*exp(-post_rho*(tot_years[i] - tot_years[j])^2 )
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
  
  pred_mu = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
  pred_Sigma = Sigma11 - Sigma12 %*% Sigma22_inv %*% Sigma21
  
  pred_low = pred_mu - 1.95 * sqrt(diag(pred_Sigma))
  pred_high = pred_mu + 1.95 * sqrt(diag(pred_Sigma))
  
  plot_years = c(plot_years,pred_x)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_k = c(plot_k,rep(treatments[k],length(pred_x)))
}

plot_years = plot_years + first_year
df = data.frame("Years" = plot_years,
                "Mean" = plot_pred_mu,
                "Low" = plot_pred_low,
                "High" = plot_pred_high,
                "K" = plot_k)

ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
  geom_line() 

ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
  geom_line() + 
  facet_grid(~K) + 
  geom_ribbon(aes(ymin = Low, ymax = High), 
              alpha = 0.2, fill = "deepskyblue4")
