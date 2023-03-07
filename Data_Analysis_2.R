library(igraph)
library(ggplot2)
library(R2jags)
library(scales)

dat = read.csv("Data/data_tbnma.csv")

# #lets get cheeky?
# dat = dat[dat$Year < 2010,]

#######
# EDA #
#######

table(dat$Treatment)
#68 datapoints for 15 treatments
#mostly VAN and LIN
#some evidence that VAN is vecoming less effective oer time?
#so question of interest is really owh do VAN and LIN compare

mean(dat$Year[dat$Treatment == "VAN"])
mean(dat$Year[dat$Treatment == "LIN"])
#distributions of these are about the same

########
# make data structures

#find number of studies/treatments
uniq_studies = unique(dat$Study)
I = length(uniq_studies)
K = length(table(dat$Treatment))

#find years in which the studies took place
#normalize so that the first year is 0
years = NULL
for(i in 1:length(uniq_studies)){
  temp_dat = dat[dat$Study == uniq_studies[i],]
  years = c(years,temp_dat$Year[1])
}
first_year = min(years)
years = years - first_year
T = max(years)
hist(years) #these appear decently well spread out

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

#find num_treat - number of arms in each study
num_treat = rep(0,I)
for(i in 1:I){
  temp_dat = dat[which(dat$Study == study_map$Study[i]),]
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
# 
# #EDA - generate list of all treatment comparisons
# #all basically make t_mat into list of pairwise comparisons
# treat_mat = NULL
# for(i in 1:nrow(t_mat)){
#   treat_mat = cbind(treat_mat,combn(t_mat[i,1:num_treat[i]],2))
# }
# treat_mat = t(treat_mat)  
# 
# #change order to numeric
# treat_mat = treat_mat[order(treat_mat[,1]),]
# for(i in unique(treat_mat[,1])){
#   if(length(which(treat_mat[,1] == i))>1){
#     temp_treat = treat_mat[which(treat_mat[,1] == i),]
#     temp_treat[order(temp_treat[,2]),]
#     treat_mat[which(treat_mat[,1] == i),] = temp_treat
#   }
# }
# 
# uniq_treat_mat = unique(treat_mat)
# uniq_treat_mat = uniq_treat_mat[1:22,]
# 
# g2 = graph_from_data_frame(d = uniq_treat_mat, directed = FALSE)
# node_vec = unique(as.vector(uniq_treat_mat))
# size_vec = rep(0,length(node_vec))
# for(i in 1:length(node_vec)){
#   size_vec[i] = 7*sum(treat_mat == node_vec[i])^(1/3)
# }
# V(g2)$size = size_vec
# E(g2)$weight = c(13,4,2,3,2,6,4,4,2,4,2,
#                  1,1,1,1,1,1,5,1,1,
#                  1,2)
# l = layout.circle(g2)
# l = layout.fruchterman.reingold(g2)
# l = layout.davidson.harel(g2)
# plot(g2, 
#      layout = l,
#      vertex.label = treatment_map$Treatment[node_vec],
#      vertex.label.dist = c(0,0,2,-1.5,
#                            1.75,1.5,1.5,-1.5,
#                            -1.5,-1.5,1.75,1.75,
#                            -1.5,-1.5,-1.5,-1.5,
#                            -1.5,1.75,-1.5),
#      vertex.label.color = "blue",
#      vertex.color = "aquamarine",
#      edge.color = "black",
#      edge.width = E(g2)$weight)
#make time structures

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

#find how many trials/success in each study/arm
n_ik = matrix(0,nrow = I, ncol = max(num_treat))
y_ik = matrix(0,nrow = I, ncol = max(num_treat))
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
    y_ik[i,treat_ind] = temp_dat$Success[j]
  }
}

#Make data structures
plot_years = NULL
plot_pred_mu = NULL
plot_pred_low = NULL
plot_pred_high = NULL
plot_mod = NULL

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

}

##############
# Plot Stuff #
##############

df = data.frame("Years" = plot_years,
                "Effect" = plot_pred_mu,
                "Low" = plot_pred_low,
                "High" = plot_pred_high,
                "Model" = plot_mod
                )

df = df[df$Model == "tBNMA",]

ggplot(data = df, aes(x = Years, y = Effect)) +
  geom_line(linewidth = 2, color = "#619CFF") +
  # facet_wrap(~Model) +
  geom_ribbon(aes(ymin = Low, ymax = High),
              alpha = 0.2, fill = "deepskyblue4",
              color = "#619CFF") +
  ggtitle("Credible Intervals") + 
  theme(plot.title = element_text(size = 40, face = "bold"),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_blank(),
        legend.position = "none"
  )
# 
# ##########3
# # Look at probabilities VAN > LIN?
# 
# #BNMA
# sum((jags_fit$BUGSoutput$sims.list$d[,2] > 0))/length(jags_fit$BUGSoutput$sims.list$d[,2])
# 
# probs = rep(0,length(pred_x))
# for(t in 1:length(pred_x)){
#   probs[t] = sum((out_ntk[,t,2]>0))/length(out_ntk[,t,2])
# }
# plot(pred_x,probs,type = "l")
# 

#####
# Mean and CI

model = NULL
treatment = NULL
post_mean = NULL
up_bound = NULL
low_bound = NULL

##########
# Recover BNMA estimates

for(k in 1:K){
  #recover values
  d_vec = jags_fit$BUGSoutput$sims.list$d[,k]
  pred_mu = mean(d_vec)
  quant_mat = find_quant(d_vec)
  pred_low = quant_mat[1]
  pred_high = quant_mat[2]

  #put into data structures
  model = c(model,"BNMA")
  treament = c(treatment, treatment_map$Treatment[k])
  post_mean =c(post_mean,pred_mu)
  up_bound = c(up_bound,pred_high)
  low_bound = c(low_bound,pred_low)
}

##########
# Recover Meta-BNMA estimates

for(k in 1:K){

  obs_x = years_kt[[k]]
  pred_x = 20

  d_vec = jags_fit_meta$BUGSoutput$sims.list$d[,k]

  #find credible intervals for d_kt
  if(k == 2){
    b_vec = jags_fit_meta$BUGSoutput$sims.list$b[,k]
  } else {
    b_vec = rep(0,length(d_vec))
  }

  pred_vec = d_vec + b_vec*pred_x

  pred_mu = mean(pred_vec)

  quant_mat = find_quant(pred_vec)
  pred_low = quant_mat[1]
  pred_high = quant_mat[2]

  #put into data structures
  model = c(model,"Meta-BNMA")
  treament = c(treatment, treatment_map$Treatment[k])
  post_mean =c(post_mean,pred_mu)
  up_bound = c(up_bound,pred_high)
  low_bound = c(low_bound,pred_low)
}

##########
# Recover tBNMA estimates

for(k in 1:K){

  obs_x = years_kt[[k]]
  pred_x = max(years)

  #find credible intervals for d_kt
  if(k == 2){
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
      cond_Sigma = max(0.001,Sigma11 - Sigma12 %*% (Sigma22_inv)%*% Sigma21) #numeric instability

      # # A = chol(cond_Sigma)
      # A = svd(cond_Sigma)
      # A1 = A$u %*% diag(sqrt(A$d))

      out_ntk[n,,k] = sqrt(as.numeric(cond_Sigma)) %*% rnorm(length(pred_x)) + cond_mu

      # if(n %% 500 ==0){
      #   print(n)
      # }
    }

    pred_mu = mean(out_ntk[,,k])
    quant_mat = find_quant(out_ntk[,,2])
    pred_low = quant_mat[1]
    pred_high = quant_mat[2]


  } else {
    #recover values
    d_vec = jags_fit_time$BUGSoutput$sims.list$d[,k]
    pred_mu = mean(d_vec)
    quant_mat = find_quant(d_vec)
    pred_low = quant_mat[1]
    pred_high = quant_mat[2]
  }

  #put into data structures
  model = c(model,"tBNMA")
  treament = c(treatment, treatment_map$Treatment[k])
  post_mean =c(post_mean,pred_mu)
  up_bound = c(up_bound,pred_high)
  low_bound = c(low_bound,pred_low)
}


df = data.frame(Mean = post_mean,
                Low = low_bound,
                Up = up_bound,
                Treatment = as.factor(treatment_map$Treatment),
                Model = as.factor(model)
)

ggplot(df) +
  geom_errorbar(aes(x = Model, y = Mean, ymin=Low, ymax=Up, color = Model),
                position = position_dodge(0.75),
                size = 1.5) +
  geom_point( aes(x = Model, y = Mean, color = Model),
              position = position_dodge(0.75),
              size = 3, shape = 21, fill = "white") +
  facet_wrap(~Treatment) +
  coord_flip() +
  xlab("Model") +
  ylab("Treatment Effect") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "none")
