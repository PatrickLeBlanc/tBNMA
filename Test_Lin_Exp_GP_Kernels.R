library(igraph)
library(ggplot2)
library(R2jags)
library(scales)
library(gridExtra)
library(tseries)

dat = read.csv("Data/data_date_tbnma.csv")


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
# hist(years) #these appear decently well spread out

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

###########
# GP BNMA #
###########

m_mu = 0
prec_mu = 1

#JAGS Now

#save data for jags
jags_data <- list("I", "n_ik","y_ik",
                  "num_treat","t_mat",
                  "m_mu","prec_mu",
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

print(jags_fit_time)

d_kt_post = jags_fit_time$BUGSoutput$mean$d_kt

d_kt_list = NULL
years_kt = NULL
for(k in 1:K){
  print(mean(d_kt_post[k,1:tsize[k]]))
  d_kt_list[[k]] = d_kt_post[k,1:tsize[k]]
  years_kt[[k]] = short_year_ki[k,1:tsize[k]]
}
#very slightly different results - found telehealth more effective?

post_rho = jags_fit_time$BUGSoutput$mean$rho
post_phi = jags_fit_time$BUGSoutput$mean$phi
post_psi = jags_fit_time$BUGSoutput$mean$psi
post_sb = jags_fit_time$BUGSoutput$mean$sb
post_sl = jags_fit_time$BUGSoutput$mean$sl
post_d = jags_fit_time$BUGSoutput$mean$d

jags_fit_time$BUGSoutput$summary[1:5,c(1,3,7)]

post_sd = jags_fit_time$BUGSoutput$mean$sd

plot_years = NULL
plot_pred_mu = NULL
plot_pred_low = NULL
plot_pred_high = NULL
plot_k = NULL

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

    out_ntk[n,,k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)

    if(n %% 50 ==0){
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
  plot_k = c(plot_k,rep(treatment_map$Treatment[k],length(pred_x)))
  
}

df = data.frame("Years" = plot_years,
                "Mean" = plot_pred_mu,
                "Low" = plot_pred_low,
                "High" = plot_pred_high,
                "K" = as.factor(plot_k))

# ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
#   geom_line()

p1 = ggplot(data = df, aes(x = Years, y = Mean)) +
  geom_line() +
  # facet_wrap(~K) +
  geom_ribbon(aes(ymin = Low, ymax = High),
              alpha = 0.2, fill = "deepskyblue4") +
  ggtitle("Posterior Mean Credible Intervals")  +
  scale_y_continuous(limits = c(-1.2, 1))


for(k in unique(df$K)){
  print(k)
  print( max(df$Low[df$K == k]) - min(df$High[df$K == k]))
}

plot_years = NULL
plot_pred_mu = NULL
plot_pred_low = NULL
plot_pred_high = NULL
plot_k = NULL

# for(k in c(2,8,9,10,15)){
for(k in c(2)){
  
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
      Sigma[i,i] = post_phi[k]^2 + post_sb[k]^2 + post_sl[k]^2*tot_years[i]*tot_years[i]
    } else {
      Sigma[i,i] = post_psi^2 + post_phi[k]^2 + post_sb[k]^2 + post_sl[k]^2*tot_years[i]*tot_years[i]
    }

    if(i < nrow(Sigma)){
      for(j in (i+1):ncol(Sigma)){
        Sigma[i,j] = post_phi[k]^2*exp(-post_rho[k]*abs(tot_years[i] - tot_years[j]) ) + post_sb[k]^2 + post_sl[k]^2*tot_years[i]*tot_years[j]
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

  pred_low = pred_mu - 1.96 * sqrt(diag(pred_Sigma))
  pred_high = pred_mu + 1.96 * sqrt(diag(pred_Sigma))

  plot_years = c(plot_years,first_year + pred_x)
  plot_pred_mu = c(plot_pred_mu,pred_mu)
  plot_pred_low = c(plot_pred_low,pred_low)
  plot_pred_high = c(plot_pred_high,pred_high)
  plot_k = c(plot_k,rep(treatment_map$Treatment[k],length(pred_x)))
  
}

df = data.frame("Years" = plot_years,
                "Mean" = plot_pred_mu,
                "Low" = plot_pred_low,
                "High" = plot_pred_high,
                "K" = as.factor(plot_k))

# ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
#   geom_line()

p2 = ggplot(data = df, aes(x = Years, y = Mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = Low, ymax = High),
              alpha = 0.2, fill = "deepskyblue4") + 
  ggtitle("Postior Predictive Distribution")  +
  scale_y_continuous(limits = c(-1.2, 1))

for(k in unique(df$K)){
  print(k)
  print( max(df$Low[df$K == k]) - min(df$High[df$K == k]))
}

grid.arrange(p1,p2,ncol = 2)


#Test if time-series is stationary?
#Dickey-Fuller test

plot(short_year_ki[2,]+1:tsize[2]*0.01,jags_fit_time$BUGSoutput$mean$d_kt[2,])
adf.test(jags_fit_time$BUGSoutput$mean$d_kt[2,short_year_ki[2,]])
adf.test(df$Mean)
