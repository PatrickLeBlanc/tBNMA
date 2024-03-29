  library(ggplot2)
  library(R2jags)
  library(gridExtra)
  
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
  
  set.seed(1234)
  # 
  # #set parameters
  # I = 58
  # K = 19
  # T = 19
  # 
  # #list studies and treatments
  # studies = 1:I
  # treatments = 1:K
  # 
  # #simulate years
  # years = sample(1:T,I,replace = TRUE)
  # first_year = 0
  # 
  # #simulate treatment matrix with at most 3 treatments per study
  # t_mat = matrix(0,nrow = I, ncol = 3)
  # num_treat_prob = c(1/2,1/2) #50% chance of 2 or 3 treatments
  # num_treat = rep(0,I)
  # for(i in 1:I){
  #   
  #   num_treat[i] = sample(2:3,1,replace = TRUE,prob = num_treat_prob)
  #   
  #   #simulate treatmetn 1
  #   t_mat[i,1] = sample(1:K,1)
  #   
  #   prob_treat = c(0.36,0.22,rep(0.025,17))
  #   t_mat[i,1] = sample(1:K,1,prob = prob_treat)
  #   
  #   #simulate other treatments - can't repeat treatments
  #   for(j in 1:num_treat[i]){
  #     pos_treat = 1:K
  #     pos_treat = pos_treat[-t_mat[i,1:(j-1)]]
  #     pos_prob = prob_treat[-t_mat[i,1:(j-1)]]
  #     t_mat[i,j] = sample(pos_treat,1,prob = pos_prob)
  #   }
  #   
  #   
  #   #want lower numbered treatments first
  #   t_mat[i,1:num_treat[i]] = sort(t_mat[i,1:num_treat[i]])
  #   
  # }
  
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
  m_mu = 0
  prec_mu = 1
  
  #save data for jags
  jags_data <- list("I","n_ik","y_ik",
                    "num_treat", "t_mat",
                    "m_mu","prec_mu",
                    "K")
  
  #note which params to save
  jags_params <- c("d","sd",
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
  
  print(jags_fit)
  
  post_d = jags_fit$BUGSoutput$mean$d
  plot(d,post_d, xlab = "True d_{1k}", ylab = "Posterior Mean d_{1k}", main = "BNMA")
  lines(-5:5,-5:5)
  #seems to give accurate inference to mean
  #bad at estimating end-time effect
  
  
  # #calculate sucra
  # 
  # #extract mcmc chains
  # d_mat = jags_fit$BUGSoutput$sims.matrix[,1:5]
  # #calc mc rank prob for each treatment
  # bnma_rank_prob = matrix(0,nrow = K, ncol = K)
  # for(n in 1:nrow(d_mat)){
  #   temp = order(d_mat[n,],decreasing = TRUE)
  #   for(i in 1:length(temp)){
  #     bnma_rank_prob[temp[i],i] = bnma_rank_prob[temp[i],i] + 1
  #   }
  # }
  # for(k in 1:K){
  #   bnma_rank_prob[k,] = bnma_rank_prob[k,]/sum(bnma_rank_prob[k,])
  # }
  # #calc sucra for each treatment
  # bnma_sucra = rep(0,K)
  # for(k in 1:K){
  #   for(i in 1:(K-1)){
  #     bnma_sucra[k] = bnma_sucra[k] + sum(bnma_rank_prob[k,1:i])
  #   }
  # }
  # bnma_sucra = bnma_sucra/(K-1)
  
  
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
  #                   "num_treat","t_mat",
  #                   "m_mu","prec_mu",
  #                   "K",
  #                   "years")
  # 
  # #note which params to save
  # jags_params <- c("b","d","sd","z_k")
  # 
  # #define inititailization values
  # jags_inits <- function(){
  #   list("sd" = runif(1,0,5)
  #   )
  # }
  # 
  # dir = getwd()
  # jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm_Time_Z.bug")
  # 
  # #fit model
  # set.seed(123)
  # 
  # jags_fit_meta = R2jags::jags(data = jags_data,
  #                              inits = jags_inits,
  #                              parameters.to.save = jags_params,
  #                              n.chains = 3, n.iter = 9000,n.burnin = 1000,
  #                              model.file =  jags_file
  # )
  # 
  # print(jags_fit_meta)
  # 
  # # #calculate sucra
  # # 
  # # #extract mcmc chains
  # # d_mat = matrix(0,nrow = nrow(jags_fit$BUGSoutput$sims.matrix), ncol = K)
  # # #find each d_k^T (might be easier to do this in JAGS in future?)
  # # for(k in 1:K){
  # #   d_mat[,k] = jags_fit_meta$BUGSoutput$sims.matrix[,K + k] + T*jags_fit_meta$BUGSoutput$sims.matrix[,k]
  # # }
  # # 
  # # #calc mc rank prob for each treatment
  # # bnma_meta_rank_prob = matrix(0,nrow = K, ncol = K)
  # # for(n in 1:nrow(d_mat)){
  # #   temp = order(d_mat[n,],decreasing = TRUE)
  # #   for(i in 1:length(temp)){
  # #     bnma_meta_rank_prob[temp[i],i] = bnma_meta_rank_prob[temp[i],i] + 1
  # #   }
  # # }
  # # for(k in 1:K){
  # #   bnma_meta_rank_prob[k,] = bnma_meta_rank_prob[k,]/sum(bnma_meta_rank_prob[k,])
  # # }
  # # #calc sucra for each treatment
  # # bnma_meta_sucra = rep(0,K)
  # # for(k in 1:K){
  # #   for(i in 1:(K-1)){
  # #     bnma_meta_sucra[k] = bnma_meta_sucra[k] + sum(bnma_meta_rank_prob[k,1:i])
  # #   }
  # # }
  # # bnma_meta_sucra = bnma_meta_sucra/(K-1)
  # 
  # ##################
  # # BNMA Sigmoidal #
  # ##################
  # 
  # pi = 0.1
  # m_mu = 0
  # prec_mu = 1
  # 
  # jags_data <- list("I","n_ik","y_ik",
  #                   "m_mu","prec_mu",
  #                   "num_treat", "t_mat",
  #                   "K",
  #                   "short_year_ki",
  #                   "tsize","TS",
  #                   "lts_ind", "T",
  #                   "pi")
  # 
  # #note which params to save
  # jags_params <- c("z_k","d",
  #                  "height_k",
  #                  "scale_k", "center_k",
  #                  "psi")
  # 
  # #define inititailization values
  # jags_inits <- function(){
  #   list("sd" = runif(1,0,5)
  #   )
  # }
  # 
  # dir = getwd()
  # jags_file = paste0(dir,"/jags/BNMA_Like_Bin_Trial_Multi_Arm_Time_Sigmoidal_Z.bug")
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
  # ggplot(data = df, aes(x = Years, y = Mean, group = K, color = K)) +
  #   geom_line() +
  #   facet_grid(~K) +
  #   geom_ribbon(aes(ymin = Low, ymax = High),
  #               alpha = 0.2, fill = "deepskyblue4")
  # 
  # 
  # # post_delta = jags_fit_sig$BUGSoutput$mean$delta
  # # plot(delta,post_delta)
  # # lines(-5:5,-5:5)
  # #again, accurate inference to mean
  # 
  # 
  # # #calculate sucra
  # # 
  # # #extract mcmc chains
  # # d_mat = matrix(0,nrow = nrow(jags_fit_sig$BUGSoutput$sims.matrix), ncol = K)
  # # #find each d_k^T (might be easier to do this in JAGS in future?)
  # # for(k in 1:K){
  # # 
  # #   d_k = jags_fit_sig$BUGSoutput$sims.matrix[,K + k]
  # #   a_k = jags_fit_sig$BUGSoutput$sims.matrix[,2*K + 1 + k]
  # #   b_k = jags_fit_sig$BUGSoutput$sims.matrix[,3*K + 2 + k]
  # #   c_k = jags_fit_sig$BUGSoutput$sims.matrix[,k]
  # #   psi = jags_fit_sig$BUGSoutput$sims.matrix[,17]
  # # 
  # #   z_k = jags_fit_sig$BUGSoutput$sims.matrix[,22 + k]
  # # 
  # #   for(it in 1:nrow(d_mat)){
  # #     sig_draw = rnorm(1,(d_k[it] - a_k[it]) + 2*(d_k[it] + a_k[it])/(1+exp(-b_k[it]*(T - c_k[it]))),psi^2)
  # #     d_mat[it,k] = z_k[it]*sig_draw + (1-z_k[it])*d_k[it]
  # #   }
  # # 
  # # }
  # # 
  # # #calc mc rank prob for each treatment
  # # bnma_sig_rank_prob = matrix(0,nrow = K, ncol = K)
  # # for(n in 1:nrow(d_mat)){
  # #   temp = order(d_mat[n,],decreasing = TRUE)
  # #   for(i in 1:length(temp)){
  # #     bnma_sig_rank_prob[temp[i],i] = bnma_sig_rank_prob[temp[i],i] + 1
  # #   }
  # # }
  # # for(k in 1:K){
  # #   bnma_sig_rank_prob[k,] = bnma_sig_rank_prob[k,]/sum(bnma_sig_rank_prob[k,])
  # # }
  # # #calc sucra for each treatment
  # # bnma_sig_sucra = rep(0,K)
  # # for(k in 1:K){
  # #   for(i in 1:(K-1)){
  # #     bnma_sig_sucra[k] = bnma_sig_sucra[k] + sum(bnma_sig_rank_prob[k,1:i])
  # #   }
  # # }
  # # bnma_sig_sucra = bnma_sig_sucra/(K-1)
  
  ############
  # GP-BNMA #
  ###########
  
  m_mu = 0
  prec_mu = 1
  
  #JAGS Now
  
  mu = mu_i
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
  
  find_quant = function(x){
    return(quantile(x,c(0.025,0.975)))
  }
  
  # for(k in c(2,8,9,10,15)){
  for(k in c(2)){
    
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
    ggtitle("Posterior Predictive Distribution")  +
    scale_y_continuous(limits = c(-1.2, 1))
  
  for(k in unique(df$K)){
    print(k)
    print( max(df$Low[df$K == k]) - min(df$High[df$K == k]))
  }
  
  grid.arrange(p1,p2,ncol = 2)
  # 
  # # #extract mcmc chains
  # # d_mat = matrix(0,nrow = nrow(jags_fit$BUGSoutput$sims.matrix), ncol = K)
  # # #find each d_k^T (might be easier to do this in JAGS in future?)
  # # for(k in 2:K){
  # #   #find
  # #   obs_x = years_kt[[k]]
  # #   y = d_kt_list[[k]]
  # #   pred_x =T
  # #   mu1 = rep(post_d[k],length(pred_x))
  # #   mu2 = rep(post_d[k],length(obs_x))
  # #   one_ind = 1:length(pred_x)
  # #   two_ind = (length(pred_x)+1):(length(pred_x) + length(obs_x))
  # # 
  # #   tot_years = c(pred_x,obs_x)
  # # 
  # #   phi_vec = jags_fit_time$BUGSoutput$sims.list$phi[,k]
  # #   rho_vec = jags_fit_time$BUGSoutput$sims.list$rho[,k]
  # #   psi_vec = jags_fit_time$BUGSoutput$sims.list$psi
  # # 
  # # 
  # #   for(n in 1:nrow(d_mat)){
  # # 
  # #     #find covariance
  # #     Sigma = matrix(0,nrow = length(tot_years),ncol = length(tot_years))
  # #     for(i in 1:nrow(Sigma)){
  # # 
  # #       if(i <= length(pred_x)){
  # #         Sigma[i,i] = phi_vec[n]^2
  # #       } else {
  # #         Sigma[i,i] = psi_vec[n]^2 + phi_vec[n]^2
  # #       }
  # # 
  # #       if(i < nrow(Sigma)){
  # #         for(j in (i+1):ncol(Sigma)){
  # #           Sigma[i,j] = phi_vec[n]^2*exp(-rho_vec[n]*(tot_years[i] - tot_years[j])^2 )
  # #           Sigma[j,i] = Sigma[i,j]
  # #         }
  # #       }
  # # 
  # #     }
  # # 
  # #     Sigma11 = Sigma[one_ind,one_ind]
  # #     Sigma12 = Sigma[one_ind,two_ind]
  # #     Sigma22 = Sigma[two_ind,two_ind]
  # #     Sigma22_inv = solve(Sigma22)
  # #     Sigma21 = Sigma[two_ind,one_ind]
  # # 
  # #     d_mat[n,k] = mu1 + Sigma12 %*% Sigma22_inv %*%(y - mu2)
  # # 
  # #     if(n %% 50 ==0){
  # #       print(n)
  # #     }
  # #   }
  # # }
  # # #calc mc rank prob for each treatment
  # # gpbnma_rank_prob = matrix(0,nrow = K, ncol = K)
  # # for(n in 1:nrow(d_mat)){
  # #   temp = order(d_mat[n,],decreasing = TRUE)
  # #   for(i in 1:length(temp)){
  # #     gpbnma_rank_prob[temp[i],i] = gpbnma_rank_prob[temp[i],i] + 1
  # #   }
  # # }
  # # for(k in 1:K){
  # #   gpbnma_rank_prob[k,] = gpbnma_rank_prob[k,]/sum(gpbnma_rank_prob[k,])
  # # }
  # # #calc sucra for each treatment
  # # gpbnma_sucra = rep(0,K)
  # # for(k in 1:K){
  # #   for(i in 1:(K-1)){
  # #     gpbnma_sucra[k] = gpbnma_sucra[k] + sum(gpbnma_rank_prob[k,1:i])
  # #   }
  # # }
  # # bnma_gp_sucra = gpbnma_sucra/(K-1)
  # 
  # ############
  # # Plots ?! #
  # ############
  # 
  # greater_than_zero = function(x){
  #   return(sum((x>0))/length(x))
  # }
  # 
  # #####
  # # Mean and CI
  # 
  # post_prob =NULL
  # model = NULL
  # treatment = NULL
  # post_mean = NULL
  # up_bound = NULL
  # low_bound = NULL
  # line = NULL
  # 
  # post_d = jags_fit$BUGSoutput$sims.list$d
  # temp = apply(post_d,2,greater_than_zero)
  # post_prob = c(post_prob,temp)
  # post_mean = c(post_mean,apply(post_d,2,mean))
  # quant_mat = apply(post_d,2,find_quant)
  # low_bound = c(low_bound,quant_mat[1,])
  # up_bound = c(up_bound,quant_mat[2,])
  # model = c(model,rep("BNMA",length(temp)))
  # treatment = c(treatment,treatment_map$Treatment)
  # line = c(line,rep(0,length(temp)))
  # 
  # post_d = jags_fit_meta$BUGSoutput$sims.list$d
  # temp = apply(post_d,2,greater_than_zero)
  # post_prob = c(post_prob,temp)
  # post_mean = c(post_mean,apply(post_d,2,mean))
  # quant_mat = apply(post_d,2,find_quant)
  # low_bound = c(low_bound,quant_mat[1,])
  # up_bound = c(up_bound,quant_mat[2,])
  # model = c(model,rep("Meta-BNMA",length(temp)))
  # treatment = c(treatment,treatment_map$Treatment)
  # line = c(line,rep(0,length(temp)))
  # 
  # post_d = jags_fit_sig$BUGSoutput$sims.list$d
  # temp = apply(post_d,2,greater_than_zero)
  # post_prob = c(post_prob,temp)
  # post_mean = c(post_mean,apply(post_d,2,mean))
  # quant_mat = apply(post_d,2,find_quant)
  # low_bound = c(low_bound,quant_mat[1,])
  # up_bound = c(up_bound,quant_mat[2,])
  # model = c(model,rep("Sig-BNMA",length(temp)))
  # treatment = c(treatment,treatment_map$Treatment)
  # line = c(line,rep(0,length(temp)))
  # 
  # post_d = jags_fit_time$BUGSoutput$sims.list$d
  # temp = apply(post_d,2,greater_than_zero)
  # post_prob = c(post_prob,temp)
  # post_mean = c(post_mean,apply(post_d,2,mean))
  # quant_mat = apply(post_d,2,find_quant)
  # low_bound = c(low_bound,quant_mat[1,])
  # up_bound = c(up_bound,quant_mat[2,])
  # model = c(model,rep("GP-BNMA",length(temp)))
  # treatment = c(treatment,treatment_map$Treatment)
  # line = c(line,rep(0,length(temp)))
  # 
  # post_prob = c(post_prob,rep(0,length(temp)))
  # post_mean = c(post_mean,d)
  # low_bound = c(low_bound,d)
  # up_bound = c(up_bound,d)
  # model = c(model,rep("Truth",length(temp)))
  # treatment = c(treatment, treatment_map$Treatment)
  # line = c(line,rep(0,length(temp)))
  # 
  # df = data.frame(Probability = post_prob,
  #                 Mean = post_mean,
  #                 Low = low_bound,
  #                 Up = up_bound,
  #                 Treatment = as.factor(treatment_map$Treatment),
  #                 Model = as.factor(model),
  #                 Line = line)
  # 
  # ggplot(df) + 
  #   geom_errorbar(aes(x = Model, y = Mean, ymin=Low, ymax=Up, color = Model),
  #                 position = position_dodge(0.75),
  #                 size = 1.5) +
  #   geom_point( aes(x = Model, y = Mean, color = Model), 
  #               position = position_dodge(0.75), 
  #               size = 3, shape = 21, fill = "white") +
  #   facet_wrap(~Treatment) + 
  #   coord_flip() + 
  #   xlab("Model") + 
  #   ylab("Treatment Effect") + 
  #   theme(axis.text=element_text(size=14),
  #         axis.title=element_text(size=16,face="bold"),
  #         legend.position = "none")
  # 
  # # 
  # # #####
  # # # Boxplots
  # # 
  # # model = NULL
  # # treatment = NULL
  # # post_mean = NULL
  # # up_bound = NULL
  # # low_bound = NULL
  # # line = NULL
  # # 
  # # post_d = jags_fit$BUGSoutput$sims.list$d
  # # for(k in 1:ncol(post_d)){
  # #   post_mean = c(post_mean,post_d[,k])
  # #   model = c(model,rep("BNMA",nrow(post_d)))
  # #   treatment = c(treatment,rep(treatment_map$Treatment[k],nrow(post_d)))
  # # }
  # # 
  # # 
  # # post_d = jags_fit_meta$BUGSoutput$sims.list$d
  # # for(k in 1:ncol(post_d)){
  # #   post_mean = c(post_mean,post_d[,k])
  # #   model = c(model,rep("Meta-BNMA",nrow(post_d)))
  # #   treatment = c(treatment,rep(treatment_map$Treatment[k],nrow(post_d)))
  # # }
  # # 
  # # post_d = jags_fit_sig$BUGSoutput$sims.list$d
  # # for(k in 1:ncol(post_d)){
  # #   post_mean = c(post_mean,post_d[,k])
  # #   model = c(model,rep("Sig-BNMA",nrow(post_d)))
  # #   treatment = c(treatment,rep(treatment_map$Treatment[k],nrow(post_d)))
  # # }
  # # 
  # # post_d = jags_fit_time$BUGSoutput$sims.list$d
  # # for(k in 1:ncol(post_d)){
  # #   post_mean = c(post_mean,post_d[,k])
  # #   model = c(model,rep("GP-BNMA",nrow(post_d)))
  # #   treatment = c(treatment,rep(treatment_map$Treatment[k],nrow(post_d)))
  # # }
  # # 
  # # df = data.frame(Mean = post_mean,
  # #                 Treatment = as.factor(treatment_map$Treatment),
  # #                 Model = as.factor(model))
  # # 
  # # ggplot(df,aes(x = Model, y = Mean)) + 
  # #   geom_boxplot() + 
  # #   facet_wrap(~Treatment)
  # # 
  # # 
