  library(igraph)
  library(ggplot2)
  library(R2jags)
  library(scales)
  
  dat = read.csv("Data/data_tbnma.csv")
  
  #lets get cheeky?
  dat = dat[dat$Year <2009,]
  
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
  years = dat$Year[2*1:I]
  first_year = min(years)
  years = years - first_year
  T = max(years)
  hist(years) #these appear decently well spread out
  
  #list studies and treatments
  studies = 1:I
  treatments = 1:K
  
  #make map between treatment number and treatmetns
  treatment_map = data.frame("Number" = treatments,
                             "Treatment" = unique(dat$Treatment))
  
  # #make lin 1 and van 2
  # treatment_map$Treatment[1] = "LIN"
  # treatment_map$Treatment[2] = "VAN"
  
  #find num_treat
  num_treat = rep(0,I)
  for(i in 1:I){
    temp_dat = dat[which(dat$Study == uniq_studies[i]),]
    num_treat[i] = nrow(temp_dat)
  }
  
  #find treatment matrix and num_treat
  t_mat = matrix(0,nrow = I, ncol = max(num_treat))
  for(i in 1:I){
    
    temp_dat = dat[which(dat$Study == uniq_studies[i]),]
    
    for(j in 1:num_treat[i]){
      t_mat[i,j] = which(treatment_map$Treatment == temp_dat$Treatment[j])
    }
    #want lower numbered treatments first
    t_mat[i,1:num_treat[i]] = sort(t_mat[i,1:num_treat[i]])
  }
  # 
  # #EDA - generate list of all treatment comparisons
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
  #   size_vec[i] = 5*sum(treat_mat == node_vec[i])^(1/3)
  # }
  # V(g2)$size = size_vec
  # E(g2)$weight = c(13,4,2,3,2,6,4,4,2,4,2,
  #                  1,1,1,1,1,1,5,1,1,
  #                  1,2)
  # l = layout.circle(g2)
  # radian.rescale <- function(x, start=0, direction=1) {
  #   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  #   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  # }
  # lab.locs <- radian.rescale(x=1:max(node_vec), direction=-1, start=0)
  # plot(g2, 
  #      layout = l,
  #      vertex.label = treatment_map$Treatment[node_vec],
  #      # vertex.label.degree=lab.locs,
  #      vertex.label.dist = c(0,0,1.5,1.5,
  #                            1.5,1.5,-1.5,-1.5,
  #                            -1.5,-1.5,-1.5,-1.5,
  #                            -1.5,-1.5,-1.5,-1.5,
  #                            -1.5,-1.5,1.5),
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
  long_year_ki = matrix(NA,nrow = K, ncol = I)
  for(k in 1:K){
    for(i in 1:I){
      if(sum(t_mat[i,] ==k) > 0 ){
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
  
  #find how many trials/success in each study/arm
  n_ik = matrix(0,nrow = I, ncol = max(num_treat))
  y_ik = matrix(0,nrow = I, ncol = max(num_treat))
  count = 1
  for(i in 1:I){
    for(j in 1:num_treat[i]){
      n_ik[i,j] = dat$Trials[count]
      y_ik[i,j] = dat$Success[count]
      count = count + 1
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
  jags_params <- c("d","sd")
  
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
  
  post_d = jags_fit$BUGSoutput$sims.list$d
  greater_than_zero = function(x){
    return(sum((x>0))/length(x))
  }
  post_prob = apply(post_d,2,greater_than_zero)
  df = data.frame(Probability = post_prob,
                  Treatment = as.factor(treatment_map$Treatment))
  
  ggplot(data = df, aes(x = Treatment, y = Probability)) +
    geom_bar(stat = "identity", fill = "navy", color = "black") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))