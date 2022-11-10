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

t1 = unlist(lapply(TherapyFormats$treat1,replace_treat))
t2 = unlist(lapply(TherapyFormats$treat2,replace_treat))

y = TherapyFormats$TE

se = TherapyFormats$seTE


#save data for jags
jags_data <- list("I","se","y",
                  "t2","t1","K")

#note which params to save
jags_params <- c("d","sd")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/NMA/jags/BNMA_Like_Norm_Trial_Two_Arm.bug"

#fit model
library(R2jags)
set.seed(123)

jags_fit <- R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 9000,n.burnin = 1000,
                         model.file =  jags_file
)

print(jags_fit)

d_post = data.frame("index" = 1:7,
                    "treatment" = treatments,
                    "d" = jags_fit$BUGSoutput$mean$d)


#add year effects?
extract_year = function(y){
  names = strsplit(y, split = " ")
  return(names[[1]][length(names[[1]])])
}

years = unlist(lapply(studies,extract_year))
years = as.numeric(years)
years = years - min(years)

#do some EDA
dat_plot = data.frame("TE" = y,
                      "year" = years,
                      "comp" = paste(t1,t2))
#dat_plot = dat_plot[1:61,]

dat_sum = data_summary(dat_plot,varname = "TE",
                       groupnames = c("comp","year"))

ggplot(dat_sum, aes(x=year, y=TE, group=comp, color=comp)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=TE-sd, ymax=TE+sd), width=.2,
                position=position_dodge(0.05))

lm_mod = lm(TE ~ as.numeric(year)*comp, data = dat_plot)
summary(lm_mod)
#doesn't look like years have an effect

#impelment a time model
#save data for jags
jags_data <- list("I","se","y",
                  "t2","t1","K",
                  "years")

#note which params to save
jags_params <- c("b","d","sd")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/NMA/jags/BNMA_Like_Norm_Trial_Two_Arm_Time.bug"

#fit model
set.seed(123)

jags_fit_meta = R2jags::jags(data = jags_data, 
                         inits = jags_inits,
                         parameters.to.save = jags_params,
                         n.chains = 3, n.iter = 20000,n.burnin = 10000,
                         model.file =  jags_file
)

print(jags_fit_meta)

d_post_time = data.frame("index" = 1:7,
                    "treatment" = treatments,
                    "d" = jags_fit_time$BUGSoutput$mean$d)

b_post_time = data.frame("index" = 1:7,
                    "treatment" = treatments,
                    "b" = jags_fit_time$BUGSoutput$mean$b)

#figure out where studies of treatmetns with largest effects are located
hist(years[t2 == 1 | t1 == 1])
hist(years[t2 == 2 | t1 == 2])
hist(years[t2 == 3 | t1 == 3])
hist(years[t2 == 4 | t1 == 4])
hist(years[t2 == 5 | t1 == 5])
hist(years[t2 == 6 | t1 == 6])
hist(years[t2 == 7 | t1 == 7])

# delta_post_time =matrix(0,nrow = K,ncol = K)
# for(i in 1:K){
#   for(j in i:K){
#     delta_post_time[i,j] = d_post_time$d[j] - d_post_time$d[i]
#   }
# }
# 
# beta_post_time = matrix(0,nrow = K,ncol = K)
# for(i in 1:K){
#   for(j in i:K){
#     beta_post_time[i,j] = b_post_time$b[j] - b_post_time$b[i]
#   }
# }
# 
# beta_lwr_time = matrix(0,nrow = K,ncol = K)
# beta_upr_time = matrix(0,nrow = K,ncol = K)
# for(i in 1:K){
#   for(j in i:K){
#     vec = jags_fit_time$BUGSoutput$sims.array[,,j] - jags_fit_time$BUGSoutput$sims.array[,,i]
#     beta_lwr_time[i,j] = quantile(vec,0.025)
#     beta_upr_time[i,j] = quantile(vec,0.75)
#   }
# }
# 
# line_mat = matrix(0,nrow = K*K,ncol= max(years))
# time = 1:max(years)
# for(i in 1:K){
#   for(j in 1:K){
#     line_mat[K*(i-1)+j,] = delta_post_time[i,j] + beta_post_time[i,j]*time
#   }
# }
# matplot(t(line_mat),type= "l")

#be more constrained
# if t2 is waitlist or care as usual (i think this is placeob?), do time effects
time_ind = as.numeric(t2 == 1) + 1


#impelment a time model
#save data for jags
jags_data <- list("I","se","y",
                  "t2","t1","K",
                  "years", "time_ind")

#note which params to save
jags_params <- c("beta","d","sd")

#define inititailization values
jags_inits <- function(){
  list("sd" = runif(1,0,5)
  )
}

jags_file = "C:/Users/psmit/Desktop/Research/NMA/jags/BNMA_Like_Norm_Trial_Two_Arm_Time_Placebo.bug"

#fit model
set.seed(123)

jags_fit_time = R2jags::jags(data = jags_data, 
                             inits = jags_inits,
                             parameters.to.save = jags_params,
                             n.chains = 3, n.iter = 20000,n.burnin = 10000,
                             model.file =  jags_file
)

print(jags_fit_time)

d_post_placebo = data.frame("index" = 1:7,
                         "treatment" = treatments,
                         "d" = jags_fit_time$BUGSoutput$mean$d)

#very similar to no time effects - 
#whatever effects we're picking up above are coming from somewher else
#most likley between gsh and tel

#87% chance of being positive in posterior