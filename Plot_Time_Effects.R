library(meta)
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

# #read in DeSutter et al 2012
# data_file = "Data/CD004976StatsDataOnly.rm5"
# data = read.rm5(data_file)
# 
# #add some weirdly missing data
# data$year[55] = "2007" #add weirdly missing data
# 
# #no real time scalings, at least in an lm


# #read in Spinks et al 2013
# data_file = "Data/CD000023StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #scales with time, but isn't actually suitable for an NMA
# #isn't actually a network of treatments - always 'antibiotics' vs control
# # different studies give response of different symptoms

# #read in Fisher et al 2021
# data_file = "Data/CD013437StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #wont' open

# #read in Bjelakovic et al 2017
# data_file = "Data/CD011564StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #not actually a network - just one level of treatment

# #read in franco et al 2021
# data_file = "Data/CD013656StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #won't open

# #read in Virgili et al 2018
# data_file = "Data/CD007419StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #has like 3 studies???

# #read in GAllos et al 2018
# # https://pubmed.ncbi.nlm.nih.gov/29693726/
# data_file = "Data/CD011689StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #big dataset, not terribly strong evidence of time effects?

# #read in Weibel et al 2020
# # https://pubmed.ncbi.nlm.nih.gov/33075160/
# data_file = "Data/CD012859StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #big dataset, not terribly strong evidence of time effects?

# #read in Wang et al 2019
# # https://pubmed.ncbi.nlm.nih.gov/31486548/
# data_file = "Data/CD012692StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #smaller dataset, almost evidence of time effects?

# #read in Mhaskar et al 2017
# # https://pubmed.ncbi.nlm.nih.gov/31486548/
# data_file = "Data/CD003188StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #smaller dataset, almost evidence of time effects?

# #read in Tramacere et al 2015
# # https://pubmed.ncbi.nlm.nih.gov/26384035/
# data_file = "Data/CD011381StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #maybe if you squint at it?  lots of confounding variable though

# #read in Palmer et al 2014
# # https://pubmed.ncbi.nlm.nih.gov/25486075/
# data_file = "Data/CD010590StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #this might have some time effects - if it does, it isn't linear

# #read in Iheozor-Ejiofor et al 2014
# # https://pubmed.ncbi.nlm.nih.gov/31513295/
# data_file = "Data/CD013210StatsDataOnly.rm5"
# data = read.rm5(data_file)
# #no time effects

# #read in Walsh et al 2019
# # https://pubmed.ncbi.nlm.nih.gov/30829399/
# data_file = "Data/CD007868StatsDataOnly.rm5"
# data = read.rm5(data_file)
# data$year[5] =2008
# data$year[185] = 2003
# #there might be small time effects?  e.g. Salanti et al 2009
# #worth investigating with an actual BNMA?

#read in Leibovici-Weissman et al, 2014
# https://pubmed.ncbi.nlm.nih.gov/30829399/
data_file = "Data/CD008625StatsDataOnly.rm5"
data = read.rm5(data_file)
#there might be somethign here?

#for when we have to do fancy footwork with labels

split_lab = function(y){
  names = strsplit(y, split = " ")
  out_name = names[[1]][2]
  for(i in 3:length(names[[1]])){
    out_name = paste(out_name,names[[1]][i])
  }

  return(out_name)
}

data$complab = unlist(lapply(data$grplab,split_lab))

dat_plot = data.frame("TE" = data$TE,
                "year" = data$year,
                "comp" = data$complab)

dat_sum = data_summary(dat_plot,varname = "TE",
                       groupnames = c("comp","year"))


ggplot(dat_sum, aes(x=year, y=TE, group=comp, color=comp)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=TE-sd, ymax=TE+sd), width=.2,
                position=position_dodge(0.05))+
  theme(legend.position = "none")

split_name = function(y){
  names = strsplit(y, split = " ")
  return(names[[1]][1])
}

dat_plot$comp = unlist(lapply(dat_plot$comp,split_name))

lm_mod = lm(TE ~ as.numeric(year)*comp, data = dat_plot)
summary(lm_mod)
