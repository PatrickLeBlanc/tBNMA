#follwing https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/netwma.html
#chapter 12 and 13

#set wd
setwd("C:/Users/psmit/Desktop/Research/NMA")

library(dmetar)

data("TherapyFormats")

head(TherapyFormats)

#check for multi-arm studies
as.matrix(table(TherapyFormats$author))

#analyze using netmeta
library(netmeta)
m.netmeta <- netmeta(TE = TE,
                     seTE = seTE,
                     treat1 = treat1,
                     treat2 = treat2,
                     studlab = author,
                     data = TherapyFormats,
                     sm = "SMD",
                     fixed = TRUE,
                     random = FALSE,
                     reference.group = "cau",
                     details.chkmultiarm = TRUE,
                     sep.trts = " vs ")
summary(m.netmeta)

#calcualte total inconsistenty base don the
#full design-by-treatment interaction random effects model
decomp.design(m.netmeta)

#rank using netrank
netrank(m.netmeta)

#On to Bayes
library(gemtc)
library(rjags)

#load different data
data(TherapyFormatsGeMTC)

network <- mtc.network(data.re  = TherapyFormatsGeMTC$data,
                       treatments = TherapyFormatsGeMTC$treat.codes)
#run model
model <- mtc.model(network,
                   likelihood = "normal",
                   link = "identity",
                   linearModel = "random",
                   n.chain = 4)
mcmc1 <- mtc.run(model, n.adapt = 50, n.iter = 1000, thin = 10)
mcmc2 <- mtc.run(model, n.adapt = 5000, n.iter = 1e5, thin = 10)
# 
# #assess inconsistency using the node-splitting method
# nodesplit <- mtc.nodesplit(network, 
#                            linearModel = "random", 
#                            likelihood = "normal",
#                            link = "identity",
#                            n.adapt = 5000, 
#                            n.iter = 1e5, 
#                            thin = 10)
# summary(nodesplit)

#rank treatmetns using SUCRA
rank.probability <- rank.probability(mcmc2)
sucra <- dmetar::sucra(rank.probability, lower.is.better = TRUE)
