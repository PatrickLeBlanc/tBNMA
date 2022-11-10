#follows along with http://www.jkarreth.net/files/bayes-cph_Tutorial-JAGS.pdf
# An example model file is given in:
model.file <- system.file(package = "R2jags", "model", "schools.txt")

# data
J <- 8.0
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)

jags.data <- list("y","sd","J")
jags.params <- c("mu","sigma","theta")
jags.inits <- function(){
  list("mu"=rnorm(1),
       "sigma"=runif(1),
       "theta"=rnorm(J))
}

# Fit the model
jagsfit <- jags(data=list("y","sd","J"), inits = jags.inits,
                 jags.params, n.iter = 10, model.file = model.file)
 
print(jagsfit)
 
#fitting using R2jags

#simulate data - covariates/residuals
n.sim <- 100; set.seed(123)
x1 <- rnorm(n.sim, mean = 5, sd = 2)
x2 <- rbinom(n.sim, size = 1, prob = 0.3)
e <- rnorm(n.sim, mean = 0, sd = 1)
#make response
b1 <- 1.2
b2 <- -3.1
a <- 1.5
y <- a + b1 * x1 + b2 * x2 + e
#collect
sim.dat <- data.frame(y, x1, x2)

#fit and summarize frequentist lm
freq.mod <- lm(y ~ x1 + x2, data = sim.dat)
summary(freq.mod)

#do it in jags
jags_file = "C:/Users/psmit/Desktop/Research/NMA/jags/Tutorial/test.bug"

#define vectors of data matrix
y <- sim.dat$y
x1 <- sim.dat$x1
x2 <- sim.dat$x2
N <- nrow(sim.dat)

#save data for jags
jags_data <- list("y", "x1", "x2", "N")

#note which params to save
jags_params <- c("alpha", "beta1", "beta2")

#define inititailization values
jags_inits <- function(){
 list("alpha" = rnorm(1),
      "beta1" = rnorm(1),
      "beta2" = rnorm(1))
}

#fit model
library(R2jags)
set.seed(123)

jags_fit <- R2jags::jags(data = jags_data, 
                 inits = jags_inits,
                 parameters.to.save = jags_params,
                 n.chains = 3, n.iter = 9000,n.burnin = 1000,
                 model.file =  jags_file
                 )

#can update using this
jags_fit_upd <- update(jags_fit, n.iter=1000)
jags_fit_upd <- autojags(jags_fit_upd)

#diagnostics
print(jags_fit)
plot(jags_fit)
traceplot(jags_fit)

#more diagnostics with mcmc shenanigangs
library(mcmcplots)
jags_fit_mcmc <- as.mcmc(jags_fit)
summary(jags_fit_mcmc)
denplot(jags_fit_mcmc)
denplot(jags_fit_mcmc,parms = c("alpha","beta1","beta2"))
traplot(jags_fit_mcmc,parms = c("alpha","beta1","beta2"))

###########
# runjags #
###########

library(runjags)
#seems worse than R2jags - need to write the model in an R function using quotation marks??