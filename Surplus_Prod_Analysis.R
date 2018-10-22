#Code for running the rjags surplus production code
#require packages needed
library(rjags)
library(dplyr)
library(mcmcplots)
library(data.table)
library(ggplot2)

#function to summarize posterior distributions 
source("../post_summ_function.R")

GOM.dat <- read.csv("1950_start_jags.csv", stringsAsFactors = F)
head(GOM.dat)

#Years is the unique years in the data set, 1984 to 2016
Years <- unique(GOM.dat$Year)

#Loop.index is the number of 1:nrows, used to loop through all data points
GOM.dat$index <- 1:nrow(GOM.dat)
Loop.index <- GOM.dat$index

#Start.vec is the row of the first year for each state
Start.vec <- GOM.dat$index[which(GOM.dat$Year.Ind==1)]

#The model will now loop through all data points except the first year for each state
Loop.index <- Loop.index[-Start.vec]

data = list(N.obs = nrow(GOM.dat),
            Year = GOM.dat$Year.Ind, 
            Catch = GOM.dat$Pounds, 
            IOA = GOM.dat$Mean_IOA,
            State = GOM.dat$Region.index,  
            N.State = length(unique(GOM.dat$Region.index)),
            Loop.index = Loop.index,
            Start.vec = Start.vec,
            N.Year = length(Years))

library(rjags)
n.chains = 5
n.adapt = 10000
n.update = 10000
n.iter = 50000
params =  c("q", "K", "r", "tau.proc", "index.new", "pvalue.fit", "msy", "fmsy", "bmsy", "total.biomass", "tau.obs", "IOA", "biomass", "p", "k.gulf", "sigma.obs", "bp.mean", "bp.var", "bp.cv")

jm.pooled = jags.model(file = "Surplus_Production_juvenile_jag_model_1.R", data = data, n.adapt = n.adapt, n.chain = n.chains)
update(jm.pooled, n.inter = n.update)
base = coda.samples(jm.pooled, variable.names = params, n.iter = n.iter, thin = 10)
#repeat update and coda.samples until model has converged, using large n.iters would cause model to crash so just do chunks at a time
update(jm.pooled, n.inter = n.update)
base = coda.samples(jm.pooled, variable.names = params, n.iter = n.iter, thin = 10)
update(jm.pooled, n.inter = n.update)
base = coda.samples(jm.pooled, variable.names = params, n.iter = n.iter, thin = 10)
update(jm.pooled, n.inter = n.update)
base = coda.samples(jm.pooled, variable.names = params, n.iter = n.iter, thin = 10)

mcmc <- as.data.frame(rbindlist(lapply(base, as.data.frame)))
save(mcmc, file = "base_MCMC.csv")




