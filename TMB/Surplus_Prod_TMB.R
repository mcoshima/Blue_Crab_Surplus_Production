########################################################
#set up the data

setwd("~/Dropbox/Blue Crab Surplus Production/Code/Blue_Crab_Surplus_Production")

dat <- read.csv("1950_start_jags.csv")
head(dat)
setwd("~/Dropbox/Blue Crab Surplus Production/Code/Blue_Crab_Surplus_Production/TMB")

library(dplyr)
library(tidyr)
sub.dat <- dat %>% filter(Year.Ind > 34)
sub.dat <- sub.dat %>% 
  group_by(Region.index) %>% 
  mutate(grouped_id = row_number())
sub.dat
catch <- sub.dat %>% select(grouped_id, Region.index, Pounds) %>% spread(key = Region.index, value = Pounds) %>%  select(-grouped_id) 
catch

index <- sub.dat %>% select(grouped_id, Region.index, Mean_IOA) %>% spread(key = Region.index, value = Mean_IOA) %>%  select(-grouped_id)
index
write.csv(index, file = "index.mat.csv", row.names = F)
write.csv(catch, file = "catch.mat.csv", row.names = F)
#################################################################################
setwd("~/Dropbox/Blue Crab Surplus Production/Code/Blue_Crab_Surplus_Production/TMB")
catch <- read.csv("catch.mat.csv", header = T)
index <- read.csv("index.mat.csv", header = T)
library(TMB)
compile("surplus_prod_v1_GA.cpp")
dyn.load(dynlib("surplus_prod_v1_GA"))

data <- list(Pounds = as.matrix(catch[,c(1,2)]),
             IOA = as.matrix(index[,c(1,2)]),
             n = nrow(catch),
             NRegion = 2)
#rep(-9.0, data$NRegion)
parameters <- list(logq = rep(-9.0,2),
                   logk = 20.0,
                   logr = -1.1,
                   logSigmaO = 0.0,
                   dummy = 0,
                   logSig_q = 0.0,
                   logq_mu = 1)

#Map for testing
map<-list(logq=rep(factor(NA),2),logk=factor(NA),logr=factor(NA),logSigmaO = factor(NA), logSig_q = factor(NA), logq_mu = factor(NA))
#Map for fitting model, gets rid of dummy
map <- list(dummy = factor(NA))
#Setting bounds
#logq = 0,1, logk = c(16592354,Inf), logr = c(0,2), logSigma0 = (-10,10), dummy = c(-Inf,Inf)
#lowbnd = c(rep(1,12), max(catch), 1, 0.007)
#uppbnd = c(rep(2.72,12), Inf, 7.3, 148.4)
sp <- MakeADFun(data, parameters,DLL="surplus_prod_v1_GA",silent=T, random = "logq",map = map)
fit <- nlminb(sp$par, sp$fn, sp$gr) #, lower=lowbnd, upper=uppbnd)
rep <- sdreport(sp);rep
sp$report()

par(mfrow = c(1,2))
for(i in 1:data$NRegion){
  plot(seq(1,33), index[,i], pch = 16, main = i)
  lines(seq(1,33), sp$report()$Ifit[,i])
}

