#devtools::install_github("dill/beyonce")
#require packages needed
library(rjags)
library(data.table)
library(dplyr)
library(beyonce)
library(mcmcplots)

GOM.juv <- read.csv("C:/Users/w986430/Dropbox/Blue Crab Surplus Production/Code/Surplus Production_GitHub/juv.ioa.jags.csv")
head(GOM.juv)
#function to summarize posterior distributions 
source("../post_summ_function.R")

#Year.ind is the years as 1-33
#GOM.juv$Year.ind <- rep(1:length(Years), 4)
Years <- unique(GOM.juv$Year)

#Loop.index is the number of rows, used to loop through all data points
#GOM.juv$index <- 1:nrow(GOM.SP)
Loop.index <- GOM.juv$index

#Start.vec is the row of the first year for each state
Start.vec <- GOM.juv$index[which(GOM.juv$Year.Ind==1)]

#Will loop through all data points except the first year for each state
Loop.index <- Loop.index[-Start.vec]


#Gulf wide catches
year.total.catches <- as.data.frame(GOM.juv %>% group_by(Year)%>% summarise(total.catch.years = sum(Pounds)) %>% arrange(Year))
total.catch.gom <- year.total.catches[,2]
u.year.ind <- 1:length(total.catch.gom)

col. <- c(beyonce_palette(90), beyonce_palette(15)[2:6])
plot(1, type = "n", xlim = range(GOM.juv$Year), ylim = range(GOM.juv$IOA, na.rm = T), main = "IOA by State", xlab = "Year", ylab = "Index of Abundance")
for(i in 1:13){
  points(GOM.juv[which(GOM.juv$Region.index==i),1], GOM.juv[which(GOM.juv$Region.index==i),5], col = col.[i], type = "l", lwd = 2)
}
legend("topright", legend = unique(GOM.juv$Region), col = col.[1:13], lwd = 2, bty = "n")

data = list(N.obs = length(GOM.juv$Year),
            Year = GOM.juv$Year.Ind, 
            Catch = GOM.juv$Pounds, 
            IOA = GOM.juv$IOA, 
            State = GOM.juv$Region.index,  
            N.State = length(unique(GOM.juv$Region.index)),
            Loop.index = Loop.index,
            Start.vec = Start.vec)


n.chains = 3
n.adapt = 6000
n.update = 50000
n.iter = 100000
params = c("q", "K", "r", "r.mu", "tau.proc", "index.new", "p")
jm.pooled = jags.model(file = "Surplus_Production_juvenile_jag_model.R", data = data, n.adapt = n.adapt, n.chain = n.chains)

update(jm.pooled, n.inter = n.update)
zc.pooled = coda.samples(jm.pooled, variable.names = params, n.iter = n.iter, thin = 5000)
gelman.diag(zc.pooled, multivariate = F)





