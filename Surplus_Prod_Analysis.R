#devtools::install_github("dill/beyonce")
#require packages needed
library(rjags)
library(data.table)
library(dplyr)
library(beyonce)
library(mcmcplots)

GOM.adult <- read.csv("C:/Users/w986430/Dropbox/Blue Crab Surplus Production/Code/Surplus Production_GitHub/adult.ioa.jags.csv")
head(GOM.adult)

adult.pos.cor <- filter(GOM.adult, Region.index == 1|Region.index ==3|Region.index ==6|Region.index ==7|Region.index ==8|Region.index ==11|Region.index ==13)

#function to summarize posterior distributions 
source("../post_summ_function.R")

#Year.ind is the years as 1-33
#GOM.adult$Year.ind <- rep(1:length(Years), 4)
Years <- unique(GOM.adult$Year)

#Loop.index is the number of rows, used to loop through all data points
#GOM.adult$index <- 1:nrow(GOM.SP)
Loop.index <- GOM.adult$index
adult.pos.cor$index <- 1:nrow(adult.pos.cor)
Loop.index <- adult.pos.cor$index

#Start.vec is the row of the first year for each state
Start.vec <- GOM.adult$index[which(GOM.adult$Year.Ind==1)]
Start.vec <- adult.pos.cor$index[which(adult.pos.cor$Year.Ind==1)]

#Will loop through all data points except the first year for each state
Loop.index <- Loop.index[-Start.vec]


#Gulf wide catches
year.total.catches <- as.data.frame(GOM.adult %>% group_by(Year)%>% summarise(total.catch.years = sum(Pounds)) %>% arrange(Year))
total.catch.gom <- year.total.catches[,2]
u.year.ind <- 1:length(total.catch.gom)

col. <- c(beyonce_palette(90), beyonce_palette(15)[2:6])
plot(1, type = "n", xlim = range(GOM.adult$Year), ylim = range(GOM.adult$IOA, na.rm = T), main = "IOA by State", xlab = "Year", ylab = "Index of Abundance")
for(i in 1:13){
  points(GOM.adult[which(GOM.adult$Region.index==i),1], GOM.adult[which(GOM.adult$Region.index==i),5], col = col.[i], type = "l", lwd = 2)
}
legend("topright", legend = unique(GOM.adult$Region), col = col.[1:13], lwd = 2, bty = "n")

col.bio <- beyonce_palette(45)  
plot(1, type = "n", xlim = range(GOM.SP$Year), ylim = range(GOM.SP$Total.Catch), main = "Catch by State", ylab = "Catch", xlab = "Year")
for(i in 1:4){
  points(GOM.SP[which(GOM.SP$State.ID==i),3], GOM.SP[which(GOM.SP$State.ID==i),5], col = col.bio[i], type = "l", lwd = 2)
}
legend("topright", legend = unique(GOM.SP$State), col = col.bio, lwd = 2, bty = "n")

par(mfrow = c(3,2), mar = rep(4,4))
for (i in 1:5){
  plot(GOM.SP[which(GOM.SP$State.ID==i),1], GOM.SP[which(GOM.SP$State.ID==i),2], type = "l", ylab = "Catch (10^6)", xlab = "Year", lwd = 2, col = col.[i], main = levels(GOM.SP$State)[i]) # first plot
  par(new = TRUE)
  plot(GOM.SP[which(GOM.SP$State.ID==i),1], GOM.SP[which(GOM.SP$State.ID==i),3],  axes = FALSE, bty = "n", xlab = "", ylab = "", pch =16, col = col.[i])
  axis(side=4, at = pretty(range(GOM.SP[which(GOM.SP$State.ID==i),3])))
  mtext("IOA", side=4, line=2, cex = .75)
  
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
estBetaParams(6E5, 5E6)

data = list(N.obs = length(GOM.adult$Year),
            Year = GOM.adult$Year.Ind, 
            Catch = GOM.adult$Pounds, 
            IOA = GOM.adult$IOA, 
            State = GOM.adult$Region.index,  
            N.State = length(unique(GOM.adult$Region.index)),
            Loop.index = Loop.index,
            Start.vec = Start.vec)


n.chains = 5
n.adapt = 6000
n.update = 50000
n.iter = 10000000
params = c("q", "K", "r", "r.mu", "tau.proc", "index.new", "p", "pvalue.fit", "msy.state", "fmsy", "bmsy")
jm.pooled = jags.model(file = "Surplus_Production_jag_model1.R", data = data, n.adapt = n.adapt, n.chain = n.chains)

update(jm.pooled, n.inter = n.update)
zc.pooled = coda.samples(jm.pooled, variable.names = params, n.iter = n.iter, thin = 5000)
gelman.diag(zc.pooled, multivariate = F)

plot(zc.pooled[,c("r[1]", "r[2]", "r[3]", "r[4]", "r[5]", "r[6]", "r[7]", "r[8]", "r[9]", "r[10]", "r[11]", "r[12]", "r[13]")])
plot(zc.pooled[,c("K[1]", "K[2]", "K[3]", "K[4]", "K[5]", "K[6]", "K[7]", "K[8]", "K[9]", "K[10]", "K[11]", "K[12]", "K[13]")])
par(mfrow=c(3,2))
plot(zc.pooled[,c("q[1]", "q[2]", "q[3]", "q[4]", "q[5]", "q[6]", "q[7]", "q[8]", "q[9]", "q[10]", "q[11]", "q[12]", "q[13]")])


mcmc <- as.data.frame(rbindlist(lapply(zc.pooled, as.data.frame)))

denplot(mcmc, parms = c("K[1]", "K[2]", "K[3]", "K[4]", "K[5]"))
denplot(mcmc, parms = c("r[1]", "r[2]", "r[3]", "r[4]", "r[5]", "r.mu"))
denplot(mcmc, parms = c("q[1]", "q[2]", "q[3]", "q[4]", "q[5]"))
denplot(mcmc, parms = c("tau.obs[1]", "tau.obs[2]", "tau.obs[3]", "tau.obs[4]", "tau.obs[5]"))
denplot(mcmc, parms = c("tau.proc"))
denplot(mcmc, parms = "p.total")


r.est <- post.summ(zc.pooled, "r[");r.est
rmu.est <- post.summ(zc.pooled, "r.mu");rmu.est
k.est <- post.summ(zc.pooled, "K[");k.est
q.est <- post.summ(zc.pooled, "q[");q.est
I.new <- post.summ(zc.pooled, "index.new[")
al.IOA <- I.new[3,c(1:33)]
al25 <- I.new[4,c(1:33)]
al975 <- I.new[5,c(1:33)]
flapa.IOA <- I.new[3,c(34:66)]
flapa25 <- I.new[4,c(34:66)]
flapa975 <- I.new[5,c(34:66)]
flcdk.IOA <- I.new[3,c(67:99)]
flcdk25 <- I.new[4,c(67:99)]
flcdk975 <- I.new[5,c(67:99)]
flchh.IOA <- I.new[3,c(100:132)]
flchh25 <- I.new[4,c(100:132)]
flchh975 <- I.new[5,c(100:132)]
fltmb.IOA <- I.new[3,c(133:165)]
fltmb25 <- I.new[4,c(133:165)]
fltmb975 <- I.new[5,c(133:165)]
labar.IOA <- I.new[3,c(166:198)]
labar25 <- I.new[4,c(166:198)]
labar975 <- I.new[5,c(166:198)]
lacal.IOA <- I.new[3,c(199:231)]
lacal25 <- I.new[4,c(199:231)]
lacal975 <- I.new[5,c(199:231)]
lapon.IOA <- I.new[3,c(232:264)]
lapon25 <- I.new[4,c(232:264)]
lapon975 <- I.new[5,c(232:264)]
later.IOA <- I.new[3,c(265:297)]
later25 <- I.new[4,c(265:297)]
later975 <- I.new[5,c(265:297)]
laver.IOA <- I.new[3,c(298:330)]
laver25 <- I.new[4,c(298:330)]
laver975 <- I.new[5,c(298:330)]
la.IOA <- I.new[3,c(331:363)]
la25 <- I.new[4,c(331:363)]
la975 <- I.new[5,c(331:363)]
ms.IOA <- I.new[3,c(364:396)]
ms25 <- I.new[4,c(364:396)]
ms975 <- I.new[5,c(364:396)]
tx.IOA <- I.new[3,c(397:429)]
tx25 <- I.new[4,c(397:429)]
tx975 <- I.new[5,c(397:429)]


#Observed and posterior predictive IOA plots
par(mfrow = c(2,2))
plot(Years, GOM.adult[which(GOM.adult$Region.index==1),5], main = paste(unique(GOM.adult$Region)[1]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[1])
lines(Years, al.IOA, lwd = 2, col = col.[1])
lines(Years, al25, col = col.[1])
lines(Years, al975, col = col.[1])

plot(Years, GOM.adult[which(GOM.adult$Region.index==2),5], main = paste(unique(GOM.adult$Region)[2]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[2])
lines(Years, flapa.IOA, lwd = 2, col = col.[2])
lines(Years, flapa25, col = col.[2])
lines(Years, flapa975, col = col.[2])

plot(Years, GOM.adult[which(GOM.adult$Region.index==3),5], main = paste(unique(GOM.adult$Region)[3]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[3])
lines(Years, flcdk.IOA, lwd = 2, col = col.[3])
lines(Years, flcdk25, col = col.[3])
lines(Years, flcdk975, col = col.[3])

plot(Years, GOM.adult[which(GOM.adult$Region.index==4),5], main = paste(unique(GOM.adult$Region)[4]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[4])
lines(Years, flchh.IOA, lwd = 2, col = col.[4])
lines(Years, flchh25, col = col.[4])
lines(Years, flchh975, col = col.[4])

plot(Years, GOM.adult[which(GOM.adult$Region.index==5),5], main = paste(unique(GOM.adult$Region)[5]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[5])
lines(Years, fltmb.IOA, lwd = 2, col = col.[5])
lines(Years, fltmb25, col = col.[5])
lines(Years, fltmb975, col = col.[5])

plot(Years, GOM.adult[which(GOM.adult$Region.index==6),5], main = paste(unique(GOM.adult$Region)[6]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[6])
lines(Years, labar.IOA, lwd = 2, col = col.[6])
lines(Years, labar25, col = col.[6])
lines(Years, labar975, col = col.[6])

plot(Years, GOM.adult[which(GOM.adult$Region.index==7),5], main = paste(unique(GOM.adult$Region)[7]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[7])
lines(Years, lacal.IOA, lwd = 2, col = col.[7])
lines(Years, lacal25, col = col.[7])
lines(Years, lacal975, col = col.[7])

plot(Years, GOM.adult[which(GOM.adult$Region.index==8),5], main = paste(unique(GOM.adult$Region)[8]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[8])
lines(Years, lapon.IOA, lwd = 2, col = col.[8])
lines(Years, lapon25, col = col.[8])
lines(Years, lapon975, col = col.[8])

plot(Years, GOM.adult[which(GOM.adult$Region.index==9),5], main = paste(unique(GOM.adult$Region)[9]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[9])
lines(Years, later.IOA, lwd = 2, col = col.[9])
lines(Years, later25, col = col.[9])
lines(Years, later975, col = col.[9])

plot(Years, GOM.adult[which(GOM.adult$Region.index==10),5], main = paste(unique(GOM.adult$Region)[10]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[10])
lines(Years, laver.IOA, lwd = 2, col = col.[10])
lines(Years, laver25, col = col.[10])
lines(Years, laver975, col = col.[10])

plot(Years, GOM.adult[which(GOM.adult$Region.index==11),5], main = paste(unique(GOM.adult$Region)[11]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[11])
lines(Years, la.IOA, lwd = 2, col = col.[11])
lines(Years, la25, col = col.[11])
lines(Years, la975, col = col.[11])

plot(Years, GOM.adult[which(GOM.adult$Region.index==12),5], main = paste(unique(GOM.adult$Region)[12]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[12])
lines(Years, ms.IOA, lwd = 2, col = col.[12])
lines(Years, ms25, col = col.[12])
lines(Years, ms975, col = col.[12])

plot(Years, GOM.adult[which(GOM.adult$Region.index==13),5], main = paste(unique(GOM.adult$Region)[13]), ylab = "IOA", xlab = "Years", pch = 16, col = col.[13])
lines(Years, tx.IOA, lwd = 2, col = col.[13])
lines(Years, tx25, col = col.[13])
lines(Years, tx975, col = col.[13])



p.est <- post.summ(zc.pooled, "p[")
al.pest <- select(as.data.frame(p.est), contains("p[1,"))
al.best <- apply(al.pest, 2, function(x) x*k.est[1])
flapa.pest <- select(as.data.frame(p.est), contains("p[2"))
flapa.best <- apply(flapa.pest, 2, function(x) x*k.est[2])
flcdk.pest <- select(as.data.frame(p.est), contains("p[3"))
flcdk.best <- apply(flcdk.pest, 2, function(x) x*k.est[3])
flchh.pest <- select(as.data.frame(p.est), contains("p[4"))
flchh.best <- apply(flchh.pest, 2, function(x) x*k.est[4])
fltmb.pest <- select(as.data.frame(p.est), contains("p[5"))
fltmb.best <- apply(fltmb.pest, 2, function(x) x*k.est[5])
labar.pest <- select(as.data.frame(p.est), contains("p[6"))
labar.best <- apply(labar.pest, 2, function(x) x*k.est[6])
lacal.pest <- select(as.data.frame(p.est), contains("p[7"))
lacal.best <- apply(lacal.pest, 2, function(x) x*k.est[7])
lapon.pest <- select(as.data.frame(p.est), contains("p[8"))
lapon.best <- apply(lapon.pest, 2, function(x) x*k.est[8])
later.pest <- select(as.data.frame(p.est), contains("p[9"))
later.best <- apply(later.pest, 2, function(x) x*k.est[9])
laver.pest <- select(as.data.frame(p.est), contains("p[10"))
laver.best <- apply(laver.pest, 2, function(x) x*k.est[10])
la.pest <- select(as.data.frame(p.est), contains("p[11"))
la.best <- apply(la.pest, 2, function(x) x*k.est[11])
ms.pest <- select(as.data.frame(p.est), contains("p[12"))
ms.best <- apply(ms.pest, 2, function(x) x*k.est[12])
tx.pest <- select(as.data.frame(p.est), contains("p[13"))
tx.best <- apply(tx.pest, 2, function(x) x*k.est[13])

bmsy <- post.summ(zc.pooled, "bmsy[")
p1 <- select(as.data.frame(p.est), contains(",1]"));p1


#Biomass plots
col.trans <- adjustcolor(col., alpha.f = 0.3)
png(file="biomass.est.png" , height=12,  width=15 , pointsize=28, units = "in", res = 800)
par(mfrow = c(3,2), mar = c(4,4,1,1))
plot(Years, al.best[3,],  type = "l", lwd = 2, col = col.[1], ylim = c(0,max(al.best[5,])), ylab = "Biomass Estimates", main = "AL")
lines(Years, al.best[4,], col = col.[1])
lines(Years, al.best[5,], col = col.[1])
polygon(c(Years,rev(Years)), c(al.best[4,],rev(al.best[5,])),  col = col.trans[1])

plot(Years, flapa.best[3,],  type = "l", lwd = 2, col = col.[2], ylim = c(0,max(flapa.best[5,])), ylab = "Biomass Estimates", main = "FL.apa")
lines(Years, flapa.best[4,], col = col.[2])
lines(Years, flapa.best[5,], col = col.[2])
polygon(c(Years,rev(Years)), c(flapa.best[4,],rev(flapa.best[5,])),  col = col.trans[2])

plot(Years, flcdk.best[3,],  type = "l", lwd = 2, col = col.[3], ylim = c(0,max(flcdk.best[5,])), ylab = "Biomass Estimates", main = "FL.cdk")
lines(Years, flcdk.best[4,], col = col.[3])
lines(Years, flcdk.best[5,], col = col.[3])
polygon(c(Years,rev(Years)), c(flcdk.best[4,],rev(flcdk.best[5,])),  col = col.trans[3])

plot(Years, flchh.best[3,],  type = "l", lwd = 2, col = col.[4], ylim = c(0,max(flchh.best[5,])), ylab = "Biomass Estimates", main = "FL.chh")
lines(Years, flchh.best[4,], col = col.[4])
lines(Years, flchh.best[5,], col = col.[4])
polygon(c(Years,rev(Years)), c(flchh.best[4,],rev(flchh.best[5,])),  col = col.trans[4])

plot(Years, fltmb.best[3,],  type = "l", lwd = 2, col = col.[5], ylim = c(0,max(fltmb.best[5,])), ylab = "Biomass Estimates", main = "FL.tmb")
lines(Years, fltmb.best[4,], col = col.[5])
lines(Years, fltmb.best[5,], col = col.[5])
polygon(c(Years,rev(Years)), c(fltmb.best[4,],rev(fltmb.best[5,])),  col = col.trans[5])

plot(Years, labar.best[3,],  type = "l", lwd = 2, col = col.[6], ylim = c(0,max(labar.best[5,])), ylab = "Biomass Estimates", main = "LA.bar")
lines(Years, labar.best[4,], col = col.[6])
lines(Years, labar.best[5,], col = col.[6])
polygon(c(Years,rev(Years)), c(labar.best[4,],rev(labar.best[5,])),  col = col.trans[6])

plot(Years, lacal.best[3,],  type = "l", lwd = 2, col = col.[7], ylim = c(0,max(lacal.best[5,])), ylab = "Biomass Estimates", main = "LA.cal")
lines(Years, lacal.best[4,], col = col.[7])
lines(Years, lacal.best[5,], col = col.[7])
polygon(c(Years,rev(Years)), c(lacal.best[4,],rev(lacal.best[5,])),  col = col.trans[7])

plot(Years, lapon.best[3,],  type = "l", lwd = 2, col = col.[8], ylim = c(0,max(lapon.best[5,])), ylab = "Biomass Estimates", main = "LA.pon")
lines(Years, lapon.best[4,], col = col.[8])
lines(Years, lapon.best[5,], col = col.[8])
polygon(c(Years,rev(Years)), c(lapon.best[4,],rev(lapon.best[5,])),  col = col.trans[8])

plot(Years, later.best[3,],  type = "l", lwd = 2, col = col.[9], ylim = c(0,max(later.best[5,])), ylab = "Biomass Estimates", main = "LA.ter")
lines(Years, later.best[4,], col = col.[9])
lines(Years, later.best[5,], col = col.[9])
polygon(c(Years,rev(Years)), c(later.best[4,],rev(later.best[5,])),  col = col.trans[9])

plot(Years, laver.best[3,],  type = "l", lwd = 2, col = col.[10], ylim = c(0,max(laver.best[5,])), ylab = "Biomass Estimates", main = "LA.ver")
lines(Years, laver.best[4,], col = col.[10])
lines(Years, laver.best[5,], col = col.[10])
polygon(c(Years,rev(Years)), c(laver.best[4,],rev(laver.best[5,])),  col = col.trans[10])

plot(Years, la.best[3,],  type = "l", lwd = 2, col = col.[11], ylim = c(0,max(la.best[5,])), ylab = "Biomass Estimates", main = "LA")
lines(Years, la.best[4,], col = col.[11])
lines(Years, la.best[5,], col = col.[11])
polygon(c(Years,rev(Years)), c(la.best[4,],rev(la.best[5,])),  col = col.trans[11])

plot(Years, ms.best[3,],  type = "l", lwd = 2, col = col.[12], ylim = c(0,max(ms.best[5,])), ylab = "Biomass Estimates", main = "MS")
lines(Years, ms.best[4,], col = col.[12])
lines(Years, ms.best[5,], col = col.[12])
polygon(c(Years,rev(Years)), c(ms.best[4,],rev(ms.best[5,])),  col = col.trans[12])

plot(Years, tx.best[3,],  type = "l", lwd = 2, col = col.[13], ylim = c(0,max(tx.best[5,])), ylab = "Biomass Estimates", main = "TX")
lines(Years, tx.best[4,], col = col.[13])
lines(Years, tx.best[5,], col = col.[13])
polygon(c(Years,rev(Years)), c(tx.best[4,],rev(tx.best[5,])),  col = col.trans[13])

dev.off()

#Parameter catterpillar plots
caterplot(mcmc, parms = c("q[1]", "q[2]", "q[3]", "q[4]", "q[5]", "q[6]", "q[7]", "q[8]", "q[9]", "q[10]", "q[11]", "q[12]", "q[13]"))
caterplot(mcmc, parms = c("r[1]", "r[2]", "r[3]", "r[4]", "r[5]", "r[6]", "r[7]", "r[8]", "r[9]", "r[10]", "r[11]", "r[12]", "r[13]"))
caterplot(mcmc, parms = c("K[1]", "K[2]", "K[3]", "K[4]", "K[5]", "K[6]", "K[7]", "K[8]", "K[9]", "K[10]", "K[11]", "K[12]", "K[13]"))

#B/BMSY plots
plot(Years, al.best[3,]/bmsy[3,1], type = "l", ylim = c(0,1.75), col = col.[1], lwd = 2, ylab = "Biomass/BMSY")
lines(Years, flapa.best[3,]/bmsy[3,2], col = col.[2], lwd = 2)
lines(Years, flcdk.best[3,]/bmsy[3,3], col = col.[3], lwd = 2)
lines(Years, flchh.best[3,]/bmsy[3,4], col = col.[4], lwd = 2)
lines(Years, fltmb.best[3,]/bmsy[3,5], col = col.[5], lwd = 2)
lines(Years, labar.best[3,]/bmsy[3,6], col = col.[6], lwd = 2)
lines(Years, lacal.best[3,]/bmsy[3,7], col = col.[7], lwd = 2)
lines(Years, lapon.best[3,]/bmsy[3,8], col = col.[8], lwd = 2)
lines(Years, later.best[3,]/bmsy[3,9], col = col.[9], lwd = 2)
lines(Years, laver.best[3,]/bmsy[3,10], col = col.[10], lwd = 2)
lines(Years, la.best[3,]/bmsy[3,11], col = col.[11], lwd = 2)
lines(Years, ms.best[3,]/bmsy[3,12], col = col.[12], lwd = 2)
lines(Years, tx.best[3,]/bmsy[3,13], col = col.[13], lwd = 2)
abline(h=1, lwd = 2)

