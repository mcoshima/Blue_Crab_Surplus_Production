# This is a hierarchically structured schaeffer surplus production model that allows estimation of Biomass from distinct subpopulations from distinct time series of relative abundance and catch data. The model assumes that each subpopulation has a unique carrying capacity, catchability, and instrinsic rate of population growth estimated from global priors.

model{
  
  # Parameter priors
  
  r ~ dunif(0, 2)
  # Prior for carrying capacity
  K ~ dunif(53000, 5.3E11)
  # Prior for catchability
  q.mu ~ dunif(1.13E-8, 5.4E-8)
  
  
  # Hierarchical variance priors
  q0_sig ~ dunif(0,1)
  q0_sig2 <- q0_sig*q0_sig
  q0_tau <- 1/q0_sig2
  
  # State-specific r, k, q 
  for(i in 1:N.State){
    q[i] ~ dnorm(q.mu, q0_tau)
  }
  
  sigma.proc ~ dunif(0,100000)
  var.proc <- sigma.proc * sigma.proc
  tau.proc <- 1/var.proc
  
for(i in 1:N.State){
  sigma.obs[i] ~ dunif(0,100)
  var.obs[i] <- sigma.obs[i] * sigma.obs[i]
  tau.obs[i] <- 1/var.obs[i]
  
}

  #Initial state-specific biomasses
  
  p.mu[State[Start.vec[1]],Year[Start.vec[1]]] <- 0.75
  p[State[Start.vec[1]],Year[Start.vec[1]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[2]],Year[Start.vec[2]]] <- 0.75
  p[State[Start.vec[2]],Year[Start.vec[2]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[3]],Year[Start.vec[3]]] <- 0.75
  p[State[Start.vec[3]],Year[Start.vec[3]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[4]],Year[Start.vec[4]]] <- 0.75
  p[State[Start.vec[4]],Year[Start.vec[4]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[5]],Year[Start.vec[5]]] <- 0.75
  p[State[Start.vec[5]],Year[Start.vec[5]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[6]],Year[Start.vec[6]]] <- 0.75
  p[State[Start.vec[6]],Year[Start.vec[6]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[7]],Year[Start.vec[7]]] <- 0.75
  p[State[Start.vec[7]],Year[Start.vec[7]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[8]],Year[Start.vec[8]]] <- 0.75
  p[State[Start.vec[8]],Year[Start.vec[8]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[9]],Year[Start.vec[9]]] <- 0.75
  p[State[Start.vec[9]],Year[Start.vec[9]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[10]],Year[Start.vec[10]]] <- 0.75
  p[State[Start.vec[10]],Year[Start.vec[10]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[11]],Year[Start.vec[11]]] <- 0.75
  p[State[Start.vec[11]],Year[Start.vec[11]]] ~ dbeta(1,1)
  p.mu[State[Start.vec[12]],Year[Start.vec[12]]] <- 0.75
  p[State[Start.vec[12]],Year[Start.vec[12]]] ~ dbeta(1,1)


  # Forward projection SBM
  # Loop.index is the row number
  # Calculate mean biomass for each state[t] and year[t] from the biomass from state[t] and the previous year, the state specific rate of population growth (r), and the state specific carrying capacity (K) 
  # Biomass is estimated from a normal distribution with mean of state and year specific biomass and error tau.b
  for (t in Loop.index){
    
    p.mu[State[t], Year[t]] <- max((p[State[t],Year[t-1]] + p[State[t], Year[t-1]] * r * (1-p[State[t], Year[t-1]])-Catch[Year[t-1]]/K),0)
    
    p[State[t], Year[t]] ~ dnorm(p.mu[State[t], Year[t]], tau.proc)
    
  }
  
  # Likelihood
  # Estimate Index of Abundance from a normal distribution, state and year specific biomass* global catachability coefficient (q)  
  for (i in 1:N.obs){
    I[i] <- (p[State[i], Year[i]]*q[State[i]]*K)
    IOA[i] ~ dnorm(I[i], tau.obs[State[i]])
    
    #Posterior predictions
    index.new[i] ~ dnorm(I[i], tau.obs[State[i]])
    
  }
  
  #Derived quantities
  
  for(i in 1:N.obs){
    biomass[State[i], Year[i]]<- p.mu[State[i], Year[i]]*K
  }
  
  for(i in 1:N.Year){
    total.biomass[Year[i]] <- biomass[State[1], Year[i]] + biomass[State[2], Year[i]] + biomass[State[3], Year[i]] + biomass[State[4], Year[i]] + biomass[State[5], Year[i]] + biomass[State[6], Year[i]] + biomass[State[7], Year[i]] + biomass[State[8], Year[i]] + biomass[State[9], Year[i]] + biomass[State[10], Year[i]] + biomass[State[11], Year[i]] + biomass[State[12], Year[i]]
    
  }
  
  k.gulf <- K*12
  
  msy <- (r*k.gulf)/4
  fmsy <- r/2
  bmsy <- k.gulf/2

  
  #Bp-value
 
  #mean
  mean.data <- mean(IOA[])
  mean.sim <- mean(index.new[])
  bp.mean <- step(mean.sim - mean.data)
  
  #variance
  var.data <- sd(IOA[])^2
  var.sim <- sd(index.new[])^2
  bp.var <- step(var.sim - var.data)
  
  #cv
  cv.data <- sd(IOA[])/mean.data
  cv.sim <- sd(index.new[])/mean.sim
  bp.cv <- step(cv.sim - cv.data)
  
   for(i in 1:N.obs){
    sq[i] <- (IOA[i] - I[i])^2
    sq.new[i] <- (index.new[i] - I[i])^2
  }
  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  pvalue.fit <- step(fit.new - fit)
  
  
}