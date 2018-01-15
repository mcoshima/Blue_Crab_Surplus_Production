# This is a hierarchically structured schaeffer surplus production model that allows estimation of Biomass from distinct subpopulations from distinct time series of relative abundance and catch data. The model assumes that each subpopulation has a unique carrying capacity, catchability, and instrinsic rate of population growth estimated from global priors.

model{
  
  # Parameter priors
  
  # Prior for global intrinsic rate of population growth
  r.mu ~ dunif(0,2)
  # Prior for global carrying capacity
  K.mu ~ dunif(0, 5.3E10)
  # Prior for catchability
  q.mu <- 0.00005
  
  
  
  # Hierarchical variance priors
  r0_sig ~ dunif(0,10)
  r0_sig2 <- r0_sig*r0_sig
  r0_tau <- 1/r0_sig2
  k0_sig ~ dunif(0,1E10)
  k0_sig2 <- k0_sig*k0_sig
  k0_tau <- 1/k0_sig2
  q0_sig ~ dunif(0,10)
  q0_sig2 <- q0_sig*q0_sig
  q0_tau <- 1/q0_sig2
  
  # State-specific r, k, q and variance priors
  for(i in 1:N.State){
    K[i] ~ dnorm(K.mu, k0_tau)
    r[i] ~ dnorm(r.mu, r0_tau)T(0,2)
    q[i] ~ dnorm(q.mu, q0_tau)T(0,)
  }
  
  sigma.proc ~ dunif(0,100000)
  var.proc <- sigma.proc * sigma.proc
  tau.proc <- 1/var.proc
  
  sigma.obs ~ dunif(0,100)
  var.obs <- sigma.obs * sigma.obs
  tau.obs <- 1/var.obs
  
  #Initial state-specific biomasses
  p.mu[State[Start.vec[1]],Year[Start.vec[1]]] <- 0.75
  p.mu[State[Start.vec[2]],Year[Start.vec[2]]] <- 0.75
  p.mu[State[Start.vec[3]],Year[Start.vec[3]]] <- 0.75
  p.mu[State[Start.vec[4]],Year[Start.vec[4]]] <- 0.75
  p.mu[State[Start.vec[5]],Year[Start.vec[5]]] <- 0.75
  p.mu[State[Start.vec[6]],Year[Start.vec[6]]] <- 0.75
  p.mu[State[Start.vec[7]],Year[Start.vec[7]]] <- 0.75
  p.mu[State[Start.vec[8]],Year[Start.vec[8]]] <- 0.75
  p.mu[State[Start.vec[9]],Year[Start.vec[9]]] <- 0.75
  p.mu[State[Start.vec[10]],Year[Start.vec[10]]] <- 0.75
  p.mu[State[Start.vec[11]],Year[Start.vec[11]]] <- 0.75
  p.mu[State[Start.vec[12]],Year[Start.vec[12]]] <- 0.75
  p.mu[State[Start.vec[13]],Year[Start.vec[13]]] <- 0.75
  
  p[State[Start.vec[1]],Year[Start.vec[1]]] ~ dgamma(5,10)
  p[State[Start.vec[2]],Year[Start.vec[2]]] ~ dgamma(5,10)
  p[State[Start.vec[3]],Year[Start.vec[3]]] ~ dgamma(5,10)
  p[State[Start.vec[4]],Year[Start.vec[4]]] ~ dgamma(5,10)
  p[State[Start.vec[5]],Year[Start.vec[5]]] ~ dgamma(5,10)
  p[State[Start.vec[6]],Year[Start.vec[6]]] ~ dgamma(5,10)
  p[State[Start.vec[7]],Year[Start.vec[7]]] ~ dgamma(5,10)
  p[State[Start.vec[8]],Year[Start.vec[8]]] ~ dgamma(5,10)
  p[State[Start.vec[9]],Year[Start.vec[9]]] ~ dgamma(5,10)
  p[State[Start.vec[10]],Year[Start.vec[10]]] ~ dgamma(5,10)
  p[State[Start.vec[11]],Year[Start.vec[11]]] ~ dgamma(5,10)
  p[State[Start.vec[12]],Year[Start.vec[12]]] ~ dgamma(5,10)
  p[State[Start.vec[13]],Year[Start.vec[13]]] ~ dgamma(5,10)
  
  
  # Forward projection SBM
  # Loop.index is the row number
  # Calculate mean biomass for each state[t] and year[t] from the biomass from state[t] and the previous year, the state specific rate of population growth (r), and the state specific carrying capacity (K) 
  # Biomass is estimated from a normal distribution with mean of state and year specific biomass and error tau.b
  for (t in Loop.index){
    
    p.mu[State[t], Year[t]] <- (p[State[t],Year[t-1]] + p[State[t], Year[t-1]] * r[State[t]] * (1-p[State[t], Year[t-1]])-Catch[Year[t-1]]/K[State[t]])
    
    p[State[t], Year[t]] ~ dnorm(p.mu[State[t], Year[t]], tau.proc)
    
  }
  
  # Likelihood
  # Estimate Index of Abundance from a normal distribution, state and year specific biomass* global catachability coefficient (q)  
  for (i in 1:N.obs){
    I[i] <- (p[State[i], Year[i]]*q[State[i]]*K[State[i]])
    IOA[i] ~ dnorm(I[i], tau.obs)
    
    #Posterior predictions
    index.new[i] ~ dnorm(I[i], tau.obs)
    
  }
  
  #Derived quantities
  
  # for(i in 1:N.obs){
  #   biomass[State[i], Year[i]]<- p[State[i], Year[i]]*K[State[i]]
  # }
  # 
  for(i in 1:N.State){
    msy.state[i] <- r[i]*K[i]/4
    fmsy[i] <- r[i]/2
    bmsy[i] <- K[i]/2
  }
  ###ADD 5 back in!
  #K.gulf <- K[1] + K[2] + K[3] + K[4] + K[5]
  
  # for (a in U.year){
  #   
  #   #p.gom[a] <- p[1,a] + p[2,a] +p[3,a] + p[4,a] + p[5,a]
  #   
  # }
  Bp-value
  for(i in 1:N.obs){
    sq[i] <- (IOA[i] - I[i])^2
    sq.new[i] <- (index.new[i] - I[i])^2
  }
  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  pvalue.fit <- step(fit.new - fit)
}