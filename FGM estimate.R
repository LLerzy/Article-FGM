


#Density of the copula function FGM.
dfgm=function(u,v,phi){
  c=1+phi*(1-2*u)*(1-2*v)
  return(c)
}

#Derivative of the log-likelihood function.

dflv1=function(phi,u,v,der=1){
  dln=NA;  n=length(u)
  for (i in 1:n) {
    num=(1-2*u[i])*(1-2*v[i])
    den=1+phi*num
    if(der==1){
      dln[i]=(num/den)
    }else if(der==2){
      dln[i]=-(num/den)^2
    }
  }
  sum(dln)
}


####################################################################
## Data simulation
## Vector simulation of the FGM distribution.
###########################

rfgm=function(phi){
  u1=runif(n = 1,0,1) #value of u1
  v2=runif(n = 1,0,1) #value of conditional probability
  A=phi*(2*u1-1)-1
  B=sqrt(1-2*phi*(2*u1-1)+phi^2*(2*u1-1)^2+4*phi*v2*(2*u1-1))
  u2=2*v2/(B-A)
  return(c(u1,u2))
}

###########################
## Generate samples of size n
###########################
rnfgm=function(phi,n){
  datafgm=matrix(data = NA,nrow = 0,ncol = 2,dimnames = list(c(),c("u","v")))
  for (i in 1:n) {
    datafgm=rbind(datafgm,rfgm(phi))
  }
  return(datafgm)
}

####################################################################
## Definition of models to be use in OpenBugs
###########################

# Triangular distribution as prior distribution
sink("modelo1.txt")
cat("model
      {  
      const <-10000
      for(i in 1 : n) {
      z[i] <- 1
      z[i] ~ dbern(p[i])
      p[i] <- ( 1 + phi*(1-2*u[i])*(1-2*v[i]) ) / const
      }
      y <- 0
      y ~ dloglik(rho)
      indicator1 <- step(phi-c); indicator2 <- step(1+phi); indicator3 <- step(1-phi)
      LogP <- indicator2*(1-indicator1)*(phi+1)/(c+1) + indicator1*indicator3*(1-phi)/(1-c)
      phi ~ dflat()
      rho <- log(LogP)
      }",fill=TRUE)
sink()

# Beta distribution as prior distribution.
sink("modelo2.txt")
cat("model
      {  
      for(i in 1 : n) {
      z[i] <- 0
      z[i] ~ dloglik(logL[i])
      logL[i] <- log( 1 + phi*(1-2*u[i])*(1-2*v[i]) )
      }
      theta ~ dbeta(a,b)
      phi <- 1-2*theta
      }",fill=TRUE)
sink()

# Uniform distribution as prior distribution.
sink("modelo3.txt")
cat("model
      {  
      for(i in 1 : n) {
      z[i] <- 0
      z[i] ~ dloglik(logL[i])
      logL[i] <- log( 1 + phi*(1-2*u[i])*(1-2*v[i]) )
      }
      phi ~ dunif(-1,1)
      }",fill=TRUE)
sink()

####################################################################
# Algorithmic to estimate the parameter phi using Bayesian Method
####################################################################

posterior=function(Datas,modelo,a=0,b=0,c=0){
  parametro=c("phi")
  if(modelo=="modelo3.txt"){ # Uniforme distribution
    Data <- list (u=Datas[,1], v=Datas[,2], n=nrow(Datas))
    Inits=NULL
  }else if(modelo=="modelo2.txt"){ # Beta distribution
    Inits=function(){list(phi=rbeta(1,1,1))}
    Data <- list (u=Datas[,1], v=Datas[,2], n=nrow(Datas),a=a, b=b)
  }else if(modelo=="modelo1.txt"){ # Triangular  distribution
    Inits=function(){list(phi=runif(1,0,1))}
    Data <- list (u=Datas[,1], v=Datas[,2], n=nrow(Datas), c=c)
  }
  
  ideb.sim <- bugs(data=Data, inits=Inits, parameters = parametro,
                   model.file=modelo, n.iter=5000, n.burnin = 1000,
                   n.chains=1)
  return(ideb.sim)
}

####################################################################
# Moments estimate
####################################################################

me=function(Data){
  phipearson=cor(Data,method = "pearson")[1,2]*3
  phikendall=cor(Data,method = "kendall")[1,2]*9/2
  phispearman=cor(Data,method = "spearman")[1,2]*3
  return(list(phipearson=phipearson,phikendall=phikendall,phispearman=phispearman))
}

####################################################################
# Moments estimate
####################################################################



####################################################################
# Moments estimate
####################################################################
