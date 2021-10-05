####################################################################
# Library
####################################################################
library(R2OpenBUGS)
library(copula)

####################################################################
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
# First cycle to find estimates, confidence intervals and credibility regions.
# Here phi is dependence, n is size of sample and confidence is to the intervals
####################################################################
results1=function(phi,n,nboot=200,confidence,a,b,c){
  Datas=rnfgm(phi,n) # Generate of sample
  ####################
  # Obtaining estimates for the parameter
  ####################
  # Using log-likelihood function maximum
  ELM=safeUroot(f=function(ph){mapply(function(ph){
    dflv1(ph,Datas[,1],Datas[,2])},ph)},c(-1,1))$root
  # Using moments
  EM=me(Datas)
  
  # Selection of samples that have estimates in the FGM parameter space
  while((ELM < -1.0 | ELM > 1.0)
        | (EM$phikendall< -1.0 | EM$phikendall > 1.0)
        | (EM$phispearman < -1.0 | EM$phispearman > 1.0)){
    Datas=rnfgm(phi,n)
    ELM=safeUroot(f=function(ph){mapply(function(ph){
      dflv1(ph,Datas[,1],Datas[,2])},ph)},c(-1,1))$root
    EM=me(Datas)
  }
  
  # Estimates using Bayes approach
  
  ResTriangular=posterior(Datas,"modelo1.txt", c=c)
  ResBeta=posterior(Datas,"modelo2.txt",a=a, b=b)
  ResUnif=posterior(Datas,"modelo3.txt")
  
  ####
  # Confidence intervals
  ####
  alpha=1-confidence
  #Obtaining resamples to use Bootstrap
  vphitau=c();vphispe=c();vphilm=c();i=1;ie1=0
  while(i<=nboot) {
    cho=sample(1:n,n,replace = T) #resample of index
    rdatan=Datas[cho,] #selection of sample items
    #estimates
    #philm=Bisection(-1,1,rdatan) #ELM
    philm=safeUroot(f=function(ph){mapply(function(ph){
      dflv1(ph,rdatan[,1],rdatan[,2])},ph)},c(-1,1))$root
    if(philm>=-1 & philm<=1){
      vphilm[i]=philm
      vphim=me(rdatan)
      if(vphim$phikendall>=-1 & vphim$phikendall<=1 & 
         vphim$phispearman>=-1 & vphim$phispearman<=1){
        vphitau=cbind(vphitau,vphim$phikendall)
        vphispe=cbind(vphispe,vphim$phispearman)
        i=i+1
      }}else{ie1=ie1+1}
  }
  IntervalsLM=quantile(vphilm,probs=c(alpha/2,1-alpha/2))#Quantile to estimates using ML
  IntervalsTau=quantile(vphitau,probs=c(alpha/2,1-alpha/2))#Quantile to estimates using moment
  IntervalsSpe=quantile(vphispe,probs=c(alpha/2,1-alpha/2))#Quantile to estimates using moment
  IntervalsBT=c(ResTriangular$summary[1,3],ResTriangular$summary[1,7])
  IntervalsBB=c(ResBeta$summary[1,3],ResBeta$summary[1,7])
  IntervalsBU=c(ResUnif$summary[1,3],ResUnif$summary[1,7])
  
  return(list(ELM=ELM,EM=EM,EBT=ResTriangular$mean$phi,EBB=ResBeta$mean$phi,EBU=ResUnif$mean$phi,
              ILM=IntervalsLM,ITau=IntervalsTau,ISpe=IntervalsSpe,IBT=IntervalsBT,
              IBB=IntervalsBB, IBU=IntervalsBU, vphilm=vphilm,
              vphitau=vphitau,vphispe=vphispe,Error1=ie1))
}

####################################################################
# Second cycle, in N repetition find squared mean error, minimum, maximum, 
# average length, coverage percentage
####################################################################
# Here phi is dependence, n is the sample size, confidence is to the intervals,
# N is the number of times that the first cycle is repeated

results2=function(phi,n,nboot=200, N, confidence,a,b,c){
  PhiML=NA; i=0; PhiTau=NA; PhiSpe=NA
  PhiBT=NA; PhiBB=NA; PhiBU=NA
  IntervalLM=matrix(data=NA,nrow=0,ncol=2,dimnames = list(c(),c("L","R")))
  IntervalTau=IntervalLM; IntervalSpe=IntervalLM
  IntervalBT=IntervalLM;IntervalBB=IntervalLM;IntervalBU=IntervalLM
  while (i<=N) {
    results=results1(phi,n,nboot,confidence,a,b,c)
    #Estimates
    PhiML[i]=results$ELM; PhiTau[i]=results$EM$phikendall; PhiSpe[i]=results$EM$phispearman
    PhiBT[i]=results$EBT; PhiBB[i]=results$EBB; PhiBU[i]=results$EBU
    #Confidence intervals
    IntervalLM=rbind(IntervalLM,results$ILM)
    IntervalTau=rbind(IntervalTau,results$ITau)
    IntervalSpe=rbind(IntervalSpe,results$ISpe)
    IntervalBT=rbind(IntervalBT,results$IBT)
    IntervalBB=rbind(IntervalBB,results$IBB)
    IntervalBU=rbind(IntervalBU,results$IBU)
    i=i+1}
  #Descriptive
  DesML=c(min(PhiML),mean(PhiML),sd(PhiML),max(PhiML))
  DesMTau=c(min(PhiTau),mean(PhiTau),sd(PhiTau),max(PhiTau))
  DesMSpe=c(min(PhiSpe),mean(PhiSpe),sd(PhiSpe),max(PhiSpe))
  DesBT=c(min(PhiBT),mean(PhiBT),sd(PhiBT),max(PhiBT))
  DesBB=c(min(PhiBB),mean(PhiBB),sd(PhiBB),max(PhiBB))
  DesBU=c(min(PhiBU),mean(PhiBU),sd(PhiBU),max(PhiBU))
  #Sesgo
  Bias=c(phi-DesML[2],phi-DesMTau[2],phi-DesMSpe[2],phi-DesBT[2],phi-DesBB[2],phi-DesBU[2])
  
  #Average length 
  meanlenght=function(Interval){
    lenint=sum(Interval[,2]-Interval[,1])/nrow(Interval)
    return(lenint=lenint)
  }
  
  lenghtILM=meanlenght(IntervalLM)
  lenghtITau=meanlenght(IntervalTau)
  lenghtISpe=meanlenght(IntervalSpe)
  lenghtIBT=meanlenght(IntervalBT)
  lenghtIBB=meanlenght(IntervalBB)
  lenghtIBU=meanlenght(IntervalBU)
  
  #Coverage probability
  covprob=function(Interval){
    indicator=function(i){
      if(Interval[i,1]<=phi & Interval[i,2]>=phi){
        return(1)
      }else{return(0)}
    }
    resind=NA
    for(i in 1:nrow(Interval)){
      resind[i]=indicator(i)
    }
    covint=sum(resind)/nrow(Interval)
    
    return(covint=covint)
  }
  covprobILM=covprob(IntervalLM)
  covprobITau=covprob(IntervalTau)
  covprobISpe=covprob(IntervalSpe)
  covprobIBT=covprob(IntervalBT)
  covprobIBB=covprob(IntervalBB)
  covprobIBU=covprob(IntervalBU)
  
  Descritive=matrix(data=c(c(DesML,Bias[1],lenghtILM,covprobILM),
                           c(DesMTau,Bias[2],lenghtITau,covprobITau),
                           c(DesMSpe,Bias[3],lenghtISpe,covprobISpe),
                           c(DesBT,Bias[4],lenghtIBT,covprobIBT),
                           c(DesBB,Bias[5],lenghtIBB,covprobIBB),
                           c(DesBU,Bias[6],lenghtIBU,covprobIBU)),nrow=7,ncol=6,
                    dimnames = list(list("Min","Mean","MSE","Max","Bias","Length","Coverage"),
                                    list("LM","Tau","Spe","Triangular","Beta","Unif")))
  return(Descritive)
}

####################################################################
# Running and saving results
####################################################################
size=c(10,seq(50,1000,50))
name=as.character(size)

for (i in 1 : 21){
  example=NA
  example=results2(phi=0.2,n=size[i],nboot=1000, N=1000, confidence=0.95, a=3071.6, b=4607.4, c=0.2)
  write.xlsx(x=example, file="Resultados02.xlsx",sheetName =name[i],col.names = TRUE,
             row.names = TRUE,append = TRUE)
  print(i)
}

for (i in 18 : 21){
  example=NA
  example=results2(phi=0.5,n=size[i],nboot=1000, N=1000, confidence=0.95, a=11.7, b=32.5, c=0.5)
  write.xlsx(x=example, file="Resultados05.xlsx",sheetName =name[i],col.names = TRUE,
             row.names = TRUE,append = TRUE)
  print(i)
}

for (i in 1 : 21){
  example=NA
  example=results2(phi=0.9,n=size[i],nboot=1000, N=1000, confidence=0.95, a=18.6875, b=280.3125, c=0.9)
  write.xlsx(x=example, file="Resultados09.xlsx",sheetName =name[i],col.names = TRUE,
             row.names = TRUE,append = TRUE)
  print(i)
}

################################################################################
## Graph construction
################################################################################

# Function to construct graphs

Graph=function(m,omega){
  nameGra=c("Min","Mean","MSE","Max","Bias","Length","Coverage")
  MA=matrix(data = NA,nrow=0,ncol = 3,dimnames = list(list(),list("Size",nameGra[m],"Method")))
  namemethod=1:6
  size=c(10,seq(50,1000,50))
  for (i in 1:21) {
    if(omega==0.2){
      A=read.xlsx(file = "Resultados02.xlsx",sheetIndex = i)
    }else if(omega==0.5){
      A=read.xlsx(file = "Resultados05.xlsx",sheetIndex = i)
    }else if(omega==0.9){
      A=read.xlsx(file = "Resultados09.xlsx",sheetIndex = i)
    }
    for (j in 2:7) {
      MA=rbind(MA,c(size[i],A[m,j],namemethod[(j-1)]))
    }
  }
  MA=as.data.frame(MA)
  MA$Method=factor(MA$Method, labels = c("ML","Tau-Kendall","Spearman","Triangular","Beta","Uniform"))
  return(MA)
}

# Plot for Bias
Bias=Graph(5,omega=0.2)
Bias02=ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of 0.2")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")
Bias=Graph(5,omega=0.5)
Bias05=ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of 0.5")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")
Bias=Graph(5,omega=0.9)
Bias09=ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of 0.9")+
  scale_x_continuous(breaks=c(10,seq(50,1000,50)))+xlab("Sample size")
grid.arrange(Bias02, Bias05, Bias09, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                           c(3, 3)))

# Plot for MSE
MSE=Graph(3,omega=0.2)
Mse02=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of 0.2")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")
MSE=Graph(3,omega=0.5)
Mse05=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of 0.5")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")
MSE=Graph(3,omega=0.5)
Mse09=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of 0.5")+
  scale_x_continuous(breaks=c(10,seq(50,1000,50)))+xlab("Sample size")
grid.arrange(Mse02, Mse05, Mse09, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                        c(3, 3)))

# Plot for Average Length
Length=Graph(6,omega=0.2)
leng02=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of 0.2")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=Graph(6,omega=0.5)
leng05=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of 0.5")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=Graph(6,omega=0.9)
leng09=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of 0.9")+
  scale_x_continuous(breaks=c(10,seq(50,1000,50)))+ylab("Average length")+xlab("Sample size")

grid.arrange(leng02, leng05, leng09, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                           c(3, 3)))

# Plot for Probability Coverage
Coverage=Graph(7,omega=0.2)
Coverage02=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of 0.2")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=Graph(7,omega=0.5)
Coverage05=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of 0.5")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=Graph(7,omega=0.9)
Coverage09=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of 0.9")+
  scale_x_continuous(breaks=c(10,seq(50,1000,50)))+ylab("Probability")+xlab("Sample size")

grid.arrange(Coverage02, Coverage02, Coverage09, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                       c(3, 3)))






