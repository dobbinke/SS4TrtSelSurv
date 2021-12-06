rm(list=ls())

normalqfun = function(p) {
  qnorm(p,mean=normalmean,sd=normalsd);
}



beta1fun = function(k1, k2, FUN = normalqfun) {
  
  if (k1 == k2) { myanswer=0; }
  else {
    mynum = log(log(k1) / log(k2));
    myden = FUN(0.25)-FUN(0.75);
    myanswer = mynum/myden;
  }
  myanswer;
}

lambdafun = function(k1,k2,t0,FUN = normalqfun) {
  if (k1==k2) { myanswer = (-log(k1))/t0; }
  else {
    mynum = log(k2/k1);
    myden = t0*(exp(beta1fun(k1, k2, FUN)*FUN(0.25)) - exp(beta1fun(k1, k2, FUN)*FUN(0.75)));;
    myanswer = mynum/myden;
  }
  myanswer;
}

beta3fun = function(k1,k2,k3,k4,t0,FUN = normalqfun) {
  if (k1==k2 && k3==k4) {
    myanswer = 0;
  }
  else if (k1==k2 && k3 != k4) {
    myanswer = (k3t-k4t)/(FUN(0.25)-FUN(0.75));
  }
  else {
    mynum = log(log(k2)/log(k1)) - log(log(k4)/log(k3));
    myden = FUN(0.25) - FUN(0.75) ;
    myanswer = mynum/myden;
  }
  myanswer;
}

beta2fun = function(k1, k2, k3, k4, t0, FUN = normalqfun) {
  k1t = log(-log(k1));
  k2t = log(-log(k2));
  k3t = log(-log(k3));
  k4t = log(-log(k4));
  
  if (k1 == k2 && k3 == k4) {
    myanswer = log(log(k3) / log(k1));
  }
  else if (k1 == k2 && k3 != k4) {
    myanswer = k3t - k1t - FUN(0.25) * ((k4t - k3t) / (FUN(0.75) - FUN(0.25)));
  }
  else {
    prelognum = (-1)*log(k3*k4);
    betasum = beta1fun(k1, k2, FUN) + beta3fun(k1, k2, k3, k4, t0, FUN);
    prelogden = t0*lambdafun(k1, k2, t0, FUN)*(exp(betasum*FUN(0.25)) + exp(betasum*FUN(0.75)));
    myanswer = log(prelognum / prelogden);
  }
  myanswer;
}

getlambdasbetas <- function(k1,k2,k3,k4,t0,FUN){
  # function returns a vector of lambda,beta1,beta2,beta3
  lambda = lambdafun(k1,k2,t0,FUN = normalqfun);
  beta1 = beta1fun(k1, k2, FUN = normalqfun);
  beta2 = beta2fun(k1, k2, k3, k4, t0, FUN = normalqfun);
  beta3 = beta3fun(k1, k2, k3, k4, t0, FUN = normalqfun);
  c(lambda,beta1,beta2,beta3);
}


generatedata <- function(lambda,beta1,beta2,beta3,working.n,t0,mydesign) {
  myy = normalmean + normalsd * rnorm(working.n);
  mytrt = rbinom(working.n,1,0.5); # Assumes complete randomization
  if (length(which(mytrt==0))<5) stop("Too few in the control set.");
  if (length(which(mytrt==1))<5) stop("Too few in the treatment set.");
  hazard = lambda * exp(beta1*myy+beta2*mytrt+beta3*myy*mytrt);
  mytimes = rexp(working.n,rate=hazard);
  # Check the number of events in all 4 groups is over 5
  data.frame(y=myy,trt=mytrt,survtime = mytimes);
}



getwidthemp <- function(working.n,lambda,beta1,beta2,beta3,t0,mydesign,NumberOfBootstraps){
  # function estimates the width using MCRUNS simulations
  thewidths = rep(NA,mcruns);
  library(TreatmentSelection);
  library(survival)
  
  for(i in 1:mcruns) {
    thisdata = generatedata(lambda,beta1,beta2,beta3,working.n,t0,mydesign);
    Binary4janes = as.numeric(thisdata$survtime > t0);
    ctrlgrp = which(thisdata[,2]==0);
    trtgrp = which(thisdata[,2]==1);
    nset1 = length(which(thisdata[ctrlgrp,3]<t0));
    nset2 = length(which(thisdata[ctrlgrp,3]>t0));
    nset3 = length(which(thisdata[trtgrp,3]<t0));
    nset4 = length(which(thisdata[trtgrp,3]>t0));
    if (nset1 < 5) stop("Fewer than 5 events before t0 in the control group.")
    if (nset2 < 5) stop("Fewer than 5 events after t0 in the control group.")
    if (nset3 < 5) stop("Fewer than 5 events before t0 in the treatment group.")
    if (nset4 < 5) stop("Fewer than 5 events after t0 in the treatment group.")
    
    
    datafraforjanes = data.frame(event=Binary4janes,trt=thisdata$trt,Y1=thisdata$y);
    mytrtsel = trtsel(event ~ Y1*trt,treatment.name="trt",data=datafraforjanes,
                      study.design="RCT",family=binomial("logit"),default.trt="trt all")
    tmp = evaluate(mytrtsel,bootstraps=NumberOfBootstraps);
    thewidths[i] = tmp$conf.intervals[2,7]-tmp$conf.intervals[1,7];
  }
  c(mean(thewidths),sqrt(var(thewidths)/length(thewidths)));
}


sampsizefcn = function(targetwidth,
  k1,k2,k3,k4,t0,
  Ds=c(100,150,200,250,300,350,400),
  mcruns=100,mydesign,ExpPropOfDeaths,NumberOfBootstraps=200) {
  thisbeta1 = beta1fun(k1,k2,normalqfun);
  thisbeta2 = beta2fun(k1,k2,k3,k4,t0,normalqfun);
  thisbeta3 = beta3fun(k1,k2,k3,k4,t0,normalqfun);
  thislambda = lambdafun(k1,k2,t0,normalqfun);
  thist0 = t0;
  mcruns = mcruns;
  mydesign = mydesign;
  Ds = Ds;
  ExpPropOfDeaths = ExpPropOfDeaths;
  NumberOfBootstraps = NumberOfBootstraps;
  
  library(doParallel);
  #mycluster = makeCluster(1); # comment this out when not running pureMCcomp
  mycluster = makeCluster(length(Ds)); #Comment this out just when running pureMCcomp
  registerDoParallel(mycluster)
  getDoParWorkers()
  library(foreach)
  
  myres <- foreach(i = 1:length(Ds),
                   .combine=rbind,
                   .export=c("getwidthemp","generatedata","mcruns",
                             "normalmean","normalsd","NumberOfBootstraps")
                   ) %dopar% {
    library(TreatmentSelection);
    theseres = getwidthemp(working.n=Ds[i],thislambda,thisbeta1,thisbeta2,thisbeta3,thist0,
      mydesign,NumberOfBootstraps);
    return(c(Ds[i],theseres))
  }
  print(myres);
  stopCluster(mycluster);
  
  mywidths = myres[,2];
  mywidthsinvsq = mywidths^(-2);
  myDs = myres[,1]
  mylm = lm(myDs ~ mywidthsinvsq);
  targetinvsq = 1/(targetwidth*targetwidth);
  plot(mywidthsinvsq,myDs,xlab="Inverse width squared",ylab="Sample Size",main="Point estimates with regression")
  abline(mylm$coefficients[1],mylm$coefficients[2]);
  mypred <- predict(mylm,newdata=list(mywidthsinvsq=targetinvsq),interval="confidence")
  mydf <- data.frame(samplesize=mypred[1]/ExpPropOfDeaths,
                     lowerbound = mypred[2]/ExpPropOfDeaths,
                     upperbound = mypred[3]/ExpPropOfDeaths)
  myret <- list(mydf,mylm,myDs,mywidthsinvsq);
  if ( mypred[1] < min(Ds) )
  {
    print("Extrapolation warning: Estimated sample size below range of simulated sample sizes.")
  }
  if ( mypred[1] > max(Ds) )
  {
    print("Extrapolation warning: Estimated sample size above range of simulated sample sizes.")
  }
  print(paste("Number of deaths/events required: ", ceiling(mypred[1])))
  print(paste("Number of patients required:", ceiling(mydf$samplesize)))
  return(myret)
}


Calculate.sample.size <- function(targetwidth,
  k1,k2,k3,k4,t0,
  mydesign,ExpPropOfDeaths,NumberOfBootstraps=200) {

   approx.sample.size.list <- sampsizefcn(targetwidth,
     k1,k2,k3,k4,t0,
     Ds= 200*(1:10),
     mcruns=20,mydesign,ExpPropOfDeaths,NumberOfBootstraps=NumberOfBootstraps) 

   approx.sample.size.df = approx.sample.size.list[[1]]
   approx.sample.size = as.numeric(approx.sample.size.df$samplesize)
   print(approx.sample.size)
   return(approx.sample.size.list)

#   final.sample.size.df <- sampsizefcn(targetwidth,
#     k1,k2,k3,k4,t0,
#     Ds= c(approx.sample.size/2,approx.sample.size,approx.sample.size*2),
#     mcruns=100,mydesign,ExpPropOfDeaths,NumberOfBootstraps=NumberOfBootstraps) 
# return(final.sample.size.df$samp.size)
   
   
   
}





