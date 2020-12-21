####################################################################################################
## Bootstrap method for SM(S)N Nonlinear Mixed Effects Model ##
####################################################################################################
## last updated in 2020-07-01
## downloaded from https://github.com/fernandalschumacher/NLMECOVID19

# bootstrap method for obtaining interval information regarding 
# SMN and SMSN nonlinear mixed effects model
# package parallel must be installed 

# this code was developed to fit SMSN-LME model with generalized logistic nonlinear 
#function (here called logist3 and der3logist for its derivative),
# it needs to be adapted for applying to other nonlinear functions and datasets

########
bootNLMM <- function(fitObj,#object returned from EM.Skew.NL (if symmetrical=F) or EM.sim.NL (if symmetrical=T)
                     dat1,#data.frame with columns c('Country.Region','day','date'), used for fitting
                     dat1pred,#data.frame with columns c('Country.Region','day','date'), used for predicting
                     distr, #"sn" or "st"
                     symmetrical, #logical. Should be symmetrical distribution be used?
                     Mboot=500 #number of samples to be generated and estimated
                     ) {
  dat1Plus <- rbind(dat1,dat1pred)
  sigmafit <- sqrt(fitObj$estimates$sigma2)
  D1fit <- fitObj$estimates$D
  nufit <- fitObj$estimates$nu
  fitValues <- fitObj$fitted
  uivec <- fitObj$uhat
  if (!symmetrical) lambdai <- rep(3,nrow(D1fit))*sign(fitObj$estimates$lambda)
  #
  
  library(parallel)
  ncores<- detectCores()-1
  cl <- makeCluster(getOption("cl.cores", ncores))
  clusterExport(cl=cl,varlist = c("der3logist","EM.Skew.NL","EM.sim.NL","rgamma","logist3",
                                  "predictf.skew.NL","matrix.sqrt","revert_list","emjs",
                                  "logveros","ljnormals","ljts","ljt","dmvt",
                                  "emj","traceM","logvero","ljnormal","dmvnorm",
                                  "Dmatrix"))
  ti<-Sys.time()
  set.seed(99987)
  if (distr=='sn') {
    if (symmetrical==F) {
      bootM2res<- revert_list(parSapply(cl,1:Mboot,bootSN,mutC=fitValues,sigmai=sigmafit,
                                        D1i =D1fit,lambdai=lambdai ,dat1=dat1,dat1pred=dat1pred,simplify = F))
    } else{
      bootM2res<- revert_list(parSapply(cl,1:Mboot,bootN,mutC=fitValues,sigmai=sigmafit,
                                        D1i =D1fit ,dat1=dat1,dat1pred=dat1pred,simplify = F))
      }
  } else{
    if (symmetrical==F) {
      bootM2res<- revert_list(parSapply(cl,1:Mboot,bootSTui,mutC=fitValues,sigmai=sigmafit,
                                        D1i =D1fit,lambdai=lambdai,uivec=uivec,dat1=dat1,dat1pred=dat1pred,simplify = F))
    } else{
      bootM2res<- revert_list(parSapply(cl,1:Mboot,bootTui,mutC=fitValues,sigmai=sigmafit,
                                          D1i =D1fit,uivec=uivec,dat1=dat1,dat1pred=dat1pred,simplify = F))
    }
  }
  
  tf<-Sys.time()
  stopCluster(cl)
  #
  obj.out = list()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="mins"))
  obj.out$bootfit = as.matrix(bind_cols(bootM2res[[1]]))
  obj.out$thetabi = as.matrix(bind_cols(bootM2res[[2]]))
  obj.out
}


###############################################################
##### auxiliary functions
###############################################################

bootSN <- function(applyx=NULL,mutC,sigmai,D1i,lambdai,dat1,dat1pred) {
  library(dplyr)
  dat1$y <- mutC + rnorm(nrow(dat1),mean=0,sd=sigmai)
  xmat <- matrix(dat1$day)
  ftest<-seq(1,3,by=.1)
  for (fi in ftest) {#
    skewboot=try(EM.Skew.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                            nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(fi)),
                            sigmae=sigmai^2,
                            D1=D1i,lambda=lambdai,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                            #distr="st",nu=10,lb=1.1,lu=100,
                            precisao=1e-4,max.iter=500,showiter=F,showerroriter=F),silent = F)
    if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
  }
  if (class(skewboot)[1]=="try-error") {
    for (ai in seq(4,7,by=.1)) {#
      skewboot=try(EM.Skew.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                              nlfder=der3logist,nlf=logist3,beta1=c(ai,-4,-3,log(2.5)),
                              sigmae=sigmai^2,
                              D1=D1i,lambda=lambdai,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                              #distr="st",nu=10,lb=1.1,lu=100,
                              precisao=1e-4,max.iter=500,showiter=F,showerroriter=F),silent = F)
      if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
    }
  }
  if (class(skewboot)[1]=="try-error") {
    return(list(fitted=rep(NA,nrow(dat1)+nrow(dat1pred)),theta=rep(NA,10+nlevels(dat1$Country.Region)*2)))
  }
  else {
    ind_levels <- levels(dat1$Country.Region)
    fittedboot <- skewboot$fitted
    predboot<- try(predictf.skew.NL(yfit=dat1$y,xfit=xmat,
                                indfit=dat1$Country.Region,xpred=matrix(dat1pred$day),
                                indpred=dat1pred$Country.Region,
                                distr="sn",nlfder=der3logist,nlf=logist3,
                                theta=skewboot$theta,p=4,
                                estbi=skewboot$random.effects),silent = T)
    if (class(predboot)[1]=="try-error") {
      ypred <- rep(NA,nrow(dat1pred))
    }
    else ypred<-predboot$ypred
    #calc fitted and pred
    outl <- list()
    outl$fitted <- c(fittedboot,ypred)
    outl$theta <- c(skewboot$theta,skewboot$random.effects[,1],skewboot$random.effects[,2])
    return(outl)
  } 
}
#
bootN <- function(applyx=NULL,mutC,sigmai,D1i,dat1,dat1pred) {
  library(dplyr)
  dat1$y <- mutC + rnorm(nrow(dat1),mean=0,sd=sigmai)
  xmat <- matrix(dat1$day)
  ftest<-seq(1,3,by=.1)
  for (fi in ftest) {#
    #print(fi)
    skewboot=try(
      EM.sim.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                            nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(fi)),
                            sigmae=sigmai^2,
                            D1=D1i,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                            precisao=1e-4,max.iter=500,showiter=F,showerroriter=F)
                 ,silent = T)
    if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
  }
  if (class(skewboot)[1]=="try-error") {
    for (ai in seq(4,7,by=.1)) {#
      skewboot=try(EM.sim.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                              nlfder=der3logist,nlf=logist3,beta1=c(ai,-4,-3,log(2.5)),
                              sigmae=sigmai^2,
                              D1=D1i,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                              precisao=1e-4,max.iter=500,showiter=F,showerroriter=F),silent = T)
      if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
    }
  }
  if (class(skewboot)[1]=="try-error") {
    return(list(fitted=rep(NA,nrow(dat1)+nrow(dat1pred)),theta=rep(NA,8+nlevels(dat1$Country.Region)*2)))
  }
  else {
    ind_levels <- levels(dat1$Country.Region)
    fittedboot <- skewboot$fitted
    predboot<- try(predictf.skew.NL(yfit=dat1$y,xfit=xmat,
                                indfit=dat1$Country.Region,xpred=matrix(dat1pred$day),
                                indpred=dat1pred$Country.Region,
                                distr="sn",nlfder=der3logist,nlf=logist3,
                                theta=c(skewboot$theta,0,0),p=4,
                                estbi=skewboot$random.effects),silent = T)
    if (class(predboot)[1]=="try-error") {
      ypred <- rep(NA,nrow(dat1pred))
    }
    else ypred<-predboot$ypred
    #calc fitted and pred
    outl <- list()
    outl$fitted <- c(fittedboot,ypred)
    outl$theta <- c(skewboot$theta,skewboot$random.effects[,1],skewboot$random.effects[,2])
    return(outl)
  } 
}
#
bootSTui <- function(applyx=NULL,mutC,sigmai,D1i,lambdai,uivec,dat1,dat1pred) {
  library(dplyr)
  erroi <- numeric(length=nrow(dat1))
  for (i in 1:nlevels(dat1$Country.Region)) {
    ct = levels(dat1$Country.Region)[i]
    ui <- uivec[i] #rgamma(1,nui/2,nui/2)
    erroi[dat1$Country.Region==ct]<- as.numeric(mvtnorm::rmvnorm(1,sigma=ui^(-1)*sigmai^2*diag(sum(dat1$Country.Region==ct))))
  }
  dat1$y <- mutC + erroi
  xmat <- matrix(dat1$day)
  ftest<-seq(1,3,by=.2)
  for (fi in ftest) {#
    skewboot=try(EM.Skew.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                            nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(fi)),
                            sigmae=sigmai^2,
                            D1=D1i,lambda=lambdai,#distr="sn",nu=NULL,lb=NULL,lu=NULL,
                            distr="st",nu=5,lb=1.1,lu=100,
                            precisao=1e-2,max.iter=600,showiter=F,showerroriter=F),silent = T)
    if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
  }
  if (class(skewboot)[1]=="try-error") {
    for (ai in seq(4,6,by=.2)) {#
      skewboot=try(EM.Skew.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                              nlfder=der3logist,nlf=logist3,beta1=c(ai,-4,-3,log(2.5)),
                              sigmae=sigmai^2,
                              D1=D1i,lambda=lambdai,#distr="sn",nu=NULL,lb=NULL,lu=NULL,
                              distr="st",nu=5,lb=1.1,lu=100,
                              precisao=1e-2,max.iter=600,showiter=F,showerroriter=F),silent = T)
      if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
    }
  }
  if (class(skewboot)[1]=="try-error") {
    return(list(fitted=rep(NA,nrow(dat1)+nrow(dat1pred)),theta=rep(NA,11+nlevels(dat1$Country.Region)*2)))
  }
  else {
    if (skewboot$error>=10) return(list(fitted=rep(NA,nrow(dat1)+nrow(dat1pred)),
                                        theta=rep(NA,11+nlevels(dat1$Country.Region)*2)))
    ind_levels <- levels(dat1$Country.Region)
    fittedboot <-skewboot$fitted
    predboot<- try(predictf.skew.NL(yfit=dat1$y,xfit=xmat,
                                    indfit=dat1$Country.Region,xpred=matrix(dat1pred$day),
                                    indpred=dat1pred$Country.Region,
                                    distr="st",nlfder=der3logist,nlf=logist3,
                                    theta=skewboot$theta,p=4,
                                    estbi=skewboot$random.effects),silent = T)
    
    if (class(predboot)[1]=="try-error") {
      ypred <- rep(NA,nrow(dat1pred))
    } else ypred<-predboot$ypred
    #calc fitted and pred
    outl <- list()
    outl$fitted <- c(fittedboot,ypred)
    outl$theta <- c(skewboot$theta,skewboot$random.effects[,1],skewboot$random.effects[,2])
    return(outl)
  }
}
#
bootTui <- function(applyx=NULL,mutC,sigmai,D1i,uivec,dat1,dat1pred) {
  library(dplyr)
  erroi <- numeric(length=nrow(dat1))
  for (i in 1:nlevels(dat1$Country.Region)) {
    ct = levels(dat1$Country.Region)[i]
    ui <- uivec[i] 
    erroi[dat1$Country.Region==ct]<- as.numeric(mvtnorm::rmvnorm(1,sigma=ui^(-1)*sigmai^2*diag(sum(dat1$Country.Region==ct))))
  }
  dat1$y <- mutC + erroi
  xmat <- matrix(dat1$day)
  ftest<-seq(1,3,by=.1)
  for (fi in ftest) {#
    skewboot=try(
      EM.sim.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                           nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(fi)),
                           sigmae=sigmai^2,
                           D1=D1i,distr="st",nu=5,lb=1.1,lu=100,#nu=NULL,lb=NULL,lu=NULL,
                           precisao=1e-4,max.iter=600,showiter=F,showerroriter=F)
      ,silent = T)
    if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
  }
  if (class(skewboot)[1]=="try-error") {
    for (ai in seq(4,7,by=.1)) {#
      skewboot=try(EM.sim.NL(y=dat1$y,x=xmat,ind=dat1$Country.Region,
                             nlfder=der3logist,nlf=logist3,beta1=c(ai,-4,-3,log(2.5)),
                             sigmae=sigmai^2,
                             D1=D1i,distr="st",nu=5,lb=1.1,lu=100,#nu=NULL,lb=NULL,lu=NULL,
                             precisao=1e-4,max.iter=600,showiter=F,showerroriter=F),silent = T)
      if (class(skewboot)[1]!="try-error") if (skewboot$error<10) break
    }
  }
  if (class(skewboot)[1]=="try-error") {
    return(list(fitted=rep(NA,nrow(dat1)+nrow(dat1pred)),theta=rep(NA,9+nlevels(dat1$Country.Region)*2)))
  }
  else {
    if (skewboot$error>=10)  return(list(fitted=rep(NA,nrow(dat1)+nrow(dat1pred)),theta=rep(NA,9+nlevels(dat1$Country.Region)*2)))
    ind_levels <- levels(dat1$Country.Region)
    fittedboot <- skewboot$fitted
    predboot<- try(predictf.skew.NL(yfit=dat1$y,xfit=xmat,
                                indfit=dat1$Country.Region,xpred=matrix(dat1pred$day),
                                indpred=dat1pred$Country.Region,
                                distr="st",nlfder=der3logist,nlf=logist3,
                                theta=c(skewboot$theta[-9],0,0,skewboot$theta[9]),p=4,
                                estbi=skewboot$random.effects),silent = T)
    if (class(predboot)[1]=="try-error") {
      ypred <- rep(NA,nrow(dat1pred))
    }
    else ypred<-predboot$ypred
        #calc fitted and pred
    outl <- list()
    outl$fitted <- c(fittedboot,ypred)
    outl$theta <- c(skewboot$theta,skewboot$random.effects[,1],skewboot$random.effects[,2])
    return(outl)
  } 
}

