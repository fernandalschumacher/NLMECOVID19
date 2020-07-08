####################################################################################################
## SM(S)N Nonlinear Mixed Effects Model for COVID-19 deaths data##
####################################################################################################
## last updated in 2020-07-01
## downloaded from https://github.com/fernandalschumacher/NLMECOVID19

#loading packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(skewlmm) #optional. This package fits the linear version of the model, its documentation could be helpful to understand this code
library(mvtnorm)
library(nlme) #used for initial values

#loading SMSN-NLME functions
source("smsn-nlmm.R")
source("smn-nlmm.R")
source("bootstrap-skew-fct.R")
source("loglik-NL.R")

## creating generalized logistic function and derivative
derivlog3 <- function(x,loga,logb,logC) exp(logC)*exp(loga)*
  exp(-exp(logC)*x)/((exp(logb)+exp(-exp(logC)*x))^2) #function without random effects, used for initial values

logist3<- function(x1,beta1,beta2,beta3,beta4,b1,b2) { 
  exp(beta4+beta1+beta3+b2-exp(beta3+b2)*x1)/
    ((exp(beta2+b1)+exp(-exp(beta3+b2)*x1))^(exp(beta4)+1))
}
der3logist <- deriv( ~   exp(beta4+beta1+beta3+b2-exp(beta3+b2)*x1)/
                       ((exp(beta2+b1)+exp(-exp(beta3+b2)*x1))^(exp(beta4)+1)),
                     c("beta1","beta2","beta3","beta4","b1","b2"), function(x1,beta1,beta2,beta3,beta4,b1,b2){})

###################################################################################################
## selecting countries to be used, names must match Johns Hopkins data
countries <- c("United Kingdom","Italy","US","Mexico","Peru","Brazil","Chile","Colombia","Belgium")

#reading cumulative death data 
Y<-read.csv(url("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"),
            sep=',',stringsAsFactors = T)
Y.1 <- Y %>% filter(Country.Region %in% countries,Province.State=='')
Y.1$Country.Region

# melt dataset 
Y.2 <- melt(Y.1[-c(1,3,4)],id.vars = c("Country.Region"))
Y.2<-droplevels(Y.2)
Y.2$date.char <- sub("X","",Y.2$variable)
Y.2$date<- as.Date.character(sub("X","",Y.2$variable),format="%m.%d.%y")
max(Y.2$date)

# creating new deaths variable
for (country in levels(Y.2$Country.Region)) Y.2$newcases[Y.2$Country.Region==country] = 
  c(0,diff(Y.2$value[Y.2$Country.Region==country]))
ggplot(Y.2,aes(x=date,y=newcases,colour=Country.Region))+geom_line()

###############################################
#using data since first death for each country
###############################################
dat <- filter(Y.2,Y.2$value>=1)
scalecte <- min(tapply(dat$newcases,dat$Country.Region,sd)) #constant used for numerical stability
dat$newcasesT <- dat$newcases/scalecte #transformed variable
#
dat <- dat[order(dat$Country.Region),]
# creating variable day (different for each country)
for (country in levels(dat$Country.Region)) dat$day[dat$Country.Region==country] = 1:sum(dat$Country.Region==country)

ggplot(dat,aes(x=day,y=newcases))+geom_line()+facet_wrap(~Country.Region,scales = "free_y")+
  ylab("deaths") + xlab("days since first death")

# creating dataset to use for prediction
n.ahead = 50
newData <- data.frame(NULL)
for (i in 1:n.ahead) {
  dati <- subset(dat,date==max(date),select = c("Country.Region","date","day"))
  dati$date<-dati$date+i
  dati$day<-dati$day+i
  newData <- rbind(newData,dati)
}
newData <- newData[order(newData$Country.Region),]

#####################################################################################
# fitting
#####################################################################################
#WARNING: nonlinear mixed models are quite sensible to initial values, they might need adjustment


# fitting nlme for getting initial values
lognlme3.1=nlme(newcasesT~derivlog3(day,loga,logb,logC),data=dat,
                fixed=loga+logb+logC~1,random=logb+logC~1|Country.Region,
                start=c(loga=5,logb=-4,logC=-3),
                control = list(msMaxIter=200,maxIter=200,returnObject=T)) 
fixef(lognlme3.1)
beta1=as.numeric(lognlme3.1$coefficients$fixed)
sigmae=lognlme3.1$sigma^2
D1=var(ranef(lognlme3.1))
lambda=c(3,3)*as.numeric(sign(moments::skewness(ranef(lognlme3.1))))

##############################################################################
# fitting the model
##############################################################################
ind=dat$Country.Region;x=matrix(dat$day)
ind_levels <- levels(ind)

############################# normal fit
fitN<- EM.sim.NL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,
                  nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(1)),
                  sigmae=sigmae,
                  D1=D1,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                  precisao=1e-4,max.iter=500,showiter=T,showerroriter=T)
fitN$theta
(indpar <- exp(matrix(fitN$estimates$beta,nrow=length(countries),ncol=4,byrow = T) +
                 cbind(0,fitN$random.effects,0)))
#total deaths expected
indpar[,1]/indpar[,2]^indpar[,4]*scalecte
#
fittedN<- fitN$fitted
fitlogN<- data.frame(value=dat$newcases,day=dat$day,
                     date = dat$date,Country.Region=dat$Country.Region,
                     fitted=fittedN*scalecte)
ggplot(fitlogN,aes(x=day,y=value))+geom_line()+geom_point()+ylab("deaths") +
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~Country.Region,scales = "free_y")

#predicting
predfitN<- predictf.skew.NL(yfit=dat$newcasesT,xfit=matrix(dat$day),
                            indfit=dat$Country.Region,xpred=matrix(newData$day),
                            indpred=newData$Country.Region,
                            distr="sn",nlfder=der3logist,nlf=logist3,
                            theta=c(fitN$theta,0,0),p=4,
                            estbi=fitN$random.effects)
newData$valuePredN <- predfitN$ypred*scalecte
dat$fittedN <- fittedN*scalecte
gN<-ggplot(dat,aes(x=date,y=newcases))+geom_line()+geom_point()+
  geom_line(aes(x=date,y=fittedN,color=I('blue'))) + facet_wrap(~Country.Region,scales = "free_y")+
  geom_line(aes(x=date,y=valuePredN,color=I('blue')),data=newData) + ylab("deaths")
gN

############################# skew-normal fit
fitSN<- EM.Skew.NL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,
                  nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(1)),
                  sigmae=sigmae,
                  D1=D1,lambda=lambda,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                  precisao=1e-3,max.iter=700,showiter=T,showerroriter=T)
fitSN$theta
(indpar <- exp(matrix(fitSN$estimates$beta,nrow=length(countries),ncol=4,byrow = T) +
                 cbind(0,fitSN$random.effects,0)))
#total cases expected
indpar[,1]/indpar[,2]^indpar[,4]*scalecte
#
fittedSN <- fitSN$fitted
fitlogSN<- data.frame(value=dat$newcases,day=dat$day,
                     date = dat$date,Country.Region=dat$Country.Region,
                     fitted=fittedSN*scalecte)
ggplot(fitlogSN,aes(x=day,y=value))+geom_line()+geom_point()+ylab("deaths")+
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~Country.Region,scales = "free_y")

#prediction
predfitSN<- predictf.skew.NL(yfit=dat$newcasesT,xfit=matrix(dat$day),
                            indfit=dat$Country.Region,xpred=matrix(newData$day),
                            indpred=newData$Country.Region,
                            distr="sn",nlfder=der3logist,nlf=logist3,theta=fitSN$theta,p=4,
                            estbi=fitSN$random.effects)
newData$valuePredSN <- predfitSN$ypred*scalecte
dat$fittedSN <- fittedSN*scalecte
gSN<-ggplot(dat,aes(x=date,y=newcases))+geom_line()+geom_point()+
  geom_line(aes(x=date,y=fittedSN,color=I('blue'))) + facet_wrap(~Country.Region,scales = "free_y")+
  geom_line(aes(x=date,y=valuePredSN,color=I('blue')),data=newData) + ylab("deaths") 
gSN

############################# t fit
fitT<- EM.sim.NL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,
                 nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(1)),
                 sigmae=sigmae,
                 D1=D1,
                 distr="st",nu=5,lb=1.1,lu=100,
                 precisao=1e-4,max.iter=700,showiter=T,showerroriter=T)
fitT$theta
(indpar <- exp(matrix(fitT$estimates$beta,nrow=length(countries),ncol=4,byrow = T) +
                 cbind(0,fitT$random.effects,0)))
#total cases expected
indpar[,1]/indpar[,2]^indpar[,4]*scalecte
#
fittedT<- fitT$fitted
fitlogT<- data.frame(value=dat$newcases,day=dat$day,
                  date = dat$date,Country.Region=dat$Country.Region,
                  fitted=fittedT*scalecte)
ggplot(fitlogT,aes(x=day,y=value))+geom_line()+geom_point()+ ylab("deaths")+
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~Country.Region,scales = "free_y")

#prediction
predfitT<- predictf.skew.NL(yfit=dat$newcasesT,xfit=matrix(dat$day),
                            indfit=dat$Country.Region,xpred=matrix(newData$day),
                            indpred=newData$Country.Region,
                            distr="st",nlfder=der3logist,nlf=logist3,
                            theta=c(fitT$theta[-9],0,0,fitT$theta[9]),p=4,
                            estbi=fitT$random.effects)
newData$valuePredT <- predfitT$ypred*scalecte
dat$fittedT <- fittedT*scalecte
gT<-ggplot(dat,aes(x=date,y=newcases))+geom_line()+geom_point()+
  geom_line(aes(x=date,y=fittedT,color=I('blue'))) + facet_wrap(~Country.Region,scales = "free_y")+
  geom_line(aes(x=date,y=valuePredT,color=I('blue')),data=newData) + ylab("deaths")
gT

############################# skew-t fit
fitST <- EM.Skew.NL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,
                   nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(1)),
                   sigmae=sigmae,
                   D1=D1,lambda=lambda,
                   distr="st",nu=5,lb=1.1,lu=100,
                   precisao=1e-3,max.iter=700,showiter=T,showerroriter=T)
fitST$theta
(indpar <- exp(matrix(fitST$estimates$beta,nrow=length(countries),ncol=4,byrow = T) +
                 cbind(0,fitST$random.effects,0)))
#total cases expected
indpar[,1]/indpar[,2]^indpar[,4]*scalecte
#
fittedST <- fitST$fitted
fitlogST<- data.frame(value=dat$newcases,day=dat$day,
                      date = dat$date,Country.Region=dat$Country.Region,
                      fitted=fittedST*scalecte)
ggplot(fitlogST,aes(x=day,y=value))+geom_line()+geom_point()+ylab("deaths")+
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~Country.Region,scales = "free_y")

#prediction
predfitST<- predictf.skew.NL(yfit=dat$newcasesT,xfit=matrix(dat$day),
                             indfit=dat$Country.Region,xpred=matrix(newData$day),
                             indpred=newData$Country.Region,
                             distr="st",nlfder=der3logist,nlf=logist3,theta=fitST$theta,p=4,
                             estbi=fitST$random.effects)
newData$valuePredST <- predfitST$ypred*scalecte
dat$fittedST <- fittedST*scalecte
gST<-ggplot(dat,aes(x=date,y=newcases))+geom_line()+geom_point()+
  geom_line(aes(x=date,y=fittedST,color=I('blue'))) + facet_wrap(~Country.Region,scales = "free_y")+
  geom_line(aes(x=date,y=valuePredST,color=I('blue')),data=newData) + ylab("deaths") 
gST

##############################################
#selection criteria
(l1 <- logveroNL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,nlfder=der3logist,
                 fitObj = fitN,distr="sn")) #normal
(l2 <- logveroNL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,nlfder=der3logist,
                 fitObj = fitSN,distr="sn")) #sn
(l3 <- logveroNL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,nlfder=der3logist,
                 fitObj = fitT,distr="st")) #t
(l4 <- logveroNL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region,nlfder=der3logist,
                 fitObj = fitST,distr="st")) #st

N=nrow(dat)
cbind(c(l1,l2,l3,l4) #loglik
      ,c( #AIC
-2*l1+2*length(fitN$theta),
-2*l2+2*length(fitSN$theta),
-2*l3+2*length(fitT$theta),
-2*l4+2*length(fitST$theta)),
# BIC
c(-2*l1+log(N)*length(fitN$theta),
-2*l2+log(N)*length(fitSN$theta),
-2*l3+log(N)*length(fitT$theta),
-2*l4+log(N)*length(fitST$theta)))


#######################################################################
#bootstrap
#######################################################################
# this function was developed specifically for COVID-19 dataset, it needs adapting for using with other datasets
# it takes a while to run, depending on Mboot and distr
# package parallel must be installed

dat1 <- dat[,c('Country.Region','day','date')]
dat1pred <- newData[,c('Country.Region','day','date')]

######
# t
bootTres<-bootNLMM(fitObj=fitT,dat1=dat1,dat1pred=dat1pred,
                   distr="st",symmetrical=T,Mboot=500)
bootTres$elapsedTime
bootTres$bootfit <- bootTres$bootfit*scalecte

# st
bootSTres<-bootNLMM(fitObj=fitST,dat1=dat1,dat1pred=dat1pred,
                   distr="st",symmetrical=F,Mboot=500)
bootSTres$elapsedTime
bootSTres$bootfit <- bootSTres$bootfit*scalecte


#######################################################
# plotting bootstrap results

fitValues <- fittedT # vector with fitted values
predValues <- predfitT$ypred #vector with predicted values
bootM2res <- bootTres$bootfit
#
sum(is.na(bootM2res[nrow(bootM2res),])) #any sample with numerical or convergence problems?
#
dat1Plus <- rbind(dat1,dat1pred)
dat1Plus$fit <- c(fitValues,predValues)*scalecte
dat1Plus$yy <- c(dat$newcases,rep(NA,nrow(dat1pred))) 
#
bootM2res0 <- bootM2res[,!is.na(bootM2res[nrow(bootM2res),])]
EQM <- apply(bootM2res0[1:nrow(dat),],2,function(x) sum((x-dat$newcases)^2))
#
# trim
trimprop <- .85
bootM2res1 <- bootM2res0[,which(EQM<= quantile(EQM,probs=trimprop))]
ncol(bootM2res1) #number of samples to be used 
#
maxpred <- max(dat1$date)+50 #at most n.ahead
bootM2res2 <- bootM2res1[dat1Plus$date<=maxpred,]
dat1Plus <- dat1Plus[dat1Plus$date<=maxpred,]
#
levelIC <- .05
dat1Plus$CIinf <- apply(bootM2res2,1,quantile,probs=c(levelIC/2),na.rm=T)
dat1Plus$CIsup<-apply(bootM2res2,1,quantile,probs=c(1-levelIC/2),na.rm=T)
#
levels_ctry <- levels(dat$Country.Region)
modaCI<-matrix(ncol=2,nrow=length(levels_ctry))
modaest <- matrix(ncol=1,nrow=length(levels_ctry))
for (i in seq_along(levels_ctry)) {
  dati <- dat1Plus[dat1Plus$Country.Region==levels_ctry[i],]
  modaCI[i,] <- as.character(range(dati$date[dati$CIsup>=max(dati$CIinf)]))
  modaest[i,] <- as.character(dati$date[which.max(dati$fit)])
}
modaCI <- data.frame(modaest = as.Date(modaest),
                     CIinfm = as.Date(modaCI[,1]),
                     CIsupm = as.Date(modaCI[,2]),
                     Country.Region=levels_ctry)
modaCI
#
dattext<-apply(cbind("Peak: ",format(modaCI$modaest,"%b %d")," \nPeak 95% CI: \n(",
                     format(modaCI$CIinfm,"%b %d"),", ",
                     format(modaCI$CIsupm,"%b %d"),")"),1,paste0,collapse = "" )
dat_text <- data.frame(
  label = dattext,
  Country.Region   = modaCI$Country.Region,
  x = min(dat$date)
)

#############################
dat1Plus2 <- merge(dat1Plus,modaCI[,-1],by="Country.Region",all.x = T)
gCI<-ggplot(dat1Plus2,aes(x=date,y=yy))+
  geom_rect(aes(xmin=CIinfm,xmax=CIsupm,ymin=-Inf,ymax=Inf), fill="grey", alpha=0.05)+
  geom_line()+geom_point()+#size=1.5
  geom_line(aes(x=date,y=fit,color=I('blue'))) + facet_wrap(~Country.Region,scales = "free_y")+
  ylab("deaths") +  geom_ribbon(aes(ymin=CIinf,ymax=CIsup), fill="blue", alpha=0.2) 

gCI2<-gCI + geom_text(
  data    = dat_text,
  size    = 2.5,
  mapping = aes(x =x, y = Inf, label = label),
  hjust   = .1,
  vjust   = 1.1
)
gCI2

#############################################################
# estimate total of deaths h-day-ahead, h<=n.ahead
dat1Plus <- rbind(dat1,dat1pred)
dat1Plus$fit <- c(fitValues,predValues)*scalecte
#
maxpred <- max(dat1$date)+30
bootM2res2 <- bootM2res1[dat1Plus$date<=maxpred,]
dat1Plus <- dat1Plus[dat1Plus$date<=maxpred,]
#
ind_levels <- levels(dat$Country.Region)
totaldeath<- matrix(nrow=length(ind_levels),ncol=ncol(bootM2res2))
for (i in 1:ncol(bootM2res2)) {
  totaldeath[,i]<-tapply(bootM2res2[,i],dat1Plus$Country.Region,sum,na.rm=T)
}
esttotal<- round(cbind(tapply(dat1Plus$fit,dat1Plus$Country.Region,sum),
t(apply(totaldeath,1,quantile,probs=c(.025,.975)))),0)
colnames(esttotal) <- c(paste("Total expected by",maxpred),"2.5% CI", "97.5% CI")
esttotal
