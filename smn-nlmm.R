####################################################################################################
## SMN Nonlinear Mixed Effects Model ##
####################################################################################################
## last updated in 2020-07-01
## downloaded from https://github.com/fernandalschumacher/NLMECOVID19

#function for estimating a SMN-NLME model

EM.sim.NL<- function(y,#response vector (N x 1)
                     x,#covariates matrix (N x p1)
                     ind,#factor whose levels indicates the subjects or groups (N x 1) 
                     nlfder,#first derivative of the nonlinear function 
                     nlf,#nonlinear function with args: x1,...,xp1, beta1,...,betap, b1,...,bq
                     distr,#distribution to be used: "sn" (normal), "st" (student's-t), "ss" (slash), or "scn" (contaminated normal)
                     beta1,#p x 1 vector with initial values for beta
                     sigmae,#initial value for sigma2
                     D1,#q x q matrix with initial values for D 
                     nu=NULL,#initial values for nu, when distr != "sn"
                     lb=NULL,#lower bound to be used inside optim for updating nu, when distr != "sn"
                     lu=NULL,#upper bound to be used inside optim for updating nu, when distr != "sn"
                     precisao=1e-4,#tolarance for the convergence criterion
                     max.iter=600,#maximum number of iterations for the EM algorithm
                     showiter=T,#logical. Should the iteration message be printed?
                     showerroriter=F#logical. Should some convergence information be printed?
                     ){
  ti = Sys.time()
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p1<-ncol(x)
  p <- length(beta1)
  q1 <- dim(D1)[2]
  #
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
  #
  bi <- matrix(0,ncol=q1,nrow=m)
  #
  criterio<-10
  count<-0
  llji = 1
  
  loglikVec<- numeric(max.iter)
  while((criterio > precisao)&(count<max.iter)){
    
    count <- count + 1
    # linearization step
    Wtil<-matrix(nrow=N,ncol=p)
    Htil<-matrix(nrow=N,ncol=q1)
    ytil<-numeric(N)
    for (j in 1:nlevels(ind)) {
      jseq = ind==levels(ind)[j]
      y1=y[jseq]
      x1=matrix(x[jseq,  ],ncol=p1)
      ub1 <- matrix(bi[j,],ncol=1)
      nj = length(y1)
      HWmat <- matrix(0,nrow=nj,ncol=(p+q1))
      fxtil <- matrix(0,nrow=nj,ncol=1)
      for(i in 1:nj)
      {
        farg <- as.list(c(x1[i,],beta1,ub1))
        names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep=""))
        formals(nlfder) <- farg
        fx <- nlfder()
        fxtil[i,] <- fx[1]
        HWmat[i,] <- attr(fx,"gradient")
      }
      
      Wtil[jseq,] <- matrix(HWmat[,1:p],nrow=nj,ncol=p)            # W til
      Htil[jseq,] <- matrix(HWmat[,(p+1):(p+q1)],nrow=nj,ncol=q1)    # H til
      ytil[jseq] <- y1 - fxtil + Wtil[jseq,]%*%beta1 + Htil[jseq,]%*%ub1
    }
    
    # expectation step
    res_emj = revert_list(tapply(1:N,ind,emjs,y=ytil, x=Wtil, z=Htil, beta1=beta1, D1=D1,
                                 sigmae=sigmae, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    uj = unlist(res_emj$uj,use.names = F)
    bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
    
    # maximization step
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    D1<-sum4/m
    #
    logvero1<-function(nu){logveros(ytil, Wtil, Htil, ind, beta1, sigmae, D1, distr, nu)}
    #
    if (distr=="sn"){ nu<-NULL} else
    {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
    }
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
    loglikVec[count] <- logveros(ytil, Wtil, Htil, ind, beta1, sigmae, D1, distr, nu)
    if (count>2){
      at<- (loglikVec[count]-loglikVec[count-1])/(loglikVec[count-1]-loglikVec[count-2])
      criterio<-abs((loglikVec[count]-loglikVec[count-1])/(1-at))
      #print(loglik[count])
    }
    if (is.nan(criterio)) criterio=10
    if (all(is.nan(teta))) stop("NaN values")
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") 
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio," - loglik =",loglikVec[count],"\r") #  criterium ",criterio," or ",criterio2,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  ### update bi and ui
  Wtil<-matrix(nrow=N,ncol=p)
  Htil<-matrix(nrow=N,ncol=q1)
  ytil<-numeric(N)
  for (j in 1:nlevels(ind)) {
    jseq = ind==levels(ind)[j]
    y1=y[jseq]
    x1=matrix(x[jseq,  ],ncol=p1)
    ub1 <- matrix(bi[j,],ncol=1)
    nj = length(y1)
    HWmat <- matrix(0,nrow=nj,ncol=(p+q1))
    fxtil <- matrix(0,nrow=nj,ncol=1)
    for(i in 1:nj)
    {
      farg <- as.list(c(x1[i,],beta1,ub1))
      names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep=""))
      formals(nlfder) <- farg
      fx <- nlfder()
      fxtil[i,] <- fx[1]
      HWmat[i,] <- attr(fx,"gradient")
    }
    
    Wtil[jseq,] <- matrix(HWmat[,1:p],nrow=nj,ncol=p)            # W til
    Htil[jseq,] <- matrix(HWmat[,(p+1):(p+q1)],nrow=nj,ncol=q1)    # H til
    ytil[jseq] <- y1 - fxtil + Wtil[jseq,]%*%beta1 + Htil[jseq,]%*%ub1
  }
  
  #
  res_emj = revert_list(tapply(1:N,ind,emjs,y=ytil, x=Wtil, z=Htil, beta1=beta1, D1=D1,
                               sigmae=sigmae, distr=distr,nu=nu))
  uj = unlist(res_emj$uj,use.names = F)
  bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
  ###
  # creating object to return
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,dd,nu)
  if (distr=="sn") names(theta)<-c(paste0("beta",1:p),"sigma2",paste0("Dsqrt",1:length(dd)))
  else names(theta)<- c(paste0("beta",1:p),"sigma2",paste0("Dsqrt",1:length(dd)),paste0("nu",1:length(nu)))
  
  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           D=D1),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  
  obj.out$random.effects<- bi
  
  obj.out$loglik <-loglikVec[count]
  # fitted values
  fittedval<- numeric(length(ind))
  ind_levels <- levels(ind)
  for (i in seq_along(ind_levels)) {
    seqi <- which(ind==ind_levels[i])
    xfiti <- matrix(x[seqi,],ncol=ncol(x))
    ub1 <- matrix(bi[i,],ncol=1)
    #ub1 <- bi[i,]
    for (j in seq_along(seqi)) 
    {
      farg <- as.list(c(xfiti[j,],beta1,ub1))
      names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep=""))
      formals(nlf) <- farg
      fittedval[seqi[j]] <- nlf()
    }
  }
  obj.out$fitted <- fittedval
  
  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}


###############################################################
##### auxiliary functions
###############################################################

################################################################
#Log-likelihood - independent
################################################################
ljnormals <-function(j,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(dmvnorm(y1,med,Psi))
}

ljts <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  log(dtj)
}

ljss <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))
  resp <- integrate(Vectorize(f2),0,1)$value
  log(nu*resp)
}

ljcns <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
         (1-nu[1])*dmvnorm(y1,med,Psi)))
}
logveros = function(y,x,z,ind,beta1,sigmae,D1,distr,nu){ #ind = indicadora de individuo
  m<-n_distinct(ind)
  N<-length(ind)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormals,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljts,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljss,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcns,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  lv
}

##############################################################################
# EM - independent
##############################################################################
emjs = function(jseq, y, x, z, beta1,D1, sigmae,distr,nu) {
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(nj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if  (distr=="sn"){
    uj<-1
  }
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  }
  
  if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  }
  
  if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }
  
  bi<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)
  
  Tbj<-solve(solve(D1)+t(z1)%*%z1/sigmae)
  r<-Tbj%*%t(z1)%*%(y1-x1%*%beta1)/sigmae
  ub<-uj*r
  ub2j<-Tbj+uj*r%*%t(r)
  #
  sum1<-uj*t(x1)%*%x1 #denom beta
  sum2<-(t(x1)%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%z1%*%ub-
    t(ub)%*%t(z1)%*%(y1-x1%*%beta1)+traceM(ub2j%*%t(z1)%*%z1) #soma do sig2
  sum4<-ub2j #soma do Gamma
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,uj=uj)
  obj.out$bi=bi
  return(obj.out)
}
