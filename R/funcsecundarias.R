#############################
#funcoes utilizadas no pacote
#############################


#malha de pontos simulados
Gerarmalha <- function(nn, a = 0, b = 1) {
  x <- array(NA, dim = c(nn, nn))
  y <- array(NA, dim = c(nn, nn))
  malha <- seq(a, b, length = nn)

  for (i in 1:nn)
  {
    for (j in 1:nn)
    {
      x[j, i] <- malha[i]
      y[j, i] <- malha[j]
    }
  }

  xg <- as.matrix(as.vector(x))

  yg <- as.matrix(as.vector(y))


  res <- cbind(xg, yg)
  res
}

##################################
#função cria a matriz de covariancia
gSigma <- function(b, v, def) {
  n <- nrow(def)
  R <- exp(-b * (as.matrix(dist(def))))
  mat <- v * R
  mat
}


#######################
#gerando dados da normal multivariada

normalmulti<- function(x, PSI, b, v,pontos_simulados){

  #valores necessarios para a normal multivariada
  #vetor de médias com as covariaveis
  media<- x%*%PSI

  #matriz de covariancia
  SS<-gSigma(b,v,pontos_simulados)

  #gerando dados de uma normal multivariada
  dados<- MASS::mvrnorm(1, media, SS)

  return(dados)

}

#######################
#gerando dados da normal multivariada

normalmultiMQ<- function(medi,b,v,pontos_simulados){

  #matriz de covariancia
  SS<-gSigma(b,v,pontos_simulados)

  #gerando dados de uma normal multivariada
  dados<- MASS::mvrnorm(1, log(medi), SS)

  return(dados)

}

##############################

amostrarW=function(WW,MM,QQ,XXw,XX,ff,PPw,bb,vv,S){
  WW<-as.matrix(WW)
  PPw<-matrix(PPw,ncol = 1)


  #trocar para length (dimensão de W)
  n=nrow(WW)


  #Gerar um novo W (W proposto)
  WWprop=as.matrix(MASS::mvrnorm(1,WW,ff*diag(1,n)))

  #matriz de variancia e covariancia
  SSig=gSigma(bb,vv,S)

  postWW=likelihoodhawkes(XX,exp(WW),exp(MM),exp(QQ))-0.5*t(WW-XXw%*%PPw)%*%solve(SSig)%*%(WW-XXw%*%PPw)

  postWWprop=likelihoodhawkes(XX,exp(WWprop),exp(MM),exp(QQ))-0.5*t(WWprop-XXw%*%PPw)%*%solve(SSig)%*%(WWprop-XXw%*%PPw)

  prob=min(exp((postWWprop)-(postWW)),1)

  u=runif(1,0,1)

  if(u<prob){

    Wprox=WWprop
    rejei=1


  }else{

    Wprox=WW
    rejei=0
  }

  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}


amostrarM=function(WW,MM,QQ,XX,ff,bb,vv,S,Mpriori){
  MM<-as.matrix(MM)


  #trocar para length (dimensão de W)
  n=nrow(MM)


  #Gerar um novo W (W proposto)
  MMprop=as.matrix(MASS::mvrnorm(1,MM,ff*diag(1,n)))

  if(sum(ifelse(MMprop<as.matrix(QQ),0,1))>=1){
    res<-list(MM,0)
    return(res)

  } else{

  }
  #matriz de variancia e covariancia
  SSig=gSigma(bb,vv,S)

  postMM=likelihoodhawkes(XX,exp(WW),exp(MM),exp(QQ))-0.5*t(MM-Mpriori)%*%solve(SSig)%*%(MM-Mpriori)

  postMMprop=likelihoodhawkes(XX,exp(WW),exp(MMprop),exp(QQ))-0.5*t(MMprop-Mpriori)%*%solve(SSig)%*%(MMprop-Mpriori)

  prob=min(exp((postMMprop)-(postMM)),1)

  u=runif(1,0,1)

  if(u<prob){

    Mprox=MMprop
    rejei=1


  }else{

    Mprox=MM
    rejei=0
  }

  res=as.matrix(Mprox)
  res=list(Mprox,rejei)
  res
}



amostrarQ=function(WW,MM,QQ,XX,ff,b_beta,v_beta,S,Qpriori){
  QQ<-as.matrix(QQ)



  #trocar para length (dimensão de W)
  n=nrow(QQ)


  #Gerar um novo W (W proposto)
  QQprop=as.matrix(MASS::mvrnorm(1,QQ,ff*diag(1,n)))

  if(sum(ifelse(QQprop>as.matrix(MM),0,1))>=1){
    res<-list(QQ,0)
    return(res)

  } else{

  }
  #matriz de variancia e covariancia
  SSig=gSigma(bb,vv,S)

  postQQ=likelihoodhawkes(XX,exp(WW),exp(MM),exp(QQ))-0.5*t(QQ-Qpriori)%*%solve(SSig)%*%(QQ-Qpriori)

  postQQprop=likelihoodhawkes(XX,exp(WW),exp(MM),exp(QQprop))-0.5*t(QQprop-Qpriori)%*%solve(SSig)%*%(QQprop-Qpriori)

  prob=min(exp((postQQprop)-(postQQ)),1)

  u=runif(1,0,1)

  if(u<prob){

    Qprox=QQprop
    rejei=1


  }else{

    Qprox=QQ
    rejei=0
  }

  res=as.matrix(Qprox)
  res=list(Qprox,rejei)
  res
}




amostrarb=function(W,v_lambda,b_lambda,S,ab,bb,Xw,PSI_lambda,u1){
  W<-as.matrix(W)
  PSI_lambda<-matrix(PSI_lambda,ncol = 1)
  bprop=rgamma(1,shape=b_lambda*u1, rate = u1)

  SSigprop=gSigma(bprop,v_lambda,S)

  if((det(SSigprop)==0)|(bprop< 0.005)){
    return(list(b_lambda,0))
  }

  SSig=gSigma(b_lambda,v_lambda,S)
  SSigprop=gSigma(bprop,v_lambda,S)


  logp=-0.5*t(W-Xw%*%PSI_lambda)%*%solve(SSig)%*%(W-Xw%*%PSI_lambda)-0.5*log(det(SSig))+(ab-1)*log(b_lambda)-bb*b_lambda

  logpprop=-0.5*t(W-Xw%*%PSI_lambda)%*%solve(SSigprop)%*%(W-Xw%*%PSI_lambda)-0.5*log(det(SSigprop))+(ab-1)*log(bprop)-bb*bprop

  logprob=logpprop+log(dgamma(b_lambda,shape=bprop*u1,rate=u1))-(logp+log(dgamma(bprop,shape=b_lambda*u1,rate=u1)))
  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=bprop
    rejei=1

  }else{

    bprox=b_lambda
    rejei=0;

  }
  res=list(bprox,rejei)
  res
}


amostrarbMQ=function(MQ,v,b,S,ab,bb,MQpriori,u1){
  MQ<-as.matrix(MQ)
  bprop=rgamma(1,shape=b*u1, rate = u1)

  SSigprop=gSigma(bprop,v,S)

  if((det(SSigprop)==0)|(bprop< 0.005)){
    return(list(b,0))
  }

  SSig=gSigma(b,v,S)
  SSigprop=gSigma(bprop,v,S)


  logp=-0.5*t(MQ-MQpriori)%*%solve(SSig)%*%(MQ-MQpriori)-0.5*log(det(SSig))+(ab-1)*log(b)-bb*b

  logpprop=-0.5*t(MQ-MQpriori)%*%solve(SSigprop)%*%(MQ-MQpriori)-0.5*log(det(SSigprop))+(ab-1)*log(bprop)-bb*bprop

  logprob=logpprop+log(dgamma(b,shape=bprop*u1,rate=u1))-(logp+log(dgamma(bprop,shape=b*u1,rate=u1)))
  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=bprop
    rejei=1

  }else{

    bprox=b
    rejei=0;

  }
  res=list(bprox,rejei)
  res
}





gCorr<-function(b,def){
  n=nrow(def)
  R=exp(-b*(as.matrix(dist(def))))
  mat=R
  mat

}


########################################
# sintonizador do b (b_lambda, b_alpha, b_beta)
sintonizar=function(bar,taxa,tau,mat,i){

  mat=as.matrix(mat)

  mater=(1/50)*sum(mat[(i-49):i,1])

  if(mater>=taxa){
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)-delta
    temp5=exp(temp4)
    return(temp5)
  }else{
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)+delta
    temp5=exp(temp4)
    return(temp5)
  }



}

