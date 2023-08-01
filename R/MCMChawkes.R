MCMChawkes<-function(point,xdata,xcov,iter=100000,bar=90000,pul=2,AMpri=matrix(rep(1,nrow(point)),ncol=1),BMpri=100*diag(1,nrow(point)),AQpri=matrix(rep(1,nrow(point)),ncol=1),BQpri=100*diag(1,nrow(point)),AW=as.matrix(rep(0,ncol(xcov))),BW=100*diag(1,ncol(xcov)),
                     avW=1,bvW=1,avM=1,bvM=1,avQ=1,bvQ=1,abW=.01,bbW=.01,abM=.01,bbM=.01,abQ=.01,bbQ=.01){
  x<-xdata
  Xw<-xcov
  S=point
  n=nrow(S)

  if(iter<=bar){
    print("the number of iterations must be greater than the burn in")
  }
  else{

  }
  if(bar<=pul){
    print("the number of burn in must be greater than the jump")
  }
  else{

  }

  #Parametros necessários com chutes iniciais
  W=rep(1,n)
  M=rep(1,n)
  Q=rep(1,n)
  fw=.5
  PSI_lambda=rep(1,sqrt(n))
  b_lambda=1
  v_lambda=1
  fm=.5
  b_alpha=1
  v_alpha=1
  fq=.5
  b_beta=1
  v_beta=1
  uW=1

  aceitaW<-NULL
  matW<-NULL
  aceitaM<-NULL
  matMeanM<-NULL
  matM<-NULL
  matMeanQ<-NULL
  aceitaQ<-NULL
  matQ<-NULL
  psi_W<-NULL
  psi_M<-NULL
  psi_Q<-NULL
  VetV_lambda<-VetV_alpha<-VetV_beta<-NULL
  aceitaWB<-NULL
  vetWB<-NULL
  aceitaMB<-NULL
  vetMB<-NULL
  aceitaQB<-NULL
  vetQB<-NULL

  #varios MCMC
  for(j in 1:iter){

    if(j<=bar){

      #Para W
      tempW=amostrarW(W,M,Q,Xw,x,fw,PSI_lambda,b_lambda,v_lambda,Xw[,-1])
      W<-as.numeric(tempW[[1]])
      aceitaW<-c(aceitaW,tempW[[2]])

      #Mean M
      varMpriori<-solve(solve(gSigma(b_alpha,v_alpha,S))+solve(BMpri))
      meanMpriori<-varMpriori%*%(solve(gSigma(b_alpha,v_alpha,S))%*%M+solve(BMpri)%*%AMpri)
      MeanM<-mvrnorm(1,meanMpriori,varMpriori)

      #Para M
      tempM=amostrarM(W,M,Q,x,fm,b_alpha,v_alpha,S,MeanM)
      M<-as.numeric(tempM[[1]])
      aceitaM<-c(tempM[[2]],aceitaM)

      #Mean Q
      varQpriori<-solve(solve(gSigma(b_beta,v_beta,S))+solve(BQpri))
      meanQpriori<-varQpriori%*%(solve(gSigma(b_beta,v_beta,S))%*%Q+solve(BQpri)%*%AQpri)
      MeanQ<-mvrnorm(1,meanQpriori,varQpriori)

      #Para Q
      tempQ=amostrarQ(W,M,Q,x,fq,b_beta,v_beta,S,MeanQ)
      Q<-as.numeric(tempQ[[1]])
      aceitaQ<-c(tempQ[[2]],aceitaQ)

      #Para PSI_lambda
      varPsiW=solve(solve(BW)+t(Xw)%*%solve(gSigma(b_lambda,v_lambda,S))%*%Xw)
      medPsiW=(t(AW)%*%solve(BW)+t(W)%*%solve(gSigma(b_lambda,v_lambda,S))%*%Xw)%*%varPsiW
      PSI_lambda=as.matrix(mvrnorm(1,medPsiW,varPsiW))

      #Para b_lambda
      tempBW=amostrarb(W,v_lambda,b_lambda,S,abW,bbW,Xw,PSI_lambda,uW)
      b_lambda<-as.numeric(tempBW[[1]])
      aceitaWB<-c(aceitaWB,tempBW[[2]])

      #Para b_alpha
      tempBM=amostrarbMQ(M,v_alpha,b_alpha,S,abM,bbM,MeanM,uM)
      b_alpha<-as.numeric(tempBM[[1]])
      aceitaMB<-c(aceitaMB,tempBM[[2]])
      #Para b_beta
      tempQB=amostrarbMQ(Q,v_beta,b_beta,S,abQ,bbQ,MeanQ,uQ)
      b_beta<-as.numeric(tempQB[[1]])
      aceitaQB<-c(aceitaQB,tempQB[[2]])

      # #Para v_lambda
      RRW=gCorr(b_lambda,S)
      aW=(n/2)+avW
      bW=0.5*t(W-Xw%*%as.matrix(PSI_lambda))%*%solve(RRW)%*%(W-Xw%*%as.matrix(PSI_lambda))+bvW
      v_lambda=1/rgamma(1,shape=aW, rate = bW)
      #
      #Para v_alpha
      RRM=gCorr(b_alpha,S)
      aM=(n/2)+avM
      bM=0.5*t(M-MeanM)%*%solve(RRM)%*%(M-MeanM)+bvM
      v_alpha=1/rgamma(1,shape=aM, rate = bM)

      #Para v_beta
      RRQ=gCorr(b_beta,S)
      aQ=(n/2)+avQ
      bQ=0.5*t(Q-MeanQ)%*%solve(RRQ)%*%(Q-MeanQ)+bvQ
      v_beta=1/rgamma(1,shape=aQ, rate = bQ)


      if((j%%50)==0){

        uW=sintonizar(bar,0.44,uW,aceitaWB,j)
        uM=sintonizar(bar,0.44,uM,aceitaMB,j)
        uQ=sintonizar(bar,0.44,uQ,aceitaQB,j)

      }else{

      }

    }else{

      if(j==(bar+1)){

        #rm(MDefT,MkappaT,MPhiT)
        rm(aceitaW)
        aceitaW=NULL
        rm(aceitaM)
        aceitaM=NULL
        rm(aceitaQ)
        aceitaQ=NULL
      }else{

      }

      if((j%%pul)==0){

        #Para W
        tempW=amostrarW(W,M,Q,Xw,x,fw,PSI_lambda,b_lambda,v_lambda,Xw[,-1])
        W<-as.numeric(tempW[[1]])
        matW<-rbind(matW,t(tempW[[1]]))
        aceitaW<-c(aceitaW,tempW[[2]])

        #MeanM
        varMpriori<-solve(solve(gSigma(b_alpha,v_alpha,S))+solve(BMpri))
        meanMpriori<-varMpriori%*%(solve(gSigma(b_alpha,v_alpha,S))%*%M+solve(BMpri)%*%AMpri)
        MeanM<-mvrnorm(1,meanMpriori,varMpriori)
        matMeanM<-rbind(matMeanM,t(MeanM))

        #Para M
        tempM=amostrarM(W,M,Q,x,fm,b_alpha,v_alpha,S,MeanM)
        M<-as.numeric(tempM[[1]])
        matM<-rbind(matM,t(tempM[[1]]))
        aceitaM<-c(aceitaM,tempM[[2]])

        #MeanQ
        varQpriori<-solve(solve(gSigma(b_beta,v_beta,S))+solve(BQpri))
        meanQpriori<-varQpriori%*%(solve(gSigma(b_beta,v_beta,S))%*%Q+solve(BQpri)%*%AQpri)
        MeanQ<-mvrnorm(1,meanQpriori,varQpriori)
        matMeanQ<-rbind(matMeanQ,t(MeanQ))

        #Para Q
        tempQ=amostrarQ(W,M,Q,x,fq,b_beta,v_beta,S,MeanQ)
        Q<-as.numeric(tempQ[[1]])
        matQ<-rbind(matQ,t(tempQ[[1]]))
        aceitaQ<-c(aceitaQ,tempQ[[2]])

        # #Para PSI_lambda
        varPsiW=solve(solve(BW)+t(Xw)%*%solve(gSigma(b_lambda,v_lambda,S))%*%Xw)
        medPsiW=(t(AW)%*%solve(BW)+t(W)%*%solve(gSigma(b_lambda,v_lambda,S))%*%Xw)%*%varPsiW
        PSI_lambda=as.matrix(mvrnorm(1,medPsiW,varPsiW))
        psi_W<-rbind(psi_W,t(PSI_lambda))

        # #Para b_lambda
        tempBW=amostrarb(W,v_lambda,b_lambda,S,abW,bbW,Xw,PSI_lambda,uW)
        b_lambda<-as.numeric(tempBW[[1]])
        aceitaWB<-c(aceitaWB,tempBW[[2]])
        vetWB<-c(vetWB,tempBW[[1]])
        #
        #Para b_alpha
        tempBM=amostrarbMQ(M,v_alpha,b_alpha,S,abM,bbM,MeanM,uM)
        b_alpha<-as.numeric(tempBM[[1]])
        aceitaMB<-c(aceitaMB,tempBM[[2]])
        vetMB<-c(vetMB,tempBM[[1]])

        #Para b_beta
        tempQB=amostrarbMQ(Q,v_beta,b_beta,S,abQ,bbQ,MeanQ,uQ)
        b_beta<-as.numeric(tempQB[[1]])
        aceitaQB<-c(aceitaQB,tempQB[[2]])
        vetQB<-c(vetQB,tempQB[[1]])

        # #Para v_lambda
        RRW=gCorr(b_lambda,S)
        aW=(n/2)+avW
        bW=0.5*t(W-Xw%*%as.matrix(PSI_lambda))%*%solve(RRW)%*%(W-Xw%*%as.matrix(PSI_lambda))+bvW
        v_lambda=1/rgamma(1,shape=aW, rate = bW)
        VetV_lambda<-c(VetV_lambda,v_lambda)
        #
        #Para v_alpha
        RRM=gCorr(b_alpha,S)
        aM=(n/2)+avM
        bM=0.5*t(M-MeanM)%*%solve(RRM)%*%(M-MeanM)+bvM
        v_alpha=1/rgamma(1,shape=aM, rate = bM)
        VetV_alpha<-c(VetV_alpha,v_alpha)

        #Para v_beta
        RRQ=gCorr(b_beta,S)
        aQ=(n/2)+avQ
        bQ=0.5*t(Q-MeanQ)%*%solve(RRQ)%*%(Q-MeanQ)+bvQ
        v_beta=1/rgamma(1,shape=aQ, rate = bQ)
        VetV_beta<-c(VetV_beta,v_beta)

      }else{

        #Para W
        tempW=amostrarW(W,M,Q,Xw,x,fw,PSI_lambda,b_lambda,v_lambda,Xw[,-1])
        W<-as.numeric(tempW[[1]])

        #MeanM
        varMpriori<-solve(solve(gSigma(b_alpha,v_alpha,S))+solve(BMpri))
        meanMpriori<-varMpriori%*%(solve(gSigma(b_alpha,v_alpha,S))%*%M+solve(BMpri)%*%AMpri)
        MeanM<-mvrnorm(1,meanMpriori,varMpriori)

        #Para M
        tempM=amostrarM(W,M,Q,x,fm,b_alpha,v_alpha,S,MeanM)
        M<-as.numeric(tempM[[1]])

        #MeanQ
        varQpriori<-solve(solve(gSigma(b_beta,v_beta,S))+solve(BQpri))
        meanQpriori<-varQpriori%*%(solve(gSigma(b_beta,v_beta,S))%*%Q+solve(BQpri)%*%AQpri)
        MeanQ<-mvrnorm(1,meanQpriori,varQpriori)

        #Para Q
        tempQ=amostrarQ(W,M,Q,x,fq,b_beta,v_beta,S,MeanQ)
        Q<-as.numeric(tempQ[[1]])

        # #Para PSI_lambda
        varPsiW=solve(solve(BW)+t(Xw)%*%solve(gSigma(b_lambda,v_lambda,S))%*%Xw)
        medPsiW=(t(AW)%*%solve(BW)+t(W)%*%solve(gSigma(b_lambda,v_lambda,S))%*%Xw)%*%varPsiW
        PSI_lambda=as.matrix(mvrnorm(1,medPsiW,varPsiW))

        # #Para b_lambda
        tempBW=amostrarb(W,v_lambda,b_lambda,S,abW,bbW,Xw,PSI_lambda,uW)
        b_lambda<-as.numeric(tempBW[[1]])
        #

        #Para b_alpha
        tempBM=amostrarbMQ(M,v_alpha,b_alpha,S,abM,bbM,MeanM,uM)
        b_alpha<-as.numeric(tempBM[[1]])

        #Para b_beta
        tempQB=amostrarbMQ(Q,v_beta,b_beta,S,abQ,bbQ,MeanQ,uQ)
        b_beta<-as.numeric(tempQB[[1]])

        #Para v_lambda
        RRW=gCorr(b_lambda,S)
        aW=(n/2)+avW
        bW=0.5*t(W-Xw%*%as.matrix(PSI_lambda))%*%solve(RRW)%*%(W-Xw%*%as.matrix(PSI_lambda))+bvW
        v_lambda=1/rgamma(1,shape=aW, rate = bW)
        #

        #Para v_alpha
        RRM=gCorr(b_alpha,S)
        aM=(n/2)+avM
        bM=0.5*t(M-MeanM)%*%solve(RRM)%*%(M-MeanM)+bvM
        v_alpha=1/rgamma(1,shape=aM, rate = bM)


        #Para v_beta
        RRQ=gCorr(b_beta,S)
        aQ=(n/2)+avQ
        bQ=0.5*t(Q-MeanQ)%*%solve(RRQ)%*%(Q-MeanQ)+bvQ
        v_beta=1/rgamma(1,shape=aQ, rate = bQ)

      }


    }

    print(j)
  }
  resul=list(matW,aceitaW,matMeanM,matM,aceitaM,matMeanQ,matQ,aceitaQ,psi_W,psi_M,psi_Q,uW,uM,uQ,
             vetWB,aceitaWB,vetMB,aceitaMB,vetQB,aceitaQB,VetV_lambda,VetV_alpha,VetV_beta)

  names(resul)<-c("matW","aceitaW","matMeanM","matM","aceitaM","matMeanQ","matQ","aceitaQ","psi_W","psi_M","psi_Q","uW_sint","uM_sint","uQ_sint",
                  "vetWB","aceitaWB","vetMB","aceitaMB","vetQB","aceitaQB","VetV_lambda","VetV_alpha","VetV_beta")
  return(resul)
}


