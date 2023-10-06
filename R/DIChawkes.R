#' Deviance Information Criterion
#' The deviance information criterion (DIC) is a hierarchical
#' modeling generalization of the Akaike information criterion (AIC).
#' It is particularly useful in Bayesian model selection problems where the posterior distributions of the models have been obtained by Markov chain Monte Carlo (MCMC) simulation. DIC is an asymptotic approximation as the sample size becomes large, like AIC. It is only valid when the posterior distribution is approximately multivariate normal.
#'
#' @param xdata data from a hawkes process
#' @param matW W parameter matrix resulting from MCMC
#' @param matM M parameter matrix resulting from MCMC
#' @param matQ Q parameter matrix resulting from MCMC
#'
#' @return DIC result
#' @export
DIChawkes<- function(xdata,matW,matM,matQ){

EmatW<-exp(matW)
EmatM<-exp(matM)
EmatQ<-exp(matQ)

D_theta<-NULL

for (i in 1:nrow(EmatQ)) {
  D_theta[i]<- -2*likelihoodhawkes(xdata,EmatW[i,],EmatM[i,],EmatQ[i,])
}

EmatW_hat<-colMeans(EmatW)
EmatM_hat<-colMeans(EmatM)
EmatQ_hat<-colMeans(EmatQ)

D_thetahat<--2*likelihoodhawkes(xdata,EmatW_hat,EmatM_hat,EmatQ_hat)

pD=mean(D_theta)-D_thetahat
dic<- pD+mean(D_theta)

return(dic)
}

