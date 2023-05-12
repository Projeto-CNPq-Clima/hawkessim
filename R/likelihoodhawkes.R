#'Estimação de Verossimilhança do Simulador do processo de Hawkes
#'
#' Esta é a função para estimar a função de Verossimilhança dos dados a partir de um processo de Hawkes
#'
#' @param dados2 dados do processo
#' @param M Vetor de processos Gaussianos relacionado a lambda
#' @param W Vetor de processos Gaussianos relacionado a alpha
#' @param Q Vetor de processos Gaussianos relacionado a beta
#' @param horizon número positivo, horizonte dos dados.
#'
#' @export
likelihoodhawkes<-function(dados2,M,W,Q,horizon){

  termo<-rep(0,ncol(dados2))
  k=2
  for (k in 1:ncol(dados2)) {
    lambda<-exp(W[k])
    alpha<-exp(M[k])
    beta<-exp(Q[k])
    dadosteste<-as.numeric(na.omit(dados2[,k]))
    termossoma<-rep(0,(length(dadosteste)-1))

    for (i in 2:(length(dadosteste))) {

      termossoma[i-1]<-lambda*alpha*(sum(exp(-beta*(dadosteste[i]-dadosteste[1:(i-1)]))))

    }

    termo1<-sum(log(termossoma))

    termo2<-(alpha/beta)*sum(1-exp(-beta*(horizon-dadosteste)))

    termo[k]<-termo1-lambda*horizon-termo2
  }
  return(termo)
}
