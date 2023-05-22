#'Estimação de Verossimilhança do Simulador do processo de Hawkes
#'
#' Esta é a função para estimar a função de Verossimilhança dos dados a partir de um processo de Hawkes
#'
#' @param data0 dados do processo
#' @param alpha0 Vetor de parâmetros do processo
#' @param beta0 Vetor de parâmetros do processo
#' @param lambda0 Vetor de parâmetros do processo
#'
#' @export
likelihoodhawkes<-function(data0,lambda0,alpha0,beta0){

  termo<-rep(0,ncol(data0))
  for (k in 1:ncol(data0)) {
    lambda<-lambda0[k]
    alpha<-alpha0[k]
    beta<-beta0[k]
    dadosteste<-as.numeric(na.omit(data0[,k]))
    termossoma<-rep(0,(length(dadosteste)-1))

    for (i in 2:(length(dadosteste))) {

      termossoma[i-1]<-lambda+alpha*(sum(exp(-beta*(dadosteste[i]-dadosteste[1:(i-1)]))))

    }

    termo1<-sum(log(termossoma))+log(lambda)

    termo2<-(alpha/beta)*sum(1-exp(-beta*(dadosteste[length(dadosteste)]-dadosteste)))

    termo[k]<-termo1-lambda*dadosteste[length(dadosteste)]-termo2
  }
  termofinal<-sum(termo)
  return(termofinal)
}


