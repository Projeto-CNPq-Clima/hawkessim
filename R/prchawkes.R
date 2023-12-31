#'Simulador do processo de Hawkes
#'
#' Esta é a função para simular os dados a partir de um processo de Hawkes
#'
#' @param n um número inteiro positivo, número de pontos na malha
#' @param PSI_lambda vetor de covariaveis px1
#' @param b_lambda número inteiro positivo
#' @param v_lambda número inteiro positivo, variância do parâmetro lambda
#' @param mu_alpha vetor de covariaveis px1
#' @param b_alpha número inteiro positivo
#' @param v_alpha número inteiro positivo, variância do parâmetro alpha
#' @param mu_beta vetor de covariaveis px1
#' @param b_beta número inteiro positivo
#' @param v_beta número inteiro positivo, variância do parâmetro beta
#' @param horizon número positivo, horizonte dos dados.
#'
#' @export
prchawkes<-function(n, PSI_lambda, b_lambda,v_lambda,mu_alpha,b_alpha,v_alpha,
                     mu_beta,b_beta,v_beta,horizon){

  nn<-n^2

  #gerar a malha de pontos simulados (que serão os dados depois - x)
  pontos_simulados<-Gerarmalha(n)

  #incluindo intercepto
  xx<-cbind(rep(1,nn),pontos_simulados)


  ##############lambda###################

  #Vetor W
  W<-normalmulti(xx,PSI_lambda,b_lambda,v_lambda,pontos_simulados)

  #vetor lambda
  lambda<- exp(W)

  ##############alpha#################

  #Vetor M
  M<-normalmultiMQ(mu_alpha,b_alpha,v_alpha,pontos_simulados)

  #vetor alpha
  alpha<- exp(M)


  ###############beta##################

  #Vetor Q
  Q<-normalmultiMQ(mu_beta,b_beta,v_beta,pontos_simulados)

  while (sum(exp(M)>=exp(Q))) {
    Q<-normalmultiMQ(mu_beta,b_beta,v_beta,pontos_simulados)
  }

  #vetor beta
  beta<- exp(Q)




  ################################

  #simulando processo de hawkes
  #vetor vazio para n^2 listas
  lista <- vector("list", length = nn)

  #processo de hawkes
  for (i in 1:nn) {
    lista[[i]]<-hawkes::simulateHawkes(lambda[i],alpha[i],beta[i],horizon)
  }

  #encontrando o maior tamanho
  maior<-0
  for (i in 1:nn) {
    if (length(unlist(lista[[i]])) > maior)
    {
      maior<- length(unlist(lista[[i]]))
    }
  }


  #matriz com todos os dados simulados (cada linha é um processo de Hawkes)
  #matriz vazia com maior tamanho sendo o processo com mais dados
  proc_hawkes<- matrix(NA,nn, nrow=maior)

  #alocando cada processo em uma coluna diferente
  for (i in 1:nn) {
    X<-unlist(lista[[i]])
    tamanho<- maior-length(X)
    proc_hawkes[,i]<-c(X,rep(NA,tamanho))
  }

  #retornando valores
  res <- list(proc_hawkes, W, M, Q, xx)
  names(res)<-c("data process","W","M","Q","X")
  return(res)
}


