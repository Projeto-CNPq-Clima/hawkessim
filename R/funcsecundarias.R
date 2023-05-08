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
gSigma <- function(b, v, deff) {
  n <- nrow(deff)
  R <- exp(-b * (as.matrix(dist(deff))))
  mat <- v * R
  mat
}


#######################
#gerando dados da normal multivariada
normalmulti<- function(x, PSI, b, v){

  #valores necessarios para a normal multivariada
  #vetor de médias com as covariaveis
  media<- x%*%PSI

  #matriz de covariancia
  SS<-gSigma(b,v,x)

  #gerando dados de uma normal multivariada
  dados<- MASS::mvrnorm(1, media, SS)

  return(dados)

}

