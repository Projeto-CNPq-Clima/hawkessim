test_that("prchawkes works", {
  #numero de pontos na malha
  n<-3
  #vetor de covariaveis nxp
  PSI_lambda<-as.matrix(c(1,-3,-4))
  #correlação (constante positiva)
  b_lambda<-1
  #variância (constante positiva)
  v_lambda<-2
  #vetor de covariaveis nxp
  mu_alpha<-rep(.1,9)
  #correlação (constante positiva)
  b_alpha<-1
  #variância (constante positiva)
  v_alpha<-2
  #vetor de covariaveis nxp
  mu_beta<-rep(2,9)
  #correlação (constante positiva)
  b_beta<-1
  #variância (constante positiva)
  v_beta<-2
  #hriznte/temp
  horizon<-30


  hawkessim::prchawkes(n,PSI_lambda, b_lambda,v_lambda,mu_alpha,b_alpha,v_alpha,mu_beta,b_beta,v_beta,horizon)
})
