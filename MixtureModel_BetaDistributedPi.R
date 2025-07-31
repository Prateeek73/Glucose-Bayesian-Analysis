model {
  for (i in 1:n) {
    X[i] ~ dcat(c(pi, 1 - pi))
    y[i] ~ dnorm(
      ifelse(X[i] == 1, theta1, theta2), 
      ifelse(X[i] == 1, 1 / s21, 1 / s22))  
  }
  pi ~ dbeta(alpha, beta)
  theta1 ~ dnorm(mu0, 1 / s21)
  theta2 ~ dnorm(mu0, 1 / s22)
  tau1 ~ dgamma(v0/2, v0 * s20 / 2)
  tau2 ~ dgamma(v0/2, v0 * s20 / 2)
  s21 <- 1 / tau1
  s22 <- 1 / tau2
}