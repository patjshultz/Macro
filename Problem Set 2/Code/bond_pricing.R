bond_price <- function(gamma, g, sigma, j, beta){
  return(beta^j * (exp(j*(-gamma *g + 0.5 * gamma^2 *sigma^2))))
}

bond_yield <- function(gamma, g, sigma, j, beta){
  Qb <- beta^j * (exp(j*(-gamma *g + 0.5 * gamma^2 *sigma^2)))
  return(Qb^(-1/j))
}

maturities <- seq(0.1, 10, by = 0.1)
gamma <- 2
sigma <- 0.1
beta <- 0.97
g <- 0.05

plot(maturities, bond_yield(gamma = 10, g = 0.03, sigma = 0.1, maturities, beta), type= "l")
plot(maturities, bond_yield(gamma = 10, g = 0.05, sigma = 0.1, maturities, beta), type= "l")



# pricing with an AR(1) growth process 
maturity <- 25
phi <- 0.9
gbar = 0.03
sigma_e <- 0.1
sigma_y <- 0.05
g_t <- 0.04

bond_price_ar1 <- function(maturity, gamma,beta,  phi, gbar, sigma_e, sigma_y, g_t){
  sum_mean_term <- 0
  sum_epsilon_term <- 0

  for(k in 1:maturity) {
      sum_mean_term <- sum_mean_term + phi^k * (1 - phi) * gbar
  }
    
  for(k in 1:maturity) {
    sum_epsilon_term <- sum_epsilon_term + phi^(2*(k - 1)) + sigma_e^2
  }
    
    
  X <- -gamma * sum_mean_term - gamma  * phi^maturity * g_t + 0.5 * gamma^2 * sum_epsilon_term #+ 0.5 * gamma^2 * sigma_y^2
  bond_price <- beta^maturity * exp(X)
  yield <- bond_price^(-1/maturity)
  return(yield)
}

yields <- bond_price_ar1(maturity = 1:20, gamma = 3, beta = 0.5, phi = 0.4, gbar= -0.1, sigma_e = 0.3, sigma_y = 0.1, g_t = 0.03)
plot(1:20, yields, type = "l")
