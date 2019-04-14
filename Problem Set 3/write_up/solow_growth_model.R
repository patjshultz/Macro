library(Cairo)

#===============================================================================
# Macro PSET 2 Question 3
# Solow growth model
#===============================================================================

# parameters 
delta <- 0.1;
alpha <- 1/3;
A <- 2/5
k0 <- 0.5

#===============================================================================
# Parts 1 and 2
#===============================================================================
# steady state 
capital_ss <- function(A, delta, alpha){
  return((A / (4 * delta))^(1/(1 - alpha)))
}

# verify that A = 2/5 gives kss = 1
capital_ss(A = 2/5, delta = 0.1, alpha = 1/3)

# capital differenc eequation
get_kp <- function(A, delta, k, alpha){
  kp <- 0.25 * (A * k^alpha) + (1 - delta) * k
  return(kp)
}

# set up simulation
k <- k0
periods <- 1:100
capital <- rep(NA, length(periods))

# run simulation
for(i in periods){
 capital[i] <- get_kp(A, delta, k, alpha)
 k <- capital[i]
}

tau <- min(which(capital > 0.9))



Cairo("write_up/solow_growth_alpha_low.png",type = "png", width =800, height = 600)
par(mar = c(5.4, 5.1, 1, 1))
par(ps = 30)
plot(periods, capital, type = "l", lwd = 4)
abline(h = 0.9, v = tau)
text(paste(tau, "periods to k =  0.9"), x = tau + 40, y =0.8)
dev.off()

#===============================================================================
# Parts 3 and 4
#===============================================================================

# recompute steady state and normalize A = 1
capital_ss(A = 2/5, delta = 0.1, alpha = 2/3)

k <- 0.5
# run simulation
for(i in periods){
  capital[i] <- get_kp(A, delta, k, alpha = 2/3)
  k <- capital[i]
}
tau <- min(which(capital > 0.9))



Cairo("write_up/solow_growth_alpha_high.png",type = "png", width =800, height = 600)
par(mar = c(5.4, 5.1, 1, 1))
par(ps = 30)
plot(periods, capital, type = "l", lwd = 4)
abline(h = 0.9, v = tau)
text(paste(tau, "periods to k =  0.9"), x = tau + 25 , y =0.8)
dev.off()
