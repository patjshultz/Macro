#===============================================================================
# Finance 924: Assignment 4
#===============================================================================

rm(list = ls())
library(quantmod)
library(lubridate)
library(sandwich)
library(ggplot2)
library(mFilter)
theme_set(theme_bw(base_size = 20))

#===============================================================================
# Load data from FRED
#===============================================================================
alpha <- 0.67
beta <- 0.99

# first get consumption data
getSymbols(c("GDPC1", "ND000334Q", "GPDIC1", "GPDICA", "PCECC96", 
             "LFWA64TTUSM647N", "LFWA64TTUSM647S", "LREM64TTUSM156S", 
             "AWHMAN", "AVHWPEUSA065NRUG"), src='FRED') 
dates <- index(GDPC1)
ind <- which(dates>as.Date("1954-12-31"))

Y <- as.numeric(GDPC1); Y_PW <- Y[ind]
I <- as.numeric(GPDIC1); I_PW <- I[ind]
C <- as.numeric(PCECC96); C_PW <- C[ind]

population <- LFWA64TTUSM647N; employment <- LREM64TTUSM156S; hours <- AWHMAN

Qdata <- data.frame(dates, Y, I, C)
Qdata_PW <- Qdata[ind, ]

nQ <- nrow(Qdata)
nQ_PW <- nrow(Qdata_PW)

#===============================================================================
# part a 
#===============================================================================

## calculate average growth rate of economy 
gamma_Y <- (1 + mean(log(Y_PW[2:nQ_PW]/Y_PW[1:(nQ_PW-1)])))^4
gamma_C <- (1 + mean(log(C_PW[2:nQ_PW]/C_PW[1:(nQ_PW-1)])))^4
gamma_I <- (1 + mean(log(I_PW[2:nQ_PW]/I_PW[1:(nQ_PW-1)])))^4
gamma <- (gamma_Y + gamma_C + gamma_I)/3


## calculate the average depreciation rate 
Ibar <- mean(I)
Ybar <- mean(Y)
x <- Ibar/Ybar

numerator <- beta * (1 - alpha) * (gamma - 1) + (beta -1) * x
denominator <- beta *(x - alpha + 1)
delta <- numerator/denominator
delta * 100 

## construct quarterly number of hours worked ##
hours$quarterly <- NA

# calculate quarterly hours worked
for(i in 3:nrow(hours)){
  hours$quarterly[i] <- sum(4 * hours$AWHMAN[i:(i-2)])
}

# plot quarterly hours
Q_index <- which(month(as.Date(index(hours))) %in% c(3, 6, 9, 12))
hoursQ <- data.frame(date = index(hours)[Q_index], hQ = hours$quarterly[Q_index])
colnames(hoursQ) <- c("date", "H")
hoursQ_PW <- hoursQ[which(hoursQ$date>as.Date("1954-12-31")), ]

ggplot(data = hoursQ_PW, aes(x = date, y = H))+
  geom_line(color = "#00AFBB", size = 2) + 
  geom_hline(yintercept = mean(hoursQ_PW$H)) 




## construct a time series of the US capital stock since 1955  ##
capital_accumulation <- function(investment = I, delta = 0.02){
  # investment: vector of investment data
  # delta:  numeric scalar, the depreciation rate
  nobs <- length(investment)
  K <- matrix(data = NA, nrow = nobs, ncol = 1)
  K[1] <- 0
  
  for (t in 1:(nobs-1)) {
   K[t+1] <- investment[t]/4 + (1-delta) * K[t] 
   cat("I:", investment[t]/4, " Kt:", K[t], "Kt1:", K[t+1], "\n")
  }
  return(K)
}


# calculate capital stock for entire sample
K_stock <- capital_accumulation(I, delta = 0.02)
Qdata$K <-as.numeric(K_stock)

# subset to just post war data
ind <- which(Qdata$date>as.Date("1954-12-31"))
K_stock_PW <- K_stock[ind]
Qdata_PW$K <- K_stock_PW


# plot 
df <- melt(Qdata_PW[, c(-2)], id.vars = "dates")
ggplot(df, aes(x=dates, y=value, col=variable)) + geom_line() #+ geom_vline(xintercept = as.Date("2008-09-15"))

## construct solow residuals

# construct quarterly population series
popQ <- population[which(month(index(population)) %in% c(3, 6, 9, 12))]
NQ <- data.frame(N = popQ)
colnames(NQ) <- "N"
NQ$date <- as.Date(rownames(NQ)) 

# data frame with population and hours to get N*H
pop_hours <- merge(NQ, hoursQ, by = "date")
colnames(pop_hours)[1] <- "dates"
pop_hours$NH <- pop_hours$N * pop_hours$H
min(pop_hours$date)

# merge by date vector for NH and Y, I, C, K
Qdata$dates <- as.yearqtr(Qdata$dates)
pop_hours$dates <- as.yearqtr(pop_hours$dates)

solow_data <- merge(x = Qdata, y = pop_hours, by = "dates")
k <- solow_data$K / solow_data$N
solow_residual <- log(solow_data$Y) - (1-alpha)*log(k) - alpha * log(solow_data$N)

TFP_shocks <- data.frame(dates = solow_data$dates,
                         lnA = solow_residual,
                         A = exp(solow_residual),
                         dA = c(NA, diff(exp(solow_residual))))

nshocks <- nrow(TFP_shocks)
ggplot(data = TFP_shocks, aes(x = dates, y = dA))+
  geom_line(color = "#00AFBB", size = 1) 
ggplot(data = TFP_shocks, aes(x = dates, y = A))+
  geom_line(color = "#00AFBB", size = 1) 
ggplot(data = TFP_shocks, aes(x = dates, y = lnA))+
  geom_line(color = "#00AFBB", size = 1) 


lead <- TFP_shocks$lnA[1:(nshocks-1)]
lag <- TFP_shocks$lnA[2:nshocks]
abar <- rep(mean(solow_residual), length(solow_residual)-1)

x <- cbind(abar, lag)
y <- lead

solve(t(lag)%*%lag)%*%t(lag)%*%y

#===============================================================================
# b) Detrend the data using the HP filter
#===============================================================================

colnames(hoursQ_PW) <- c("dates", "H")
# first merge all data from previous section to a single dataframe with a shared date vector
hoursQ_PW$dates <- as.yearqtr(hoursQ_PW$date)
Qdata_PW$dates <- as.yearqtr(Qdata_PW$dates)
str(Qdata_PW); str(TFP_shocks); str(hoursQ_PW); str(pop_hours)

# remove everything from the environment exepct for a dataframe of all qrtly data
tfp_data <- Reduce(merge, list(hoursQ_PW, Qdata_PW, TFP_shocks, pop_hours))
str(tfp_data)
dates <- tfp_data$dates
rm(list=setdiff(ls(), c("dates", "tfp_data", "alpha"))); 
colnames(tfp_data)
nobs <- nrow(tfp_data)

# HP filter of output
hpY <- hpfilter(x = log(tfp_data$Y), type = "lambda", freq = 1600, drift = T)
tcY <- data.frame(date = dates, tY = hpY$trend, cY = hpY$cycle)
ggplot(data = tcY, aes(x = date, y = tY))+
  geom_line(color = "blue", size = 2) +
  ggtitle("Trend GDP")

ggplot(data = tcY, aes(x = date, y = cY))+
  geom_line(color = "blue", size = 2)   +
  ggtitle("Cycle GDP")

# HP filter of consumption
hpC <- hpfilter(x = log(tfp_data$C), type = "lambda", freq = 1600, drift = T)
tcC <- data.frame(date = dates, tC = hpC$trend, cC = hpC$cycle)
ggplot(data = tcC, aes(x = date, y = tC))+
  geom_line(color = "blue", size = 2) +
  ggtitle("Trend Consumption")

ggplot(data = tcC, aes(x = date, y = cC))+
  geom_line(color = "blue", size = 2)   +
  ggtitle("Cycle Consumption")

# HP filter of investment
hpI <- hpfilter(x = log(tfp_data$I), type = "lambda", freq = 1600, drift = T)
tcI <- data.frame(date = dates, tI = hpC$trend, cI = hpI$cycle)
ggplot(data = tcI, aes(x = date, y = tI))+
  geom_line(color = "blue", size = 2) +
  ggtitle("Trend Investment")

ggplot(data = tcI, aes(x = date, y = cI))+
  geom_line(color = "blue", size = 2)   +
  ggtitle("Cycle Investment")

# HP filter of hours
hpH <- hpfilter(x = log(tfp_data$H), type = "lambda", freq = 1600, drift = F)
tcH <- data.frame(date = dates, tH = hpH$trend, cH = hpH$cycle)
ggplot(data = tcH, aes(x = date, y = tH))+
  geom_line(color = "blue", size = 2) +
  ggtitle("Trend Hours")

ggplot(data = tcH, aes(x = date, y = cH))+
  geom_line(color = "blue", size = 2)   +
  ggtitle("Cycle Hours")


# create a data frame of cycle components
cycles <- data.frame(dates = dates,
                     GDP = hpY$cycle, 
                     C = hpC$cycle,
                     I = hpI$cycle)

cycles_long <- melt(cycles, id.vars = "dates")
ggplot(cycles_long, aes(x=dates, y=value, col=variable)) + 
  geom_line() +
  ggtitle("Cycle: GDP, C, I") 
  


#-------------------------------------------------------------------------------
# b1) calculcate symmay statistics of consumption, output, hours, and investment
#-------------------------------------------------------------------------------
summary_tbl <- matrix(data = NA, nrow = 2, ncol = 3)
colnames(summary_tbl) <- c("Cons.", "Inv.", "Hours")

sigmaY <- sd(hpY$cycle) 
sigmaC <- sd(hpC$cycle) 
sigmaI <- sd(hpI$cycle) 
sigmaH <- sd(hpH$cycle)


summary_tbl[1, 1] <- sigmaC/sigmaY
summary_tbl[1, 2] <- sigmaI/sigmaY
summary_tbl[1, 3] <- sigmaH/sigmaY

summary_tbl[2, 1] <- cor(hpY$cycle, hpC$cycle)
summary_tbl[2, 2] <- cor(hpY$cycle, hpI$cycle)
summary_tbl[2, 3] <- cor(hpY$cycle, hpH$cycle)


#-------------------------------------------------------------------------------
# b2/b3) Fit an AR(2)/AR(1) process to the Solow resdiual
#-------------------------------------------------------------------------------

hpa <- hpfilter(x = tfp_data$lnA, type = "lambda", freq = 1600, drift = F)
lnA_cycle <- hpa$cycle
tca <- data.frame(date = dates, ta = hpa$trend, ca = hpa$cycle)
ggplot(data = tca, aes(x = date, y = ta))+
  geom_line(color = "blue", size = 2) +
  ggtitle("Trend Solow Residual")
ggplot(data = tca, aes(x = date, y = ca))+
  geom_line(color = "blue", size = 2) +
  ggtitle("Cycle Solow Residual")


nobs <- length(lnA_cycle)
SR <- lnA_cycle[3:nobs]
SR1 <- lnA_cycle[2:(nobs-1)]
SR2 <- lnA_cycle[1:(nobs-2)]

ar(x = lnA_cycle, order.max = 2)

fit_AR2 <- summary(lm(SR ~ SR1 + SR2 + 0))
acf(fit_AR2$residuals)

ar(x = lnA_cycle, order.max = 1, method = "ols")
fit_AR1 <- summary(lm(SR ~ SR1 + 0))
acf(fit_AR1$residuals)

#===============================================================================
# c) Calibrate psi, so steady state value matches psi in the data
#===============================================================================

hrs_per_q <- 24 * 30.5 * 3
Hss <- mean(tfp_data$H/ hrs_per_q)
YCss <- mean(tfp_data$Y/tfp_data$C)
yss <- mean(tfp_data$Y*10^9/tfp_data$N)

psi_KPR <- alpha * (YCss) * ((1 - Hss)/Hss)
psi_GHH <- alpha * (yss / mean(tfp_data$H)^1.5)

#===============================================================================
# f) download data on government expenditures 
#===============================================================================
getSymbols("W068RCQ027SBEA", src='FRED') 
G <- data.frame(dates = as.yearqtr(index(W068RCQ027SBEA)), G= W068RCQ027SBEA)
colnames(G) <- c("dates", "G")

# merge government data to old data
GY <- merge(G, tfp_data, by = "dates")
nobs <- nrow(GY)
mean(GY$G/GY$Y)


# use the hp filter to calculate the cycle and trend of gov. expendtiture
hpG<- hpfilter(x = log(GY$G), type = "lambda", freq = 1600, drift = F)
tcG <- data.frame(date = GY$dates, tG = hpG$trend, cG = hpG$cycle)
ggplot(data = tcG, aes(x = date, y = tG))+
  geom_line(color = "blue", size = 2) +
  ggtitle("Trend Gov. Exp.")

ggplot(data = tcG, aes(x = date, y = cG))+
  geom_line(color = "blue", size = 2)   +
  ggtitle("Cycle Gov. Exp")

# using HP filtered cycle series, estimate an AR(1) process for detrended level of g
Glead <- hpG$cycle[2:nobs]
Glag <- hpG$cycle[1:(nobs-1)]
fit <- summary(lm(Glead ~ Glag))

sd(fit$residuals)
#
t <- matrix( data = c("a" ,    0.00711,
"y",     0.0126,
"i" ,    1.82,
"c"  ,   0.14,
"h",     0.66,
"w"  ,   0.35), ncol = 2)

