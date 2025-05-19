library(abind)
library(CovTools)
library(ks)
library(Riemann)
library(shapes)
library(StatPerMeCo)

####Loading data####
load("ytest.RData")

a1 = read.csv2('1.Approx50.csv')
a2 = read.csv2('2.Cholesky50_50PC.csv')
a3 = read.csv2('3.DCC50.csv')
a4 = read.csv('4.SPDNetM50.csv', header = T)
a5 = read.csv('5.SPDNetF50.csv')
a6 = read.csv('6.SPDNetLE50.csv')
a7 = read.csv('7.SPDNetM50_3lag.csv')
a8 = read.csv('8.SPDNetF50_3lag.csv')
a9 = read.csv('9.SPDNetLE50_3lag.csv')
a10 = read.csv('10.SPDNetM50_5lag.csv')
a11 = read.csv('11.SPDNetF50_5lag.csv')
a12 = read.csv('12.SPDNetLE50_5lag.csv')
a13 = read.csv('13.SPDNetM50_10lag.csv')
a14 = read.csv('14.SPDNetF50_10lag.csv')
a15 = read.csv('15.SPDNetLE50_10lag.csv')
a16 = read.csv('16.SPDNetMHAR.csv')
a17 = read.csv('17.SPDNetLEHAR.csv')
a18 = read.csv('18.SPDNetMHARProc.csv')
a19 = read.csv('19.SPDNetLEHARProc.csv')

predlist = list(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20)

npred = nrow(a5)
nassets = 50
####Computing losses####
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

approximant <- function(A){
  B <- (A+t(A))/2
  C <- (A-t(A))/2
  BB <- tensr::polar(B)
  H <- BB$Z
  speB <- eigen(B)
  L <- speB$values
  Z <- speB$vectors
  L[L<0] <- 0
  Xf <- Z %*% diag(L) %*% t(Z)
  return(Xf)
}

LERM <- function(A, B){
  Aspd = as.matrix(Matrix::nearPD(A)$mat)
  Bspd = as.matrix(Matrix::nearPD(B)$mat)
  lerm = shapes::distcov(Aspd, Bspd, method = 'LogEuclidean')
  return(lerm)
}

mcov = list()
m_list <- vector("list", length(predlist))
names(m_list) <- paste0("m", seq_along(predlist))
lossF = matrix(NA, ncol = length(predlist), nrow = npred)
lossE = matrix(NA, ncol = length(predlist), nrow = npred)
lossP = matrix(NA, ncol = length(predlist), nrow = npred)
lossLE = matrix(NA, ncol = length(predlist), nrow = npred)

for (i in 1:npred){
  mcov[[i]] =  invvech(as.matrix(ytest[i,]))
  for (k in seq_along(predlist)){
    m_list[[k]][[i]] <- invvech(as.matrix(predlist[[k]][i, ]))
    lossF[i,k] = Frobenius(mcov[[i]], m_list[[k]][[i]])
    lossE[i,k] = euc.dist(t(as.matrix(ytest[i,])),as.matrix(predlist[[k]][i,]))
    lossP[i,k] = procOPA(mcov[[i]], m_list[[k]][[i]], scale = F)$rmsd
    lossLE[i,k] = LERM(mcov[[i]], m_list[[k]][[i]], method = 'LogEuclidean')
  }
  print(i)
}
###Frobenius####
colnames(lossF) <- c('Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                     'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                     'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
colMeans(lossF[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                   'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                   'GeoHARLE','GeoHARLEProc',
                   'GeoHARM', 'GeoHARMProc',
                   'Cholesky','DCC', 'RW')])
###Euclidean####
colnames(lossE) <- c('Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                     'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                     'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
colMeans(lossE[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                   'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                   'GeoHARLE','GeoHARLEProc',
                   'GeoHARM', 'GeoHARMProc',
                   'Cholesky','DCC', 'RW')])
###Procrustes####
colnames(lossP) <- c('Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                     'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                     'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
colMeans(lossP[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                   'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                   'GeoHARLE','GeoHARLEProc',
                   'GeoHARM', 'GeoHARMProc',
                   'Cholesky','DCC', 'RW')])

###LogE####
colnames(lossLE) <- c('Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                      'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                      'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
colMeans(lossLE[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                    'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                    'GeoHARLE','GeoHARLEProc',
                    'GeoHARM', 'GeoHARMProc',
                    'Cholesky','DCC', 'RW')])



####MCS####
library(MCS)
set.seed(123)
lossf = lossF[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                  'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                  'GeoHARLE','GeoHARLEProc',
                  'GeoHARM', 'GeoHARMProc',
                  'Cholesky','DCC', 'RW')]
MCSF1 = MCSprocedure(lossf, alpha = 0.1,B=10000,statistic='TR')

losse = lossE[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                  'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                  'GeoHARLE','GeoHARLEProc',
                  'GeoHARM', 'GeoHARMProc',
                  'Cholesky','DCC', 'RW')]
MCSE1 = MCSprocedure(losse, alpha = 0.1,B=10000,statistic='TR')

lossproc = lossP[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                     'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                     'GeoHARLE','GeoHARLEProc',
                     'GeoHARM', 'GeoHARMProc',
                     'Cholesky','DCC', 'RW')]
MCSP1 = MCSprocedure(lossproc, alpha = 0.1,B=10000,statistic='TR')

lossle = lossLE[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                     'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                     'GeoHARLE','GeoHARLEProc',
                     'GeoHARM', 'GeoHARMProc',
                     'Cholesky','DCC', 'RW')]
MCSLE1 = MCSprocedure(lossle, alpha = 0.1,B=10000,statistic='TR')


####Portfolio optimization####
library(data.table)
library(dplyr)
library(ks)
library(lessR)
library(lubridate)
library(Matrix)
library(matrixcalc)
library(mgarch)
library(pracma)
library(quantmod)
library(reshape2)
library(rmgarch)
library(rockchalk)
library(rugarch)
library(shapes)
library(starvars)
library(StatPerMeCo)
library(tensr)
library(tidyr)
library(vars)
require(nloptr)
library(quadprog)


dairetu = read.csv('ret50.csv')
nassets = ncol(dairet)
GMV <- function(Sigma){
  D_matrix <- 2 * as.matrix(nearPD(Sigma)$mat)
  n <- nrow(Sigma)
  d_vector <- rep(0, n)
  A_matrix <- cbind(rep(1, n), diag(n))
  b_vector <- c(1, rep(0, n))
  # use solve.QP to minimize portfolio variance
  quad_prog <- solve.QP(Dmat = D_matrix, dvec = d_vector, Amat = A_matrix, bvec = b_vector, meq = 1) 
  return(quad_prog$solution)
}
mcov <- list()

dairet <- as.matrix(dairetu[(window+1):nrow(dairetu),])
dairetp = dairet/100
pesina <- rep((1/nassets),nassets)
por_ret = matrix(0, ncol = ncol(lossF), nrow = npred)
por_retlong = matrix(0, ncol = ncol(lossF), nrow = npred)
naive <- NULL
turnover = matrix(0, ncol = ncol(lossF), nrow = npred)
turnoverlong = matrix(0, ncol = ncol(lossF), nrow = npred)
turnoverna = rep(0, npred)
turnoverlongna = rep(0, npred) 
iota = as.matrix(rep(1, nassets))
for (i in 1:ntest){
  for (k in seq_along(predlist)){tryCatch({
    pesilong= GMV((m_list[[k]][[i]]))
    pesi = as.matrix(solve(Matrix::nearPD(m_list[[k]][[i]])$mat)%*%iota)/as.numeric((t(iota)%*%solve(Matrix::nearPD(m_list[[k]][[i]])$mat)%*%iota))
    por_ret[i,k] = t(dairet[i,])%*%pesi
    por_retlong[i,k] = t(dairet[i,])%*%pesilong
    if(i >1){
      for(j in 1:nassets){
        turnover[i,k] = turnover[i,k] + abs(pesi[j] - prev_pesi[j] * (1 + dairetp[i-1, j]) / 
                                  as.numeric(1 + t(prev_pesi) %*% dairetp[i-1,])) 
        turnoverlong[i,k] = turnoverlong[i,k] + abs(pesilong[j] - prev_pesilong[j] * (1 + dairetp[i-1, j]) / 
                                      as.numeric(1 + t(prev_pesilong) %*% dairetp[i-1,]))
      }
  
    }
    prev_pesi <- pesi
    prev_pesilong <- pesilong
  }, error=function(e){})
  }
  naive[i] <- t(dairet[i,])%*%pesina
  if(i >1){
    for(j in 1:nassets){
      turnoverna[i] = turnoverna[i] + abs(pesina[j] - pesina[j] * (1 + dairetp[i-1, j]) / 
                                            as.numeric(1 + t(pesina) %*% dairetp[i-1,])) 
      turnoverlongna[i] = turnoverlongna[i] + abs(pesina[j] - pesina[j] * (1 + dairetp[i-1, j]) / 
                                                    as.numeric(1 + t(pesina) %*% dairetp[i-1,]))
    }
  }
  print(i)
}
por_return = as.data.frame(cbind(naive, por_ret))
colnames(por_return) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                          'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                          'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
porreturn = por_return[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                           'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                           'GeoHARLE','GeoHARLEProc',
                           'GeoHARM', 'GeoHARMProc',
                           'Cholesky','DCC', 'RW', 'Naive')]
porreturnlong = as.data.frame(cbind(naive, por_retlong))
colnames(porreturnlong) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                             'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                             'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
porreturnlong = porreturnlong[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                                'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                                'GeoHARLE','GeoHARLEProc',
                                'GeoHARM', 'GeoHARMProc',
                                'Cholesky','DCC', 'RW', 'Naive')]
turnover = as.data.frame(cbind(turnoverna, turnover))
colnames(turnover) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                        'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                        'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
turnover = turnover[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                        'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                        'GeoHARLE','GeoHARLEProc',
                        'GeoHARM', 'GeoHARMProc',
                        'Cholesky','DCC', 'RW', 'Naive')]
turnoverlong = as.data.frame(cbind(turnoverlongna, turnoverlong))
colnames(turnoverlong) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                            'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                            'SPDNetLE10l','GeoHARM','GeoHARLE','GeoHARMProc', 'GeoHARLEProc','RW')
turnoverlong = turnoverlong[,c( 'SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                                'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                                'GeoHARLE','GeoHARLEProc',
                                'GeoHARM', 'GeoHARMProc',
                                'Cholesky','DCC', 'RW', 'Naive')]
colMeans(turnover)
colMeans(turnoverlong)

####Robustness###
library(quantmod)
library(xlsx)
library(highfrequency)
library(xts)
library(ggplot2)
start_date <- as.Date("2007-06-26")
end_date <- as.Date("2021-07-02")

data2017 = read.csv2('SPXUSD2017.csv', header = F)
data2018 = read.csv2('SPXUSD2018.csv', header = F)
data2019 = read.csv2('SPXUSD2019.csv', header = F)
data2020 = read.csv2('SPXUSD2020.csv', header = F)
data2021 = read.csv2('SPXUSD2021.csv', header = F)

data = rbind(data2017, data2018, data2019, data2020, data2021)
colnames(data) = c('date','open', 'high', 'low', 'close', 'adj')
data$date = as.POSIXct(data$date, format = "%Y-%m-%d %H:%M")

df = data[,c('date', 'close')]
df.xts = xts(df$close, order.by = df$date)
rv <- rRVar(df.xts, makeReturns = T)

datercov = as.Date(row.names(dairet), format = '%Y-%m-%d')
datercovtest = datercov[(window+1):nrow(rcovariance)]

SP500test = rv[as.Date(index(rv)) %in% datercovtest,]

thresh = quantile(SP500test, 0.90)

SP500test = as.data.frame(SP500test)
SP500test$regime[SP500test$V1 < thresh] = 'Low'
SP500test$regime[SP500test$V1 >= thresh] = 'High'
SP500test$regime[SP500test$V1 < thresh] = 'Low'

plot(datercovtest, SP500test$V1, type = "l", col = "black", ylim = range(SP500test$V1), ylab = "y2", xlab = "")
for (j in 1:nrow(SP500test)) {
  if (SP500test$regime[j] == 'High') {
    points(datercovtest[j],  SP500test$V1[j], pch = 2, col = "blue")
  } else if (SP500test$regime[j] == 'Low') {
    points(datercovtest[j],  SP500test$V1[j], pch = 15, col = "red") 
  }
}  

####Low regime###
set.seed(123)
#Frobenius
library(MCS)
lossf_low = lossf[SP500test$regime == 'Low',]
round(colMeans(lossf_low),3)
MCSFlow = MCSprocedure(lossf_low, alpha = 0.1,B=10000,statistic='TR')

#Euclidean
losse_low = losse[SP500test$regime == 'Low',]
round(colMeans(losse_low),3)
MCSElow = MCSprocedure(losse_low, alpha = 0.1,B=10000,statistic='TR')

#Procrustes
lossp_low = lossproc[SP500test$regime == 'Low',]
round(colMeans(lossp_low),3)
MCSPlow = MCSprocedure(lossp_low, alpha = 0.1,B=10000,statistic='TR')

#Log-Euclidean
lossle_low = lossle[SP500test$regime == 'Low',]
round(colMeans(lossle_low),3)
MCSLElow = MCSprocedure(lossle_low, alpha = 0.1,B=10000,statistic='TR')


####High regime####
#Frobenius
library(MCS)
set.seed(123)
lossf_high = lossf[SP500test$regime == 'High',]
round(colMeans(lossf_high),3)
MCSFhigh = MCSprocedure(lossf_high, alpha = 0.1,B=10000,statistic='TR')

#Euclidean
losse_high = losse[SP500test$regime == 'High',]
round(colMeans(losse_high),3)
MCSEhigh = MCSprocedure(losse_high, alpha = 0.1,B=10000,statistic='TR')

#Procrustes
lossp_high = lossproc[SP500test$regime == 'High',]
round(colMeans(lossp_high),3)
MCSPhigh = MCSprocedure(lossp_high, alpha = 0.1,B=10000,statistic='TR')

#Log-Euclidean
lossle_high = lossle[SP500test$regime == 'High',]
round(colMeans(lossle_high),3)
MCSLEhigh = MCSprocedure(lossle_high, alpha = 0.1,B=10000,statistic='TR')


#####Portfolio performance####
library(xts)
library(PerformanceAnalytics)
library(PortfolioAnalytics)
library(cvar)
Date1 <- seq(as.Date("2021/8/1"), by = "day", length.out = npred)
porreturn_xts <- xts(porreturn, order.by = Date1)
porreturnlong_xts <- xts(porreturnlong, order.by = Date1)
charts.PerformanceSummary(porreturn_xts, main = "Performance Summary")
performance_table <- table.Stats(porreturn_xts)
colMeans(porreturn)


# Annualized volatility
annualized_vol <- rep(NA, ncol(porreturn))
annualized_vollong <- rep(NA, ncol(porreturnlong))
for(k in 1:ncol(porreturn)){
  annualized_vol[k] = apply.yearly(porreturn_xts[k], FUN = function(x) sqrt(252 * sd(x)))
  annualized_vollong[k] = apply.yearly(porreturnlong_xts[k], FUN = function(x) sqrt(252 * sd(x)))
}
names(annualized_vol) = colnames(porreturn)
annualized_vol
names(annualized_vollong) = colnames(porreturn)
annualized_vollong