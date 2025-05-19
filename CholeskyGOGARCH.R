library(data.table)
library(dplyr)
library(fnets)
library(ks)
library(lessR)
library(lubridate)
library(matrixcalc)
library(mgarch)
library(pracma)
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
library(xts)

####Data loading####
#dairet = read.csv('SP100_return.csv')
#corr  = read.csv('SP100_corr.csv')
#variances = read.csv('SP100_var.csv')

row.names(dairet) = dairet[,1]
dairet = dairet[,-1]
row.names(corr) = corr[,1]
corr = corr[,-1]
row.names(variances) = variances[,1]
variances = variances[,-1]

ndays <- nrow(corr)
corr_list <- list()
cov_list <- list()
stock <- c('AAPL','ABT','ACN','ADBE','AMGN','AMT','AMZN','AXP', 'BAC','CAT','CMCSA',
           'COST','CRM','CSCO','CVX',
           'DHR','DIS', 'GE','GOOG', 'GS',
           'HD','IBM','INTU','ISRG','JNJ','JPM','KO', 'LLY','LOW',
           'MA','MCD','MRK','MSFT', 'NFLX',
           'NKE','NVDA','ORCL','PEP','PFE', 'PG','QCOM','T','TJX','TMO','TXN','UNH','UNP',
           'VZ','WFC','WMT')
N = ncol(variances)
rcovariance = matrix(nrow =nrow(corr), ncol = (length(stock)*(length(stock)+1)/2))
dairet = dairet[,colnames(dairet) %in% stock]

for (i in 1:ndays) {
  corrs <- as.numeric(corr[i, ])
  mat <- diag(1, N)  
  mat[upper.tri(mat)] <- corrs  
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)] 
  colnames(mat) = rownames(mat) = colnames(variances)
  corr_list[[i]] = mat
  rcov = diag(sqrt(variances[i,])) %*% mat %*% diag(sqrt(variances[i,]))
  colnames(rcov) = rownames(rcov) = colnames(variances)
  rcov = rcov[colnames(rcov) %in% stock, rownames(rcov) %in% stock]
  cov_list[[i]] = rcov
  rcovariance[i,] = vech(rcov)
  print(i)
}
N = nrow(rcov)

####Utils####
Cholmat <- function (X, tol = sqrt(.Machine$double.eps)) {
  if (!is.numeric(X)) 
    stop("argument is not numeric")
  if (!is.matrix(X)) 
    stop("argument is not a matrix")
  n <- nrow(X)
  if (ncol(X) != n) 
    stop("matrix is not square")
  if (max(abs(X - t(X))) > tol) 
    stop("matrix is not symmetric")
  D <- rep(0, n)
  L <- diag(n)
  i <- 2:n
  D[1] <- X[1, 1]
  if (abs(D[1]) < tol) 
    stop("matrix is numerically singular")
  L[i, 1] <- X[i, 1]/D[1]
  for (j in 2:(n - 1)) {
    k <- 1:(j - 1)
    D[j] <- abs(X[j, j] - sum((L[j, k]^2) * D[k]))+0.00000000001
    if (abs(D[j]) < tol) 
      stop("matrix is numerically singular")
    i <- (j + 1):n
    L[i, j] <- (X[i, j] - colSums(L[j, k] * t(L[i, k, drop = FALSE]) * 
                                    D[k]))/D[j]
  }
  k <- 1:(n - 1)
  D[n] <- abs(X[n, n] - sum((L[n, k]^2) * D[k]))
  if (abs(D[n]) < tol) 
    stop("matrix is numerically singular")
  (L %*% diag(sqrt(D)))
}
isPSD <- function(x, tol = 0.00000001) {
  if (!is.square.matrix(x)) 
    stop("argument x is not a square matrix")
  if (!isSymmetric(x)) 
    stop("argument x is not a symmetric matrix")
  if (!is.numeric(x)) 
    stop("argument x is not a numeric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < tol) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues < 0)) {
    return(FALSE)
  }
  return(TRUE)
}
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
window = 2386
ntest = 1023
eigen1 = matrix(nrow = ndays, ncol = N) 
for(l in 1:ndays){
  eigen1[l,] = eigen(invvech(as.matrix(rcovariance[l,])))$values
}
colMeans(eigen1)
ntilde = N*(N+1)/2

####Cholesky####
library(highfrequency)
cholfact <- matrix(NA, nrow = ndays, ncol = ntilde)
tmp <- matrix(nrow = ndays, ncol = ntilde)
for(i in 1:ndays){
  ch1 <- Cholmat(invvech(as.matrix(rcovariance[i,])),tol = 0.0000000000000000000000000001)
  cholfact[i,] <- t(vech(ch1))
  tmp[i,] <- vech(ch1 %*% t(ch1))
}
pred <- matrix(NA, nrow = ntest, ncol = ntilde)
pos <- rep(NA, ntest)
for(i in 1:ntest){
  pca_res <- prcomp(cholfact[i:(i+window-1),], scale. = TRUE)
  num_pc <- 50
  pc_data <- pca_res$x[, 1:num_pc] 
  var <- vars::VAR(pc_data, p = 1, type = 'const')
  foreca <- rbindlist(lapply(predict(var, n.ahead = 1)[[1]], as.data.table))
  loadings <- pca_res$rotation[, 1:num_pc]  
  means <- pca_res$center                   
  sds <- pca_res$scale                     
  # Inverting the Cholesky decomposition
  original_forecast <- foreca$fcst %*% t(loadings)  
  original_forecast <- sweep(original_forecast, 2, sds, FUN = "*")
  original_forecast <- sweep(original_forecast, 2, means, FUN = "+")
  l2 <- t(vech2mat(original_forecast, lowerOnly = T))
  pred2[i,] <- t(vech(crossprod(l2)))
  pos2[i] <- isPSD(invvech(pred2[i,]))
  print(i)
}


####GO-GARCH####
spec = gogarchspec(mean.model = list(demean = "constant"),
                   variance.model = list(model = "sGARCH", garchOrder = c(1,1), submodel = NULL),
                   ica = "fastica")

#Estimation
F_DCC <- matrix(ncol = ncol(rcovariance), nrow = ntest)
for (i in 1:ntest){
  try({
  dCC.fit = gogarchfit(spec = spec, data = dairet[i:(window-1+i),],
                   out.sample = 1, gfun = "tanh")
  dCC.fcst <- gogarchforecast(dCC.fit, n.ahead = 1)
  rolldCC <- as.matrix(as.data.frame(dCC.fcst@mforecast$U)) 
  F_DCC[i,] <- vech(rolldCC)
  print(i)
  }, silent = TRUE)
}
DCCdf = as.data.frame(F_DCC)
DCCpred = na.locf(DCCdf, na.rm = FALSE)

