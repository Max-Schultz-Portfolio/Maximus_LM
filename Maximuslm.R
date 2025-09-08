maximuslm = function(xk, y, mtype = NULL, conf_level = NULL , xfocus = NULL, xremoved = NULL, delta = NULL){
  ## initialization
  n = length(y)
  xk = cbind(xk)
  if (is.null(dim(xk)[2])){ 
    p = 1
  }
  else{
    p = dim(xk)[2]
  }
  X = cbind(1,xk)
  if (length(y) != nrow(X)) {
    stop("The length of 'y' must match the number of rows in 'xk'")
  }
  if (mtype == 'ridge' && (is.null(delta) || delta <= 0)) {
    stop("Please provide a positive value for 'delta' for Ridge regression.")
  }
  
  ## helper function 
  s <- function(x){
    return((x-mean(x))/sd(x))
  }
  
  # models 
  if (mtype == 'lasso') {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Please install 'glmnet'")
    }
    library(glmnet)
    lamideal <- cv.glmnet(xk, y, alpha = 1)$lambda.min
    m <- glmnet(xk, y, alpha = 1, lambda = lamideal)
    betahat <- coef(m)
    yhat = X%*%betahat
    residuals = y - yhat 
    se_bhat = 'Undefined for lasso regression in this course'
    betahat_s = 'Undefined for lasso in this course -> since variances are different'
    lamideal2 <- cv.glmnet(s(xk), y, alpha = 1)$lambda.min
    m2 <- glmnet(s(xk), y, alpha = 1, lambda = lamideal2)
    zetahat <- coef(m2)
    H = 'Undefined for lasso in this course'
    lev = 'Undefined for lasso in this course'
    sresid = 'Undefined for lasso in this course'
  } else if (mtype == 'ridge' && (!is.null(delta) && delta != 0)){ 
    XtX = t(X) %*% X
    ridge = XtX + delta * diag(p + 1)  
    betahat = solve(ridge) %*% t(X) %*% y 
    H = X%*%solve(ridge)%*%t(X)
    lev = diag(H)
    yhat = X%*%betahat
    residuals = y - yhat 
    SSE = sum(residuals^2)
    MSE = SSE/(n-p-1)
    sig2 = MSE  
    sresid = (residuals -0)/ (sqrt(MSE)*sqrt(1-lev)) 
    se_bhat = 'Undefined for ridge regression in this course'
    betahat_s = 'Undefined for ridge in this course -> since variances are different'
    X_s = cbind(1,s(xk))
    zetahat = solve(t(X_s)%*%X_s + delta * diag(p + 1)) %*% t(X_s) %*% y
    ps_det <- det(ridge)
  }else if (mtype == 'ls' || delta == 0){
    betahat = solve(t(X) %*% X) %*% t(X) %*% y  
    yhat = X%*%betahat
    residuals = y - yhat 
    SSE = sum(residuals^2)
    MSE = SSE/(n-p-1)
    sig2 = MSE  
    sig_bhat = sig2 * solve(t(X)%*%X) 
    se_bhat = sqrt(diag(sig_bhat)) 
    H = X%*%solve(t(X)%*%X)%*%t(X) 
    lev = diag(H)
    betahat_s = betahat / se_bhat
    X_s = cbind(1,s(xk))
    zetahat = solve(t(X_s) %*% X_s) %*% t(X_s) %*% y    
    ps_det = sqrt(abs(det(t(X)%*%X)))
    sresid = (residuals -0)/ (sqrt(MSE)*sqrt(1-lev)) 
  } else {
    stop('Model type not supported')
  }
  
  # metrics that generally stay the same 
  SSE = sum(residuals^2)
  MSE = SSE/(n-p-1)
  sig2 = MSE  
  SST = var(y) * (n-1)
  MST = SST/ (n-1)
  SSM = SST-SSE
  MSM = SSM/p
  Fstat = MSM/MSE
  RMSE = sqrt(SSE/n)
  pvalue = pf(Fstat, p, n-p-1, lower.tail = F)
  r2 = 1- SSE/SST
  r2adj = 1- MSE/MST 
  cor_matrix <- cor(cbind(y,xk))
  
  if (!is.null(conf_level) && (is.null(delta) || delta == 0) && mtype != 'lasso'){ # there could be some tiny edge case here since ls is a special case of ridge
    t_crit = qt((1 + conf_level) / 2, df = n - p - 1)
    conf_ints = matrix(0, nrow = p + 1, ncol = 2)
    for (i in 1:(p + 1)) {
      conf_ints[i, 1] = betahat[i] - t_crit * se_bhat[i]  
      conf_ints[i, 2] = betahat[i] + t_crit * se_bhat[i]
    }
  } else{
    conf_ints = 'Unable to compute... no confidence level given... or ridge with delta > 0 / lasso'
  }
  
  myvif <- c()
  mymci <- c()
  if (p == 1) {
    myvif <- append(myvif, 1)
    mymci <- append(mymci, 1)
  } else if (mtype == 'ls' || (mtype == 'ridge' && delta == 0)){
    for (k in 1:p) {
      Xp <- X[, -c(k + 1)]
      Yp <- X[, k + 1]
      betahat_k <- solve(t(Xp) %*% Xp) %*% t(Xp) %*% Yp
      yhat_k <- Xp %*% betahat_k
      sse_k <- sum((Yp - yhat_k)^2)
      sst_k <- sum((Yp - mean(Yp))^2)
      r2_k <- 1 - sse_k / sst_k  
      vif_k <- 1 / (1 - r2_k)
      mci_k <- sqrt(vif_k)
      myvif <- append(myvif, vif_k)
      mymci <- append(mymci, mci_k)
      }
  } else if (mtype == 'ridge' && delta > 0) {
    for (k in 1:p) {
      Xp <- X[, -c(k + 1)]
      Yp <- X[, k + 1]
      XtX_minus_k <- t(Xp)%*%Xp
      betahat_k <- solve(t(Xp)%*%Xp + delta * diag(ncol(Xp)))%*% t(Xp) %*% Yp
      yhat_k <- Xp %*% betahat_k  
      sse_k <- sum((Yp - yhat_k)^2)
      sst_k <- sum((Yp - mean(Yp))^2)
      r2_k <- 1 - sse_k / sst_k
      vif_k <- 1 / (1 - r2_k)
      mci_k <- sqrt(vif_k)
      myvif <- append(myvif, vif_k)
      mymci <- append(mymci, mci_k)
    }
  } else if (mtype == 'lasso'){
    myvif <- 'Undefined for lasso regression in this course'
    mymci <- 'Undefined for lasso regression in this course'
  } else{
    myvif <- 'Modeling Error - type not supported, etc.'
    mymci <- 'Modeling Error - type not supported, etc.'
  }

  
  av_m <- numeric(p)
  names(av_m) <- colnames(xk)
  if (p > 1) {
    if (mtype == 'ridge' && (!is.null(delta) && delta != 0)) {
      for (j in 1:p) {
        X_minus_j <- X[, -c(j + 1), drop = FALSE]
        x_j <- X[, j + 1]
        H_ridge <- X_minus_j %*% solve(t(X_minus_j) %*% X_minus_j + delta * diag(p)) %*% t(X_minus_j)
        e_y <- y - H_ridge %*% y
        e_x <- x_j - H_ridge %*% x_j
        av_m[j] <- sum(e_x * e_y) / sum(e_x^2)
      }
    } else if (mtype == 'lasso') {
      av_m <- 'AV slopes not computed for Lasso regression in this course.'
    } else if (delta == 0 || mtype == 'ls'){  
      for (j in 1:p) {
        X_minus_j <- X[, -c(j + 1), drop = FALSE]
        x_j <- X[, j + 1]
        H_minus_j <- X_minus_j %*% solve(t(X_minus_j) %*% X_minus_j) %*% t(X_minus_j)
        e_y <- y - H_minus_j %*% y
        e_x <- x_j - H_minus_j %*% x_j
        av_m[j] <- sum(e_x * e_y) / sum(e_x^2)
      }
    } else {
      av_m <- 'Computational error -> function arguments, etc.'
    }
  } else {
    av_m <- 'AV slopes not relevant in the case of SLR'
  }
  
  
  partial_corr <- numeric(p)
  if (mtype == 'lasso') {
    partial_corr <- "Partial correlation not supported for Lasso regression."
  } else if (p > 1 && !is.null(xfocus) && !is.null(xremoved) && mtype == 'ls') {
    xremoved = if (is.vector(xremoved)) matrix(xremoved, ncol = 1) else xremoved
    X2 <- cbind(1, xremoved)
    betahat_y <- solve(t(X2) %*% X2) %*% t(X2) %*% y
    yhat <- X2 %*% betahat_y  
    residuals_y <- y - yhat 
    betahat_xfocus <- solve(t(X2) %*% X2) %*% t(X2) %*% xfocus
    xfocus_hat <- X2 %*% betahat_xfocus  
    residuals_xfocus <- xfocus - xfocus_hat 
    partial_corr <- cor(residuals_y, residuals_xfocus)
  } else if (p > 1 && !is.null(xfocus) && !is.null(xremoved) && mtype == 'ridge'){
    xremoved <- if (is.vector(xremoved)) matrix(xremoved, ncol = 1) else xremoved
    X2 <- cbind(1, xremoved)
    ridge_matrix <- t(X2) %*% X2 + delta * diag(ncol(X2))
    betahat_y <- solve(ridge_matrix) %*% t(X2) %*% y
    yhat <- X2 %*% betahat_y  
    residuals_y <- y - yhat  
    betahat_xfocus <- solve(ridge_matrix) %*% t(X2) %*% xfocus
    xfocus_hat <- X2 %*% betahat_xfocus  
    residuals_xfocus <- xfocus - xfocus_hat  
    partial_corr <- cor(residuals_y, residuals_xfocus)
  } else {
    partial_corr <- "Partial correlation not meaningful for SLR."
  }
  
  
  
  
  
  results = list(
    '# predictors' = p, 'bhat' = betahat,'bhat_s' = betahat_s, 'zetahat' = zetahat, 'Yhat' = yhat, 'bhat SE' = se_bhat,'CI coeff.' = conf_ints,'SSE' = SSE, 'SSM' = SSM,
    'SST' = SST,'F statistic' = Fstat,'P-value' = pvalue,'r2' = r2,'r2adj' = r2adj, 'RMSE' = RMSE, 'Hat matrix' = H, 'residuals' = residuals,
    'Std. residuals' = sresid,'Leverages' = lev,'RMSE' = RMSE,'Design_Matrix' = X,'VIF' = myvif,
    'MCI' = mymci,'Correlation_Matrix (R)' = cor_matrix,'Pseudo-Determinant' = ps_det, 'Partial Correlation' = partial_corr,
    'A.V. Slopes' = av_m)
  return(results)
}


