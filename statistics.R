score_test <- function(v, y, test.set = NULL, orthogonalization = FALSE, screen.num=15, scale = FALSE) {
  if(scale) {
    v <- scale(v); y <- scale(y)
  }
  if(is.null(test.set)) {
    test.set <- 1:ncol(v)
  }
  x <- v[,test.set]; z <- v[,-test.set]
  n <- nrow(x); p.beta <- ncol(x); p.gamma <- ncol(z)
  
  ## group testing with L2 norm
  # derive epsilon
  if(p.gamma == 0){
      err <- as.numeric(y)
  }else{
    cv.model <- cv.glmnet(z, y, intercept=FALSE)
    gamma.hat <- as.numeric(coef(cv.model, s="lambda.min"))[-1]
    err <- as.numeric(y-z%*%gamma.hat)
  }

  # derive x.mat
  if((p.gamma >= 1) & orthogonalization & (screen.num > 0)){
    if (screen.num == p.beta){
      screen.index <- 1:p.beta
    }else{
      cor.xz <- cor(x, z)
      screen.index <- order(abs(cor.xz)[,1], decreasing=TRUE)[1:screen.num]
    }
    W.hat <- matrix(0, p.gamma, p.beta)
    for (i in screen.index) {
      cv.model.xz <- cv.glmnet(z, x[,i], intercept=FALSE)
      W.hat[,i] <- as.numeric(coef(cv.model.xz, s="lambda.min"))[-1]
    }
    eta <- x - z %*% W.hat
    x.mat <- eta %*% t(eta)
  }else{
    x.mat <- x %*% t(x) 
  }
  
  # construct test statistic
  err.mat <- outer(err, err, "*")
  Tn <- (sum(err.mat*x.mat) - sum(diag(err.mat*x.mat)))/n
  
  # construct p-value
  tr.Sigx2.hat <- (sum(x.mat^2) - sum(diag(x.mat^2)))/(n*(n-1))
  sigma.hat <- mean(err^2)
  Tn.std <- Tn/(sigma.hat*sqrt(2*tr.Sigx2.hat))

  pval <- 1-pnorm(Tn.std)
  return(pval)
}

QF_Test <- function(v, y, test.set, tau=1) {
  
  # outLabs <- cv.glmnet(z, y, intercept=FALSE)
  outLas <- cv.glmnet(v, y, family = "gaussian", alpha = 1,
                      intercept = T, standardize = T)
  beta.init = as.vector(coef(outLas, s = outLas$lambda.min))
  # tau = c(0, 0.5, 1, 1.5, 2)
  test.set = 1:p.beta
  Est = QF(v, y, G=test.set, A=NULL, model="linear", beta.init=beta.init, tau=tau, verbose=TRUE)
  
  # pval1 = summary(Est)$output.est$`Pr(>|z|)`[1]
  # pval2 = summary(Est)$output.est$`Pr(>|z|)`[2]
  # result <- c(pval1,pval2)
  result <- summary(Est)$output.est$`Pr(>|z|)`
  return(result)
}

ST <- function(X.f, Y.f, sub.size=nrow(X.f)*0.3,ncores = 1) {
  n <- dim(X.f)[1]
  p <- dim(X.f)[2]
  
  n1 <- sub.size
  n0 <- n-floor(n1)
  S1 <- sample(1:n, floor(n1), replace=FALSE)
  X.sub <- X.f[S1,]
  Y.sub <- Y.f[S1]
  cvfit <- cv.glmnet(X.sub, Y.sub, intercept=FALSE)
  cf <- as.numeric(coef(cvfit, s="lambda.min"))[-1]
  # cvfit <- cv.ncvreg(X.sub, Y.sub, penalty=penalty, nfolds=nfolds)
  # model.cvfit <- ncvfit(X.sub, Y.sub, penalty=penalty, lambda=cvfit$lambda.min)
  # cf <- as.numeric(model.cvfit$beta)
  set1 <- (1:p)[abs(cf)>0]
  resi <- Y.sub-X.sub%*%cf
  beta.m <- t(standardize(X.sub[,-set1]))%*%resi
  screen.set <- sort(order(abs(beta.m),decreasing=TRUE)[1:(n0-1-length(set1))])
  a <- (1:p)[-set1]
  screen.set <- union(a[screen.set],set1)
  X <- X.f[-S1,screen.set]
  Y <- Y.f[-S1]
  
  score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
  node <- score.nodewiselasso(X, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                              parallel=TRUE, ncores=ncores, oldschool = FALSE, lambdatuningfactor = 1)
  Theta <- node$out
  return(list(X=X, Y=Y, n0=n0, screen.set=screen.set, Theta=Theta))
}

ST_res <- function(test.set, X, Y, n0, screen.set, Theta, M=500) {
  Gram <- t(X)%*%X/n0
  sreg <- scalreg(X,Y)
  beta.hat <- sreg$coefficients
  sigma.sq <- sum((Y-X%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
  test.set.i <- intersect(screen.set,test.set)
  index <- screen.set%in%test.set.i
  
  Omega <- diag(Theta%*%Gram%*%t(Theta))*sigma.sq
  beta.db <- beta.hat+Theta%*%t(X)%*%(Y-X%*%beta.hat)/n0
  margin.st <- sqrt(n0)*abs(beta.db[index])/sqrt(Omega[index])
  margin.nst <- sqrt(n0)*abs(beta.db[index])
  stat.st <- max(margin.st)
  stat.nst <- max(margin.nst)
  
  stat.boot.st <- stat.boot.nst <- rep(NA,M)
  for (i in 1:M) {
    e <- rnorm(n0)
    xi.boot <- Theta[index,]%*%t(X)%*%e*sqrt(sigma.sq)/sqrt(n0)
    stat.boot.nst[i] <- max(abs(xi.boot))
    stat.boot.st[i] <- max(abs(xi.boot/sqrt(Omega[index])))
  }
  
  # if (stat.nst>quantile(stat.boot.nst,1-alpha)) rej.nst <- "reject" else rej.nst <- "fail to reject"
  # if (stat.st>quantile(stat.boot.st,1-alpha)) rej.st <- "rejct" else rej.st <- "fail to reject"
  # result <- list(stat.nst, rej.nst, stat.st, rej.st)
  # names(result) <- c("non-studentized test","non-studentized test","studentized test","studentized test")
  # rej.nst <- stat.nst > quantile(stat.boot.nst, 1-alpha)
  pval.nst <- sum(stat.boot.nst >= stat.nst)/M
  # rej.st <- stat.st > quantile(stat.boot.st, 1-alpha)
  pval.st <- sum(stat.boot.st >= stat.st)/M
  # result <- c(rej.nst, rej.st)
  # return(result)
  pvalue <- matrix(c(pval.nst, pval.st),nrow = 1,ncol = 2)
  colnames(pvalue) <- c("NST","ST")
  return(pvalue)
}

max_norm_test <- function(v,y,test.set,M=500,sub.size = nrow(v)*0.3,ncores = 1){
  ST.para <- ST(v, y, sub.size=sub.size, ncores = 1)
  return(ST_res(test.set, ST.para$X, ST.para$Y, ST.para$n0, 
                ST.para$screen.set, ST.para$Theta, M=500))
}