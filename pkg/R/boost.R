#' Implementation of the L2Boosting algorithm
#' 
#' This function implements a classical version of the L2Boosting algorithm.
#' @param X matrix of regressors
#' @param y dependent variable
#' @param iter number of iterations
#' @param beta.start initial value for the algorithm to start with
#' @return The functions returns an object of class \code{L2Boost} with the following components:
#' \itemize{
#' \item{BetaFinal}{the estimated final parameter vector.}
#' \item{BetaVector}{vector where each component gives the values of the estimated coefficient in each round.}
#' \item{Res}{matrix with the estimated function values in each round}
#' \item{Res_beta}{matrix with the estimated coefficient vectors from each round}
#'  \item{i}{Iteration when the process was terminated.}
#'  \item{S}{vector of the indices of the selected variables in each round}
#'   \item{sigma2}{estimation of the variance} 
#'   \item{stop_rule}{When the stopping rule, stopped the algorithm.}
#'   \item{iter}{maximal number of iterations of the algorithm}
#'   \item{y}{model response / dependent variable}
#'   \item{X}{model matrix}
#' }
L2Boost <- function(X,y, iter=200, beta.start = rep(0,dim(X)[2])) {
  X <- as.matrix(X)
  p <- dim(X)[2]
  n <- dim(X)[1]
  U <- y - X%*%as.vector(beta.start)
  stop_rule <- rep(FALSE,iter)
  # boosting algorithm
  f_old <- fnew <- X%*%beta.start
  S <- L <- rep(NA, iter) # selected variables
  Res <- f_old # matrix for storing results
  Res_beta <- matrix(beta.start, nrow=p, ncol=1)
  BetaVector <- rep(0, iter)
  BetaFinal <- Res_beta #<- matrix(0,nrow=p, ncol=1)
  sigma2 <- rep(NA,iter)
  sigma20 <- mean(y^2)
  for (i in 1:iter) {
    BETA <- rep(NA,p)
    SLoss <- rep(NA,p)
    for (j in 1:p) {
      BETA[j] <- sum(U*X[,j])/sum(X[,j]^2)
      SLoss[j] <- sum((U- BETA[j]*X[,j])^2)
    }
    S[i] <- k <- which.max(abs(BETA))
    L[i] <- which.min(SLoss)
    BetaVector[i] <- BETA[k]
    Res_beta <- cbind(Res_beta, Res_beta[,i])
    Res_beta[k,i+1] <- Res_beta[k,i+1] + BETA[k]
    BetaFinal[k] <- BetaFinal[k]+BETA[k]
    # updating
    f_new <- f_old + BETA[k]*X[,k]
    ind <- BetaFinal==0
    U <- y - f_new
    f_old <- f_new
    Res <- cbind(Res, f_new)
    sigma2[i] <- mean((y-f_new)^2)
    if ( i>2 && sigma2[i]/sigma2[i-1] > (1-1.1*log(p)/n)) stop_rule[i] <- TRUE # new stoping rule
  }
  res <- list(BetaFinal=BetaFinal, BetaVector=BetaVector, i=i, S=S, Res=Res, Res_beta=Res_beta,
              sigma20=sigma20, sigma2=sigma2, stop_rule=min(which(stop_rule==1)), iter=iter, y=y, X=X)
  class(res) <- "L2Boost"
  return(res)
}


#' Implementation of the orthogonal L2Boosting algorithm
#' 
#' This function implements an orthogonal projection version of the L2Boosting algorithm.
#' @param X matrix of regressors
#' @param y dependent variable
#' @param iter number of iterations
#' @return The functions returns an object of class \code{L2Boost} with the following components:
#' \itemize{
#' \item{BetaFinal}{the estimated final parameter vector.}
#' \item{BetaFinal0}{the estimated final parameter vector of the orthogonal boosting algorithm.}
#' \item{BetaVector}{vector where each component gives the values of the estimated coefficient in each round.}
#' \item{Res}{matrix with the estimated function values in each round}
#' \item{Res_beta}{matrix with the estimated coefficient vectors from each round}
#'  \item{i}{Iteration when the process was terminated.}
#'  \item{S}{vector of the indices of the selected variables in each round}
#'   \item{sigma2}{estimation of the variance} 
#'   \item{stop_rule}{When the stopping rule, stopped the algorithm.}
#'   \item{iter}{maximal number of iterations of the algorithm}
#'   \item{y}{model response / dependent variable}
#'   \item{X}{model matrix}
#' }
L2BoostOGA <- function(X,y,iter=200) {
  p <- dim(X)[2]
  n <- dim(X)[1]
  U <- y
  stop_rule <- rep(FALSE, iter)
  # boosting algorithm
  f_old <- fnew <- rep(0,n)
  S <- L <- rep(NA, iter) # selected variables
  Res <- f_old # matrix for storing results
  BetaVector <- rep(0, iter)
  BetaFinal <- rep(0,p)
  BetaVectorO <- rep(0, iter)
  BetaFinalO <- list()
  Res_beta <- matrix(0,nrow=p, ncol=1)
  sigma2 <- rep(NA,iter)
  sigma20 <- mean(y^2)
  for (i in 1:iter) {
    BETA <- rep(NA,p)
    SLoss <- rep(NA,p)
    for (j in 1:p) {
      BETA[j] <- sum(U*X[,j])/sum(X[,j]^2)
      SLoss[j] <- sum((U- BETA[j]*X[,j])^2)
    }
    S[i] <- k <- which.max(abs(BETA))
    L[i] <- which.min(SLoss)
    BetaVector[i] <- BETA[k]
    Res_beta <- cbind(Res_beta, Res_beta[,i])
    Res_beta[k,i+1] <- Res_beta[k,i+1] + BETA[k]
    BetaFinal[k] <- BetaFinal[k]+BETA[k]
    # updating
    # projection on the selected variables
    # y = X(X'X)^-1Xy
    Xproj <- X[,S[1:i]]
    BetaFinalO[[i]] <- solve(t(Xproj)%*%Xproj)%*%t(Xproj)%*%as.vector(y)
    f_new <- Xproj%*% BetaFinalO[[i]]   #solve(t(Xproj)%*%Xproj)%*%t(Xproj)%*%as.vector(y)
    ind <- BetaFinal==0
    U <- y - f_new
    f_old <- f_new
    Res <- cbind(Res, f_new)
    sigma2[i] <- mean((y-f_new)^2)
    if ( i>2 && sigma2[i]/sigma2[i-1] > (1-1.1*log(p)/n)) stop_rule[i] <- TRUE # new stoping rule
    if ((i+2)>min(p,n)) break
  }
  res <- list(BetaFinal=BetaFinal, BetaFinalO=BetaFinalO, BetaVector=BetaVector, i=i, S=S, Res=Res, 
              Res_beta=Res_beta, sigma20=sigma20, sigma2=sigma2, stop_rule=min(which(stop_rule==1)), iter=iter, X=X, y=y)
  class(res) <- "L2BoostOGA"
  return(res)
}

#' Function for calculation of different MSEs
#' 
#' This function is mainly used for analyzing the performance of the L2Boost and orthogonal L2Boost algorithms.
#' 
#' @param object object of class \code{L2Boost} or \code{L2BoostOGA}
#' @param yref reference variable with regard to which the MSE is calculated (e.g. the observed dependent variables or some oracle)
#' @param beta.true the ture beta coefficient if known
#' @return The function returns a list with different measures of th
MSECal <- function(object, yref, beta.true=NULL) {
  X <- object$X
  y <- object$y
  iter <- object$iter
  A <- rep(0,iter)
  A_post <- rep(0,iter)
  MSE <- rep(0,iter)
  if (class(object)=="L2Boost") {
  for (i in 1:iter) {
    beta <- object$Res_beta[,i+1]
    f_new <-  X%*%beta
    MSE[i] <- mean((yref-f_new)^2)
    ind <- beta==0
    MSE_post[i] <- mean((yref - predict(lm(y~X[,!ind])))^2) 
    if (!is.null(beta.true)) {
    A[i] <- mean((X%*%beta.true-f_new)^2)
    A_post[i] <- mean((X%*%beta.true-predict(lm(y~X[,!ind])))^2)
    }
  }
  }
  
  
  if (class(object)=="L2BoostOGA") {
    for (i in 1:iter) {
      beta <- object$BetaFinalO[[i]]
      Xproj <- X[,object$S[1:i]]
      f_new <- Xproj%*% beta   #solve(t(Xproj)%*%Xproj)%*%t(Xproj)%*%as.vector(y)
      MSE[i] <- mean((yref-f_new)^2)
      ind <- beta==0
      MSE_post[i] <- mean((yref - predict(lm(y~X[,!ind])))^2)
      if (!is.null(beta.true)) {
      A[i] <- mean((X%*%beta.true-f_new)^2)
      A_post[i] <- mean((X%*%beta.true-predict(lm(y~X[,!ind])))^2)
      }
    }
  }
  if (!is.null(beta.true)) return(list(A = A, MSE = MSE, A_post = A_post))
  else  return(list(MSE = MSE, A_post = A_post))
}