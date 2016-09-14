#' Implementation of the L2Boosting algorithm
#' 
#' This function implements a classical version of the L2Boosting algorithm.
#' @param X matrix of regressors
#' @param y dependent variable
#' @param iter number of iterations
#' @param beta.start initial value for the algorithm to start with
#' @param post logical, if post Boosting algorithm should be applied (default \code{FALSE}).
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
#' @export
#' @rdname L2Boost
L2Boost <- function(X,y, iter=200, beta.start = rep(0,dim(X)[2]), post=FALSE) {
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
    BetaFinal[k] <- BetaFinal[k]+BETA[k]
    if (post) {
    ind <- BetaFinal==0
      if (sum(ind)==p) {
        BetaFinalp <- rep(0,p)
        } else {
        BetaFinalp <- rep(0,p)
        BetaFinalpc <- coef(lm(y~ -1 + X[,!ind]))
        BetaFinalp[!ind] <- BetaFinalpc
        }
      Res_beta <- cbind(Res_beta, BetaFinalp)
    } else {
      Res_beta <- cbind(Res_beta, Res_beta[,i])
      Res_beta[k,i+1] <- Res_beta[k,i+1] + BETA[k]
    }
    # updating
    f_new <- f_old + BETA[k]*X[,k]
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

#' L2Boosting with multistart
#' 
#' The L2Boost algorithm is started at different values and then the intersection of all selected variables from each run 
#' is included in the final model.
#' 
#' @export
#' @rdname L2Boost
L2Boost.multistart <- function(X,y, num.start=5, iter=200) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  start.length <- min(5, p)
  Xmax <- apply(X, 2, max)
  Xmin <- apply(X, 2, min)
  finalset <- 1:p
  for (i in 1:num.start) {
    beta <- rep(0,p)
    ind <- sample(p, start.length)
    beta[ind] <- runif(start.length, min=Xmin[ind], max=Xmax[ind])
    PGA <- L2Boost(X,y,iter=iter,  beta.start=beta)
    set <- union(ind, PGA$S[1:PGA$stop_rule])
    finalset <- intersect(finalset, set)
  }
   if (length(finalset)==0) {
     message("No variables selected")
     res <- list(beta=rep(0,p), finalset = finalset)
     return(res)
     } else {
  betaEst <- coef(lm(y~X[,finalset]))[-1]
  beta <- rep(0,p)
  beta[finalset] <- betaEst
  res <- list(beta = beta, betaEst=betaEst, finalset = finalset)
  return(res)
     }
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
#' @export
#' @rdname L2BoostOGA
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
    if (i > 1 && k==S[i-1]) {
      res <- list(BetaFinal=BetaFinal, BetaFinalO=BetaFinalO, BetaVector=BetaVector, i=i, S=S, Res=Res, 
                  Res_beta=Res_beta, sigma20=sigma20, sigma2=sigma2, stop_rule=min(i-1,min(which(stop_rule==1))), iter=iter, X=X, y=y)
      class(res) <- "L2BoostOGA"
      return(res)
    }
    L[i] <- which.min(SLoss)
    BetaVector[i] <- BETA[k]
    Res_beta <- cbind(Res_beta, Res_beta[,i])
    Res_beta[k,i+1] <- Res_beta[k,i+1] + BETA[k]
    BetaFinal[k] <- BetaFinal[k]+BETA[k]
    # updating
    # projection on the selected variables
    # y = X(X'X)^-1Xy
    Xproj <- X[,S[1:i]]
    BetaFinalO[[i]] <- ginv(t(Xproj)%*%Xproj)%*%t(Xproj)%*%as.vector(y) #solve
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
#' @param Xnew design matrix of the observations for which prediction accuracy shall be evaluated.
#' @param beta.true the ture beta coefficient if known
#' @return The function returns a list with different measures of the the predictive accuracy.
#' @export
MSECal <- function(object, yref=NULL, Xnew=NULL, beta.true=NULL) {
  if (is.null(Xnew)) {
    Xnew <- object$X
  }
  if (is.null(yref)) {
    yref <- object$y
  }
  y <- object$y
  iter <- object$iter
  A <- rep(0,iter)
  #A_post <- rep(0,iter)
  MSE <- rep(0,iter)
  #MSE_post <- rep(0,iter)
  if (class(object)=="L2Boost") {
  for (i in 1:iter) {
    beta <- object$Res_beta[,i+1]
    f_new <-  Xnew%*%beta
    MSE[i] <- mean((yref-f_new)^2)
    ind <- beta==0
    #MSE_post[i] <- mean((yref - predict(lm(yref~ -1 + Xnew[,!ind])))^2) 
    if (!is.null(beta.true)) {
    A[i] <- mean((Xnew%*%beta.true-f_new)^2)
    #A_post[i] <- mean((Xnew%*%beta.true-predict(lm(yref ~ -1 + Xnew[,!ind])))^2)
    }
  }
  }
  
  if (class(object)=="L2BoostOGA") {
    for (i in 1:iter) {
      beta <- object$BetaFinalO[[i]]
      Xproj <- Xnew[,object$S[1:i]]
      f_new <- Xproj%*% beta   #solve(t(Xproj)%*%Xproj)%*%t(Xproj)%*%as.vector(y)
      MSE[i] <- mean((yref-f_new)^2)
      ind <- beta==0
      #MSE_post[i] <- mean((yref - predict(lm(y~ -1 +Xnew[,!ind])))^2)
      if (!is.null(beta.true)) {
      A[i] <- mean((Xnew%*%beta.true-f_new)^2)
      #A_post[i] <- mean((Xnew%*%beta.true-predict(lm(y~ -1 + Xnew[,!ind])))^2)
      }
    }
  }
  if (!is.null(beta.true)) return(list(A = A, MSE = MSE))
  else  return(list(MSE = MSE))
}

#' Predict fot L2Boost objects
#'
#' Predicted values based on L2Boost objects
#' @param object Object of class \code{L2Boost}
#' @param newx Matrix of x-variables to be predicted
#' @param mstop Stopping number, i.e. for which iteration of the algorithm the prediction is to be conducted
#' @return A vector with the predicted values.
#' @export
predict.L2Boost <- function(object, newx=NULL, mstop=object$stop_rule) {
  if (is.null(newx)) newx <- object$X
  beta <- object$Res_beta[,mstop+1]
  yhat <- newx%*%beta
  return(yhat)
}

#' Predict fot L2BoostOGA objects
#'
#' Predicted values based on L2BoostOGA objects
#' @param object Object of class \code{L2BoostOGA}
#' @param newx Matrix of x-variables to be predicted
#' @param mstop Stopping number, i.e. for which iteration of the algorithm the prediction is to be conducted
#' @return A vector with the predicted values.
#' @export
predict.L2BoostOGA <- function(object, newx=NULL, mstop=object$stop_rule) {
  if (is.null(newx)) newx <- object$X
  beta <- object$BetaFinalO[[mstop]]
  Xproj <- newx[,object$S[1:mstop]]
  yhat <- Xproj%*% beta
  return(yhat)
}