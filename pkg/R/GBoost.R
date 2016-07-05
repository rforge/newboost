###change evaluate function for LOSS functions
evaluate_guess_w<-function(guess,gap_temp,weight){
    if (sum(gap_temp)==0){
        return(0)
    }
    obj<-0
    n_d=length(gap_temp)
    for (i in 1:n_d){
        if (gap_temp[i]>0){
            obj<-obj+abs(guess[i]-gap_temp[i])*weight[i]
        }
    }
    obj<-obj/n_d
    return(obj)
}

#' Group Boosting
#' 
#' This functions impelements an algorithm for group L2Boosting with a L1 loss function (median loss).
#' 
#' @param X matrix of regressor variables
#' @param y depedent / repsonse variable
#' @param iter nubmer of iterations
#' @param n_per number of observations per unit (size of gropups)
#' @param weight_p weight function
#' @return A list with different components comparable to \code{L2Boost} is returned.
L1GBoost <- function(X,y, iter,n_per,weight_p){
  p <- dim(X)[2]
  n <- dim(X)[1]
  #n_per is the number of observations...
  #n_per=20 here
  J=n/n_per
  #J is the number of groups, = 66*9 in MM_MAIN
  BetaVector <- rep(0, (p*J))
  #the above record the Beta we found for each time.
  BetaFinal <- rep(c(0,rep(0,p-1)),J)
  #the above BetaFinal is the Final model, initial value is 1
  U <- y-0
  #initialize--subtract 1 first
  stop_rule <- rep(FALSE,iter)
  # boosting algorithm
  f_old <- f_new <- rep(0,n)
  S <- L <- rep(NA, iter) # selected variables
  Res <- f_old # matrix for storing results
  Res_beta <- matrix(0,nrow=p, ncol=1)
  sigma2 <- rep(NA,iter)
  #sigma2 is the loss here...
  # sigma20 <- mean(y^2)
  for (i in 1:iter) {
    #best beta initial is 1
    BETA <- rep(0,(p*J))
    temp_beta<-rep(NA,p)
    SLoss <- rep(0,p)
    temp_loss <- rep(NA,p)
    #DO BOOSTING--GREEDY ALGORITHM TO FIND BEST FEATURE
    for (ll in 1:J){
        current_U = U[((ll-1)*n_per+1):(ll*n_per)]
        current_weight<-weight_p[((ll-1)*n_per+1):(ll*n_per)]
        current_U=matrix(current_U,ncol=1)
        current_weight=matrix(current_weight,ncol=1)
            for (j in 1:p) {
        #this is the residual from the last step...
        #this is the j^tj coordinatewise temp_beta, for GROUP ll
        #solve optimal beta,j and get loss functions...temp_beta[j]
        current_X<-X[((ll-1)*n_per+1):(ll*n_per),j]
        current_X=matrix(current_X,ncol=1)
        ##if there are a lot of 0: just ignore
        ##you would rather guess 1, which is the initial value.
        if (sum(abs(current_X*current_weight)>0)<3){
            #no gain
            temp_beta[j]=0
            # temp_loss[j]=0
            #just no loss gain
            temp_loss[j]<- t(abs(current_U))%*%current_weight
        }
        else{
            # else{
            #Here change fitting with loss function: L1-quantile reg
            fit=rq(current_U~-1+current_X,tau=0.5,weight=current_weight)
            # }
            v=fit$coeff[1]
            temp_beta[j]<-v
            resi=fit$res
            temp_loss[j]<-abs(resi)%*%current_weight
        }
        BETA[[((ll-1)*p+j)]]<-temp_beta[j]
        SLoss[j]<-SLoss[j]+temp_loss[j]
        }
    }
    SLoss<-SLoss/n
    #K is the selected variable--BY Voting among groups.
    L[i] <- k<-which.min(SLoss)
    S[i] <- min(SLoss)
    #now, update Beta_Vector--the incremental at each step
    index_vector =  ((c(1:(p*J)) %% p)== (k%%p))
    for (jj in 1:(J*p)){
        if (jj%% p == (k%%p)){
                BetaVector[jj]=BETA[jj]
        }
        else{
            BetaVector[jj]=0
        }
    }
    #update the final estimator
    BetaFinal <- BetaFinal+BetaVector
    # updating---get best values from Beta_Vector
    beta_values = BetaVector[index_vector]
    update_weight = kronecker(beta_values,rep(1,n_per))
    ##you can also round up here...
    f_new <- f_old + update_weight*X[,k]
    ind <- BetaFinal==0
    U <- y - f_new
    f_old <- f_new
    Res <- cbind(Res, f_new)
    sigma2[i] <- evaluate_guess_w(f_new,y,weight_p)
    if ( i>2 && sigma2[i]/sigma2[i-1] > 0.99) stop_rule[i] <- TRUE # new stoping rule
  }
  return(list(BetaFinal=BetaFinal, BetaVector=BetaVector, L=L,i=i, S=S, Res=Res, Res_beta=Res_beta, sigma2=sigma2, stop_rule=min(which(stop_rule==1))))
}

