# DGP
n <- 100
p <- 10
s <- 5
n1 <- 50
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(rep(1,s), rep(0,p-s))
y <- X%*%beta + rnorm(n)
#
Xo <- matrix(rnorm(n1*p), ncol=p)
fo <-  Xo%*%beta
yo <- fo + rnorm(n1)

yhat_oracle <- predict(lm(y~ -1 + X[,1:s]))
beta_oracle <- coef(lm(y~ -1 + X[,1:s]))
# Estimation with L2Boosting

b1 <- L2Boost(X,y,iter=200)
b1p <- L2Boost(X,y,iter=200, post=TRUE)
b2 <- L2Boost(X,y,iter=200, beta.start=rep(1,p))
b3 <- L2BoostOGA(X,y,iter=200)
b3 <- L2BoostOGA(X,y,iter=200)
b4 <- L1GBoost(X,y, iter=200,n_per=10,weight_p=rep(1,n))

################################
MSE1 <- MSECal(b1, beta.true=beta)
MSE1p <- MSECal(b1p, beta.true=beta)
MSE1o <- MSECal(b1, yref=yo, Xnew=Xo, beta.true=beta)
MSE1po <- MSECal(b1p, yref=yo, Xnew=Xo, beta.true=beta)

# comrpasion old functions
MSE1 <- MSECal(b1, yref= yhat, beta.true=beta)
MSE1p <- MSECal(b1p, yref = yhat, beta.true=beta)
