# DGP
n <- 100
p <- 100
s <- 3 # error when s=1, p=100, n=100: sigma=0; to fix
n1 <- 50
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(rep(1,s), rep(0,p-s))
y <- X%*%beta + rnorm(n)
y_det <-  X%*%beta
#
Xo <- matrix(rnorm(n1*p), ncol=p)
fo <-  Xo%*%beta
yo <- fo + rnorm(n1)

yhat_oracle <- predict(lm(y~ -1 + X[,1:s]))
beta_oracle <- coef(lm(y~ -1 + X[,1:s]))
# Estimation with L2Boosting

b1 <- L2Boost(X,y,iter=200)
b1_det <- L2Boost(X,y_det,iter=1000)

b1p <- L2Boost(X,y,iter=200, post=TRUE)
b2 <- L2Boost(X,y,iter=200, beta.start=rep(1,p))
b3 <- L2BoostOGA(X,y,iter=200)
b3 <- L2BoostOGA(X,y,iter=200)
b4 <- L1GBoost(X,y, iter=200,n_per=10,weight_p=rep(1,n))

Xnew <- matrix(rnorm(20*p), ncol=p)
ynew <-  Xnew%*%beta + rnorm(20)
# prediction
yhat1 <- predict(b1)
yhat1new <- predict(b1, newx=Xnew)
yhat1p <- predict(b1p)
yhat1pnew <- predict(b1p, newx=Xnew)
yhat3 <- predict(b3)
yhat3new <- predict(b3, newx=Xnew)

################################
MSE1 <- MSECal(b1, beta.true=beta)
MSE1p <- MSECal(b1p, beta.true=beta)
MSE1o <- MSECal(b1, yref=yo, Xnew=Xo, beta.true=beta)
MSE1po <- MSECal(b1p, yref=yo, Xnew=Xo, beta.true=beta)

# comrpasion old functions
MSE1 <- MSECal(b1, yref= yhat, beta.true=beta)
MSE1p <- MSECal(b1p, yref = yhat, beta.true=beta)
