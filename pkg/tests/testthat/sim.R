rm(list=ls())
set.seed(12345)
# starting values
a0 <- 1 # between 0 and 1
b0 <- 50 #runif(1) + 1
m <- 100000

a <- a1 <- vector(mode="numeric", length=m)
b <- b1 <- vector(mode="numeric", length=m)
a[1] <- a1[1] <- a0
b[1] <- b1[1] <- b0

for (i in 2:m) {
  a[i] =  a[i-1]*(1-a[i-1]/b[i-1]^2)
  #a1[i] =  a1[i-1]*(1-a1[i-1]/b[i-1]^2) - runif(1)/100
  b[i] = sqrt(a[i-1]-a[i]) + b[i-1]
  #b1[i] = sqrt(a1[i-1]-a1[i]) + b1[i-1]
}

# cut off
c <- 10000
reg <- lm(log(a[c:m]) ~ log(c(c:m)))
reg
reg <- lm(log(b[c:m]) ~ log(c(c:m)))
reg

#############

d <- a/b^2
s <- d*(1:m)
plot(s)
