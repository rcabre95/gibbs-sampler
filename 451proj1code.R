# pre-clean workspace
rm(list=ls())
graphics.off()

# Read in data
setwd("~/Documents/practice data/datasets 2/")
data0 <- read.delim("pigment.dat", header = T, sep = " ")

# Separate data by sample number
require(MCMCpack)
datas1 <- subset(data0, data0$Sample==1)
datas2 <- subset(data0, data0$Sample==2)
datas1
datas2

# begin initialization of variables before simulation
vara <- 86
varb <- 58
vare <- 1
num.its <- 10000
mu <- rep(0,num.its)
a <- rep(0,num.its*15)
dim(a) <- c(num.its,15)
b <- rep(0,num.its*15)
dim(b) <- c(num.its,15)
y <- rep(0,num.its*15)
dim(y) <- c(num.its,15)
set.seed(24)
e <- rnorm(15, mean=0, sd=sqrt(vare))
n <- 30
I <- as.data.frame(c(1:15))
J <- 2
yi. <- rep(0,15)
for(i in 1:15){
  yi.[i] <- mean(subset(data0[,3], data0[,1]==i))
}
y.. <- mean(data0[,3])
v1 <- ((2/vare)+(1/vara))^(-1)
v2 <- ((1/vare)+(1/varb))^(-1)

set.seed(24)
a[1,] <- rnorm(15, mean=0, sd=sqrt(vara))
b[1,] <- rnorm(15, mean=0, sd=sqrt(varb))
mu[1] <- rnorm(1,mean=mean(data0[,3]),sd=sd(data0[,3]))

y[1,] <- mu[i]+a[i,]+b[i,]+e-mu[i-1]-a[i-1,]


# initiate gibbs sampling algorithm
set.seed(24)
for(i in 2:num.its){
  mu[i] <- rnorm(1,mean= ((mu[i]+a[i,]+b[i,]+e-mu[i-1]-a[i-1,])/n)-(2*sum(a[i-1,])/n)-((sum(b[i-1,]))/n), sd=sqrt(vare/n))
  for(j in 1:15){
    a[i,j] <- rnorm(1,mean=(J*v1/vare)*((yi.[j])-mu[i]-(1/J)*sum(b[i-1])),sd= sqrt(v1))
    b[i,] <- rnorm(1,mean=(v2/vare)*(mu[i]+a[i,]+b[i,]+e-mu[i]-a[i,j]),sd=sqrt(v2))
  }
  y[i,] <- mu[i]+a[i,]+b[i,]+e-mu[i-1]-a[i-1,]
}

#plot the convergence of each parameter
## View plots
par(mfrow=c(2,2))
plot(a[,2],type="l",ylab="Alpha",xlab="Iteration")
plot(b[,3],type="l",ylab="Beta",xlab="Iteration")
plot(y[,1],type="l",ylab="Y",xlab="Iteration")
plot(mu, type="l",ylab="Mu",xlab="Iteration")
mtext("Parameter Sample Paths", outer = TRUE, cex = 1.5, line = -2)

## Produce plots
png("sample_path.png")
par(mfrow=c(2,2))
plot(a[,2],type="l",ylab="Alpha",xlab="Iteration")
plot(b[,3],type="l",ylab="Beta",xlab="Iteration")
plot(y[,1],type="l",ylab="Y",xlab="Iteration")
plot(mu, type="l",ylab="Mu",xlab="Iteration")
mtext("Parameter Sample Paths", outer = TRUE, cex = 1.5, line = -2)
dev.off()

#Cumulative Sums (must fluctuate around 0)
# implement burn in for faster convergence
burn.in <- 1:1000
summary.a <- summary(a[-burn.in])
summary.b <- summary(b[-burn.in,1])
summary.mu <- summary(mu[-burn.in])
summary.y <- summary(y[-burn.in,1])
round(rbind(summary.a,summary.b,summary.mu,summary.y),2)

# mu burn in
temp <- mu[-burn.in]
ntemp <- length(temp)      
mtemp <- mean(temp)
temp <- cumsum(temp-mtemp)
plot(temp, type="l", xlab="iteration", ylab="mu")

# alpha burn in
temp <- a[-burn.in,1]
ntemp <- length(temp)      
mtemp <- mean(temp)
temp <- cumsum(temp-mtemp)
plot(temp, type="l", xlab="iteration", ylab="alpha")

# beta burn in 
temp <- b[-burn.in,1]
ntemp <- length(temp)      
mtemp <- mean(temp)
temp <- cumsum(temp-mtemp)
plot(temp, type="l", xlab="iteration", ylab="beta")

# predictor burn in
temp <- y[-burn.in,1]
ntemp <- length(temp)      
mtemp <- mean(temp)
temp <- cumsum(temp-mtemp)
plot(temp, type="l", xlab="iteration", ylab="y")
mtext("Cumulative Sum Plots", outer = TRUE, cex = 1.5, line = -2)

#acf plots (must converge to zero)
acf(a[-burn.in,1], lag.max=15, type="correlation", xlab="alpha");
acf(b[-burn.in,1], lag.max=15, type="correlation", xlab="beta");
acf(mu[-burn.in], lag.max=15, type="correlation", xlab="mu");
acf(y[-burn.in,1], lag.max=15, type="correlation", xlab="y");


##### reparametrization #####

## initialization of variables
data0 <- read.delim("pigment.dat", header = T, sep = " ")
datas1 <- subset(data0, data0$Sample==1)
datas2 <- subset(data0, data0$Sample==2)
vara <- 86
varb <- 58
vare <- 1
num.its <- 10000
v1 <- ((2/vare)+(1/vara))^(-1)
v2 <- ((1/vare)+(1/varb))^(-1)
v3 <- 1/((2/varb)/(1/vara))
mu <- rep(0,num.its)
y <- rep(0,num.its*3)
dim(y) <- c(num.its,3)
n <- rep(0,num.its*3)
dim(n) <- c(num.its,3)


set.seed(24)
mu[1] <- rnorm(1,mean=mean(data0[,3]),sd=sd(data0[,3]))
e <- rnorm(3,mean=0,sd=sqrt(vare))
y[1,] <- rnorm(3,mean=mu[1],sd=sqrt(vara))
n[1,] <- rnorm(3,mean=y[1,],sd=sqrt(varb))

## run reparameterized gibbs sampler
for(i in 2:num.its){
  mu[i] <- rnorm(1,(1/15)*sum(y[i-1,]),(1/15)*vara)
  y[i,] <- rnorm(3,mean=v3*((1/varb)*(sum(n[i-1,])+(mu[i]/vara))),sd=sqrt(v3))
  n[i,] <- rnorm(3,mean=v2*((mu[i]+a[i,]+b[i,]+e)/vare)+(y[i,]/varb),sd=sqrt(v2))
}

## plot sample path
plot(mu,type="l",ylab="Mu",xlab="Iteration")
plot(y[,2],type="l",ylab="Y",xlab="Iteration")
plot(n[,2],type="l",ylab="nj",xlab="Iteration")

#plot cumulative sum plots
burn.in <- 1:1000
summary.mu <- summary(mu[-burn.in])
summary.y <- summary(y[-burn.in,1])
summary.n <- summary(n[-burn.in])
round(rbind(summary.n,summary.mu,summary.y),2)

temp <- mu[-burn.in]
ntemp <- length(temp)      
mtemp <- mean(temp)
temp <- cumsum(temp-mtemp)
plot(temp, type="l", xlab="iteration", ylab="mu")

temp <- y[-burn.in,1]
ntemp <- length(temp)      
mtemp <- mean(temp)
temp <- cumsum(temp-mtemp)
plot(temp, type="l", xlab="iteration", ylab="y")

temp <- n[-burn.in,1]
ntemp <- length(temp)      
mtemp <- mean(temp)
temp <- cumsum(temp-mtemp)
plot(temp, type="l", xlab="iteration", ylab="n")

# autocorrelation plots
acf(y[-burn.in,1], lag.max=20, type="correlation", xlab="y");
acf(n[-burn.in,1], lag.max=20, type="correlation", xlab="N");
acf(mu[-burn.in], lag.max=20, type="correlation", xlab="Mu");

# Correlation matrix
cc <- cor(cbind(a,mu,b,y,n))
View(cc)

# we can see that the reparametrization of the distribution was
# mostly beneficial