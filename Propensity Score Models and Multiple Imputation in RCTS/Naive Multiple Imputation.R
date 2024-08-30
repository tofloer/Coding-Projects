


#####################################################################
### naive single imputation methods 
#####################################################################

rm(list=ls())

if (!require("rgl")) install.packages("rgl")
if (!require("VIM")) install.packages("VIM")

library(rgl)
library(VIM)

### Simulation data 1: 3 nv variables / 50% MCAR

#####################################################################
### MCAR Situation
#####################################################################

set.seed(1000)

n1 <- 2500
X1 <- rnorm(n1,8,3)
X2 <- 10 - 0.5*X1 + rnorm(n1,0,3)
X3 <- 5 + 0.6*X1 + 0.5*X2 + rnorm(n1,0,sqrt(2))

### create 50 % MCAR
misind <- sample(1:n1, n1/2)


data1 <- as.data.frame(cbind(X1,X2,X3))
is.na(data1$X3[misind]) <- TRUE
data1$X3[misind] <- NA
obsind <- which(!is.na(data1$X3))
## number of rows and columns
dim(data1)
head(data1)
### Visualization of missing data using the VIM package
barMiss(data1)
windows()
marginplot(data1[,-2])

######################################################################
### naive strategies
######################################################################
### (unconditional) mean imputation
data1.meanImp <- data1

### impute with mean(x3)
x3.bar <- mean(data1.meanImp$X3, na.rm=TRUE)
data1.meanImp$X3[misind] <- x3.bar

(mean.bef.Imp <- mean(X3))
(mean.after.Imp <- mean(data1.meanImp$X3))

(Var.bef.Imp <- var(X3))
(Var.after.Imp <- var(data1.meanImp$X3))

(corX1X3.bef.Imp <- cor(X1,X3))
(corX1X3.after.Imp <- cor(data1.meanImp$X1,data1.meanImp$X3))

summary(X3)
summary(data1.meanImp$X3)
### graphical confirmation
windows()
plot(density(X3), col = "blue", lwd=3, main = "kernel density X3",
		 xlim=range(X3),ylim=c(0,0.7))
par(new=TRUE)
plot(density(data1.meanImp$X3), col = "red", lwd=3, xlim=range(X3),ylim=c(0,0.7),
		 xlab="",ylab="",main="")
legend(x=15,y=0.25,legend = c("original","after imputation"),fill=c("blue","red"))

plot(density(X3), col = "blue", lwd=3,
		 main = "kernel density X3 -- second chance", xlim=range(X3),
		 ylim=c(0,30))
par(new=TRUE)
plot(density(data1.meanImp$X3), col = "red", lwd=3, xlim=range(X3),ylim=c(0,30),
		 xlab="",ylab="",main="")
legend(x=15,y=15,legend = c("original","after imputation"),fill=c("blue","red"))

### bivariate diagnostics look suspicious, too...
windows()
plot(X3~X1,col=rgb(0,0,1,0.3),pch=16, xlim = c(-5,20), ylim=c(4,21))
par(new=TRUE)
plot(X3~X1,data=data1.meanImp,col=rgb(1,0,0,0.3),pch=16, xlim = c(-5,20),
		 ylim = c(4,21), xlab ="", ylab= "")
legend(x=11,y=7,legend = c("original","after imputation"),fill=c("blue","red"))


################################################################################
data1.regImp <- data1

regmod <- lm(X3~X1+X2, data = data1.regImp[obsind,])

### Regression imputation
X3.HatMis <- predict(regmod, newdata=data1.regImp[misind,])

  
data1.regImp$X3[misind] <- X3.HatMis


(mean.bef.Imp <- mean(X3))
(mean.after.Imp <- mean(data1.regImp$X3))

(Var.bef.Imp <- var(X3))
(Var.after.Imp <- var(data1.regImp$X3))

(corX1X3.bef.Imp <- cor(X1,X3))
(corX1X3.after.Imp <- cor(data1.regImp$X1,data1.regImp$X3))

summary(X3)
summary(data1.regImp$X3)

windows()
plot(density(X3), col = "blue", lwd=3, main = "kernel density X3",
		 xlim=range(X3),ylim=c(0,0.2))
par(new=TRUE)
plot(density(data1.regImp$X3), col = "red", lwd=3, xlim=range(X3),
		 ylim=c(0,0.2),xlab="",ylab="",main="")
legend(x=16,y=0.15,legend = c("original","after imputation"),
       fill=c("blue","red"))


### bivariate plots do not tell us everything...
windows()
plot(X3~X1,col=rgb(0,0,1,0.3),pch=16, xlim = c(-5,20), ylim=c(4,21))
par(new=TRUE)
plot(X3~X1,data=data1.regImp,col=rgb(1,0,0,0.3),pch=16, xlim = c(-5,20),
		 ylim = c(4,21), xlab = "", ylab= "")
legend(x=12,y=7,legend = c("original","after imputation"),fill=c("blue","red"))


### Regression before and after imputation
regmodOrg <- lm(X3~X1+X2)
regmodImp <- lm(X3~X1+X2, data = data1.regImp)

summary(regmodOrg)
summary(regmodImp)


### 3D scatterplot

X1regImp <- data1.regImp$X1
X2regImp <- data1.regImp$X2
X3regImp <- data1.regImp$X3

plot3d(X3~X1+X2, col=rgb(0,0,1,0.1), xlim = c(-5,20),
											ylim=c(-6,20), zlim=c(4,21))
plot3d(X1regImp, X2regImp, X3regImp, col=rgb(1,0,0,0.2),add=TRUE)





#####################################################################
### naive single imputation methods part II
#####################################################################

rm(list=ls())

if (!require("rgl")) install.packages("rgl")
if (!require("VIM")) install.packages("VIM")

library(rgl)
library(VIM)

#####################################################################
### MCAR-Situation
#####################################################################

set.seed(1000)

n1 <- 2500
X1 <- rnorm(n1,8,3)
X2 <- 10 - 0.5*X1 + rnorm(n1,0,3)
X3 <- 5 + 0.6*X1 + 0.5*X2 + rnorm(n1,0,sqrt(2))

misind <- sample(1:n1,round(n1/2))
length(misind)

data1 <- as.data.frame(cbind(X1,X2,X3))
is.na(data1$X3[misind]) <- TRUE
obsind <- which(!is.na(data1$X3))
dim(data1)

### Visualizing missing data using the VIM package
help(package=VIM)
matrixplot(data1)
windows()
matrixplot(data1, sortby=1)

######################################################################
### further naive method
######################################################################

### simple hot deck imputation (resampling technique -- 
### NOT estimation-based)

data1.hotImp <- data1

X3.hdMis <- sample(X3[obsind],length(misind),replace=TRUE) 

data1.hotImp$X3[misind] <- X3.hdMis

### let's look at the estimates

(mean.bef.Imp <- mean(X3))
(mean.after.Imp <- mean(data1.hotImp$X3))
### fine as always

(Var.bef.Imp <- var(X3))
(Var.after.Imp <- var(data1.hotImp$X3))
### Variance has no systematic bias!

(KorrX1X3.bef.Imp <- cor(X1,X3))
(KorrX1X3.after.Imp <- cor(data1.hotImp$X1,data1.hotImp$X3))
### --> But correlations decreased!

summary(X3)
summary(data1.hotImp$X3)

windows()
plot(density(X3), col = "blue", lwd=3, main = "kernel density X3",
     xlim=range(X3),ylim=c(0,0.2))
par(new=TRUE)
plot(density(data1.hotImp$X3), col = "red", lwd=3, xlim=range(X3),
     ylim=c(0,0.2),xlab="",ylab="",main="")
legend(x=16,y=0.15,legend = c("original","after Imputation"),fill=c("blue","red"))
dev.off()
### no systematic differences for the kernel densities of X3

windows()
plot(X3~X1,col=rgb(0,0,1,0.3),pch=16, xlim = c(-5,20), ylim=c(4,21))
par(new=TRUE)
plot(X3~X1,data=data1.hotImp,col=rgb(1,0,0,0.3),pch=16, xlim = c(-5,20),
     ylim = c(4,21), xlab = "", ylab= "")
legend(x=12,y=7,legend = c("original","after Imputation"),fill=c("blue","red"))
dev.off()
### imputed values are circular -- not elliptic

### re-estimate the regression (before and after Imputation)
regmodOrg <- lm(X3~X1+X2)
regmodImp <- lm(X3~X1+X2, data = data1.hotImp)

summary(regmodOrg)
summary(regmodImp)

###==> Hot deck is "blind" towards associations between variables
### R^2 decreased!

################################################################################
### a (muuuch) less naive method
################################################################################
### Stochastic regression imputation

set.seed(2000)

n1 <- 2500
X1 <- rnorm(n1,8,3)
X2 <- 10 - 0.5*X1 + rnorm(n1,0,3)
X3 <- 5 + 0.6*X1 + 0.5*X2 + rnorm(n1,0,sqrt(2))

misind <- sample(1:n1,round(n1/2))
length(misind)

data1 <- as.data.frame(cbind(X1,X2,X3))
is.na(data1$X3[misind]) <- TRUE
obsind <- which(!is.na(data1$X3))


data1.cpdImp <- data1
regmod <- lm(X3~X1+X2, data = data1.cpdImp, subset = obsind)
sigma.hat <- summary(regmod)$sigma
beta.hat <- coef(regmod)
Design <- cbind(1,data1.cpdImp$X1,data1.cpdImp$X2)
X3.tilde <- rnorm(nrow(Design),Design%*%beta.hat,sigma.hat)
data1.cpdImp$X3[misind] <- X3.tilde[misind]


### let's see how far we'll get this time...

(mean.bef.Imp <- mean(X3))
(mean.after.Imp <- mean(data1.cpdImp$X3))
### good as usual...

(Var.bef.Imp <- var(X3))
(Var.after.Imp <- var(data1.cpdImp$X3))
### still doing fine!

(KorrX1X3.bef.Imp <- cor(X1,X3))
(KorrX1X3.after.Imp <- cor(data1.cpdImp$X1,data1.cpdImp$X3))
### --> still doing fine!

summary(X3)
summary(data1.cpdImp$X3)

windows()
plot(density(X3), col = "blue", lwd=3, main = "kernel density X3",
     xlim=range(X3),ylim=c(0,0.2))
par(new=TRUE)
plot(density(data1.cpdImp$X3), col = "red", lwd=3, xlim=range(X3),ylim=c(0,0.2),
     xlab="",ylab="",main="")
legend(x=16,y=0.15,legend = c("original","after Imputation"),fill=c("blue","red"))
dev.off()
### Kernel density almost identical as well...

windows()
plot(X3~X1,col=rgb(0,0,1,0.3),pch=16, xlim = c(-5,20), ylim=c(4,21))
par(new=TRUE)
plot(X3~X1,data=data1.cpdImp,col=rgb(1,0,0,0.3),pch=16, xlim = c(-5,20),
     ylim = c(4,21), xlab = "", ylab= "")
legend(x=12,y=7,legend = c("original","after Imputation"),fill=c("blue","red"))
dev.off()
### The bivariate scatterplot looks good

### Regression before and after Imputation
regmodOrg <- lm(X3~X1+X2)
regmodImp <- lm(X3~X1+X2, data = data1.cpdImp)

summary(regmodOrg)
summary(regmodImp)

### R^2 remains unchanged!!!

### The 3D scatterplot (the final frontier)

X1cpdImp <- data1.cpdImp$X1
X2cpdImp <- data1.cpdImp$X2
X3cpdImp <- data1.cpdImp$X3

plot3d(X3~X1+X2, col=rgb(0,0,1,0.1), xlim = c(-5,20),
       ylim=c(-6,20), zlim=c(4,21))
plot3d(X1cpdImp, X2cpdImp, X3cpdImp, col=rgb(1,0,0,0.2),add=TRUE)

### We did it! Unbiased imputations!
### Really???

###############################################################################
###
### Let's check our method based on a little simulation: We compare the
### before-deletion (BD) data with the remaining complete cases (CC) data 
### (before imputation) and the stochastic regression imputation data (cpd).
###
### Our quantities of interest now are E(X3) and the regression parameters.
### A diagnostics we will use the bias and the coverage (a=0.05).
### We repeat the simulation 1000 times (n is reduced from 2500 to 500 to
### speed up the computations -- 50 % MCAR is kept throughout the MC study).
### The true parameters can be derived theoretically: E(X2)=10-0.5*8=6 -->
### E(X3)=5+0.6*8+0.5*6=12.8;  E(alpha)=5; E(beta1)=0.6; E(beta2)=0.5
### 
###############################################################################

### We are saving results so we need a place for it (change 'path' at will)

path <- "D:/MI"
### example: path <- "c:/MI-Workshop"

if (!file.exists(file.path(path,"Data"))) dir.create(file.path(path,"Data"), 
                                                     recursive = TRUE)

set.seed(2000)

### the true values
trueVal <- c(12.8, 5, 0.6, 0.5)

modEst <- function(data, alpha=0.05, Wert=trueVal[-1],...) {
  model <- lm(X3~X1+X2, data=data, na.action = "na.omit")
  beta.hat <- model$coefficients
  se.beta.hat <- summary(model)$coefficients[ ,2]
  ## CI lower bound
  CI.low <- beta.hat - qnorm(1-alpha/2)*se.beta.hat
  ## CI upper bound
  CI.hi <- beta.hat + qnorm(1-alpha/2)*se.beta.hat
  contained <- as.numeric(CI.low <= Wert & CI.hi >= Wert)
  bias <- beta.hat-Wert
  Imp <- vector(length = 9)
  Imp[seq(1,by = 3,length.out = 3)] <- beta.hat
  Imp[seq(2,by = 3,length.out = 3)] <- contained
  Imp[seq(3,by = 3,length.out = 3)] <- bias
  return(Imp)
}

meanEst <- function(data, alpha=0.05, Wert=trueVal[1],...) {
  mean.X3 <- mean(data$X3, na.rm=TRUE)
  ## CI lower bound
  CI.low <- mean.X3 - qnorm(1-alpha/2)*sqrt(var(data$X3,
                                                na.rm=TRUE)/sum(!is.na(data$X3)))
  ## CI upper bound
  CI.hi <- mean.X3 + qnorm(1-alpha/2)*sqrt(var(data$X3,
                                               na.rm=TRUE)/sum(!is.na(data$X3)))
  contained <- as.numeric(CI.low <= Wert & CI.hi >= Wert)
  bias <- mean.X3-Wert
  Imp <- c(mean.X3, contained, bias)
  return(Imp)
}

nIter <- 1000
n <- 500
ImpBD <- ImpCC <- Impcpd <- matrix(nrow=nIter, ncol=12)

### how long does it take?
zeit1 <- Sys.time()

for (i in 1:nIter) {	
  X1 <- rnorm(n,8,3)
  X2 <- 10 - 0.5*X1 + rnorm(n,0,3)
  X3 <- 5 + 0.6*X1 + 0.5*X2 + rnorm(n,0,sqrt(2))
  data1 <- as.data.frame(cbind(X1,X2,X3))
  BD.dat <- data1 # BD (before deletion)
  misind <- sample(1:n,round(n/2))
  obsind <- which(!is.na(data1$X3))
  is.na(data1$X3[misind]) <- TRUE
  CC.dat <- data1 # CC (complete cases)
  cpd.dat <- data1
  regmod <- lm(X3~X1+X2, data = cpd.dat, subset = obsind)
  sigma.hat <- summary(regmod)$sigma
  beta.hat <- coef(regmod)
  Design <- cbind(1,cpd.dat$X1,cpd.dat$X2)
  X3.tilde <- rnorm(nrow(Design),Design%*%beta.hat,sigma.hat)
  cpd.dat$X3[misind] <- X3.tilde[misind] # stochastic regression imputation
  ##########################################################################
  BD.m <- meanEst(BD.dat)
  CC.m <- meanEst(CC.dat)
  cpd.m <- meanEst(cpd.dat)
  BD.r <- modEst(BD.dat)
  CC.r <- modEst(CC.dat)
  cpd.r <- modEst(cpd.dat)
  ImpBD[i, ] <- c(BD.m,BD.r)
  ImpCC[i, ] <- c(CC.m,CC.r)
  Impcpd[i, ] <- c(cpd.m,cpd.r)
}

save("meanEst", "modEst", "ImpBD","ImpCC","Impcpd",
     file = file.path(path,"Data/Simulation1.RData"))
zeit2 <- Sys.time()

cat("Overall run-time:",difftime(zeit2,zeit1,units="secs"),"seconds!\n")

### BD mean Coverage
BDcoverE3 <- factor(ImpBD[ ,2], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(BDcoverE3[1:i], main = "BD-Coverage E(X3)",col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=3)
}
dev.off()

### CC mean Coverage
CCcoverE3 <- factor(ImpCC[ ,2], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(CCcoverE3[1:i], main = "CC-Coverage E(X3)",col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=3)
}
dev.off()

### cpd mean Coverage

cpd.coverE3 <- factor(Impcpd[ ,2], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(cpd.coverE3[1:i], main = "cpd-Coverage E(X3)",col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=3)
}
dev.off()

### maybe a coincidence -- we have a look at the beta1:

### BD beta1 Coverage
BDcoverB1 <- factor(ImpBD[ ,8], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(BDcoverB1[1:i], main = expression(paste("BD-Coverage ",beta[1])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=3)
}
dev.off()

### CC beta1 Coverage
CCcoverB1 <- factor(ImpCC[ ,8], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(CCcoverB1[1:i], main = expression(paste("CC-Coverage ",beta[1])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=3)
}
dev.off()

### cpd beta1 Coverage
cpd.coverB1 <- factor(Impcpd[ ,8], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  Sys.sleep(0.001)
  plot(cpd.coverB1[1:i], main = expression(paste("cpd-Coverage ",beta[1])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=3)
}
dev.off()

### Tough luck...another fail!

### All results:
Overview <- rbind(colMeans(ImpBD),colMeans(ImpCC),colMeans(Impcpd))
colnames(Overview) <- c("E(x3)","C.E(X3)","B.E(X3)", "alpha", "C.alpha",
                        "B.alpha","beta1", "C.beta1", "B.beta1",
                        "beta2", "C.beta2", "B.beta2")
rownames(Overview) <- c("before deletion", "complete cases",
                        "cond. pred.-draws")

round(Overview, 4)




