# R-Code for the Poster "Non-Uniform Random Number Generation" by Tobias Flörchinger

###############################################################################
##############################   Load Packages   ##############################
###############################################################################

### Packages ###
packages =  c("Runuran"   #Variate generation (esp. RoU and TDR)
             ,"goftest"   #Anderson-Darling Test
             ,"zoo"       #Test Construction
             ,"EnvStats"  #Triangle Dirstibution
              )

### Install package if not available ###
inst <- packages %in% installed.packages()    
if(length(packages[!inst]) > 0) install.packages(packages[!inst])

### Load Packages ###
loading <- lapply(packages, require, character.only=TRUE)



##############################################################################################
##############################   Functions for the Simulation   ##############################
##############################################################################################

### CHI-Squared Test ###
chi.test <- function(x,px,shape1,shape2) {
  histo <- hist(x, breaks=c(quantile(x, probs = seq(0,1,0.25))), plot = FALSE) # Generating classes by abusing the histogram function
  breaks_cdf <- px(histo$breaks,shape1,shape2)                                 # Calculating the PDF at every break point to the next class
  null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])                # Theoretical probabilities for each class
  return(chisq.test(histo$counts, p=null.probs, rescale.p = T))
  }

### Runs Test ###
runs.test <- function(x) {
  m <- median(x)                                # serves as a threshold to divide the data
  x <- x[x!=median(x)]                          # ensuring an even number of -1s and 1s by deleting median values (if there are any)
  s <- sign(x-m)                                # x-m will be negative for any value smaller than the median
  n1 <- length(s[s==1])                         # number of 1s (values larger than the median)
  n2 <- length(s[s==-1])                        # number of -1s (values smaller than the median)
  runs <- rle(s)                                # recording the type and lengths of all runs
  r1 <- length(runs$lengths[runs$values==1])    # number of runs with 1s
  r2 <- length(runs$lengths[runs$values==-1])   # number of runs with -1s  
  n <- n1+n2                                    # total sample size
  # We now assume that the number of runs is approximately normal distributed.
  # This is reasonable for n1, n2 > 20.
  mu <- (2*n1*n2)/n + 1                         # expectation for the number of runs
  sigma2 <- (2*n1*n2*(2*n1*n2 - n))/(n^2*(n-1)) # variance for the number of runs
  r <- r1+r2
  p0 <- pnorm((r-mu)/sqrt(sigma2))              # Test statistic put into the PDF of the normal
  p <- 2*min(p0,1-p0)                           # p-value for a two-sided hypothesis test
  return(p)
  }

### Composition of KS, AD, CHI^2, RUNS Tests for continuous variate methods and target distributions ###
all.tests <- function(Method,iter) {   
  #... Norm(0,1)
  p.value <- ks.test(norm,"pnorm")$p.value
  Method[iter,2] <- p.value
  Method[iter,3] <- ifelse(p.value < 0.05,0,1)
  p.value <- ad.test(norm,"pnorm")$p.value
  Method[iter,4] <- p.value
  Method[iter,5] <- ifelse(p.value < 0.05,0,1)
  p.value <- chi.test(norm,pnorm, 0, 1)$p.value
  Method[iter,6] <- p.value
  Method[iter,7] <- ifelse(p.value < 0.05,0,1)
  Method[iter,8] <- runs.test(norm)
  #... Exp(1)
  p.value <- ks.test(exp,"pexp")$p.value
  Method[iter,10] <- p.value
  Method[iter,11] <- ifelse(p.value < 0.05,0,1)
  p.value <- ad.test(exp,"pexp")$p.value
  Method[iter,12] <- p.value
  Method[iter,13] <- ifelse(p.value < 0.05,0,1)
  p.value <- chi.test(exp,pexp, 1, 1)$p.value
  Method[iter,14] <- p.value
  Method[iter,15] <- ifelse(p.value < 0.05,0,1)
  Method[iter,16] <- runs.test(exp)
  #... Gamma(2,1)
  p.value <- ks.test(gamma,"pgamma", 2, 1)$p.value
  Method[iter,18] <- p.value
  Method[iter,19] <- ifelse(p.value < 0.05,0,1)
  p.value <- ad.test(gamma,"pgamma", 2, 1)$p.value
  Method[iter,20] <- p.value
  Method[iter,21] <- ifelse(p.value < 0.05,0,1)
  p.value <- chi.test(gamma,pgamma, 2, 1)$p.value
  Method[iter,22] <- p.value
  Method[iter,23] <- ifelse(p.value < 0.05,0,1)
  Method[iter,24] <- runs.test(gamma)
  return(Method)
  }



############################################################################
##############################   Simulation   ##############################
############################################################################
set.seed(123)
T = 1000
N = 10000   # or 1000

# Interim matrices for storing results(run-times,tests) for all iterations T and methods
NINV <- TNS.c <- R.c <- TDR <- ROU <- matrix(nrow=T,ncol=24)
TNS.d <- R.d <- ALS <- matrix(nrow=T,ncol=8)
TRI <- matrix(nrow=T,ncol=16)

strt<-Sys.time()
for(t in 1:T) { 
  if(t %% 25 == 0) cat("Iteration", t, "of", T,"...\n")
  
  ####################   Continuous Distributions   ######################################################################
    
  #######################################
  ### Numerical Inversion (Bisection) ###
  #######################################

  ### Norm(0,1) ###
  method<-Sys.time()
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  F <- function(x) (1/2)*(1+erf(x/sqrt(2)))     # cdf is directly determined
  F.inv <- function(y){uniroot(function(x){F(x)-y},interval=c(-25,25),extendInt = "yes")$root}
  F.inv <- Vectorize(F.inv)                     # c(-25,25) was faster than c(0,25)
  u <- runif(N)
  norm <- F.inv(u)                              # uniroot correspondes to bisection method that minimizes F(x)-u wir u~Unif(0,1) by finding the root x
  NINV[t,1] <- as.numeric(Sys.time()-method)
  
  ### Exp(1) ###
  method<-Sys.time()
  f <- function(x) exp(-x)
  F <- function(x) {integrate(f,0,x)$value}     # cdf follows from integration
  F.inv <- function(y){uniroot(function(x){F(x)-y},interval=c(0,25),extendInt = "yes")$root}
  F.inv <- Vectorize(F.inv)
  u <- runif(N)
  exp <- F.inv(u) 
  NINV[t,9] <- as.numeric(Sys.time()-method)
  
  ### Gamma(2,1) ###
  method<-Sys.time()
  f <- function(x) (1/gamma(2))*x*exp(-x)
  F <- function(x) {integrate(f,0,x)$value}     # cdf follows from integration
  F.inv <- function(y){uniroot(function(x){F(x)-y},interval=c(0,25),extendInt = "yes")$root}
  F.inv <- Vectorize(F.inv)
  u <- runif(N)
  gamma <- F.inv(u)
  NINV[t,17] <- as.numeric(Sys.time()-method)
  
  # KS,AD,CHI,RUNS
  NINV <- all.tests(NINV,t)
  
  ######################################
  ### Transformation and Properties  ###
  ######################################
  
  # Did not make it to the Poster as it is no universal method but rather uses transformations to get
  # closed form representation using polar coordinates, inversion, convolution and further properties.
  
   ### Norm(0,1) ###
    method<-Sys.time()
    u1 <- runif(N)
    u2 <- runif(N)
    norm <- sqrt(-2*log(u1))*cos(2*pi*u2)
    TNS.c[t,1] <- as.numeric(Sys.time()-method)
    
    ### Exp(1) ###
    method<-Sys.time()
    exp <- vector(length = N)
    u <- runif(N)
    for(i in 1:length(u1)) { exp[i] = -log(1-u[i]) }
    TNS.c[t,9] <- as.numeric(Sys.time()-method)
    
    ### Gamma(2,1) ###
    method<-Sys.time()
    gamma <- vector(length = N)
    u1 <- runif(N)
    u2 <- runif(N)
    for(i in 1:length(u1)) { gamma[i] = -log(1-u1[i])+-log(1-u2[i]) }
    TNS.c[t,17] <- as.numeric(Sys.time()-method)
    
    # KS,AD,CHI,RUNS
    TNS.c <- all.tests(TNS.c,t)
  
  ###########################
  ### Composition Method  ###
  ###########################
  
  ### Triangle(-1,1,0) ###
  ctri <- numeric(length = N)
  method<-Sys.time()
  u1 <- runif(N)
  u2 <- runif(N)
  for(i in 1:N){
    if(u1[i] < 0.5){
      ctri[i]<- sqrt(u2[i])-1
    }else{ctri[i]<- 1-sqrt(1-u2[i])}
  }
  TRI[t,9] <- as.numeric(Sys.time()-method)
  
  TRI[t,10] <- ks.test(ctri, ptri, min=-1, max=1, mode=0)$p.value
  TRI[t,11] <- ifelse(TRI[t,10] < 0.05,0,1)
  TRI[t,12] <- ad.test(ctri, ptri, min=-1, max=1, mode=0)$p.value
  TRI[t,13] <- ifelse(TRI[t,12] < 0.05,0,1)
  histo <- hist(ctri, breaks=c(quantile(ctri, probs = seq(0,1,0.25))), plot = FALSE)
  null.probs <- rollapply(ptri(histo$breaks,-1, 1, 0), 2, function(x) x[2]-x[1])
  TRI[t,14] <- chisq.test(histo$counts, p=null.probs, rescale.p = T)$p.value
  TRI[t,15] <- ifelse(TRI[t,14] < 0.05,0,1)
  TRI[t,16] <- runs.test(ctri)
  
  #########################
  ### R Default Methods ###
  #########################
  
  ### Norm(0,1) ###
  method<-Sys.time()
  norm <- rnorm(N, 0, 1)
  R.c[t,1] <- as.numeric(Sys.time()-method)
  
  ### Exp(1) ###
  method<-Sys.time()
  exp <- rexp(N, 1)
  R.c[t,9] <- as.numeric(Sys.time()-method)

  ### Gamma(2,1) ###
  method<-Sys.time()
  gamma <- rgamma(N, 2, 1)
  R.c[t,17] <- as.numeric(Sys.time()-method)

  # KS,AD,CHI,RUNS
  R.c <- all.tests(R.c,t)
  
  ### Triangle(-1,1,0) ###
  # R-Pakage EnvStats
  method<-Sys.time()
  triangle <- rtri(N, -1,1,0)
  TRI[t,1] <- as.numeric(Sys.time()-method)
  
  TRI[t,2] <- ks.test(triangle, ptri, min=-1, max=1, mode=0)$p.value
  TRI[t,3] <- ifelse(TRI[t,2] < 0.05,0,1)
  TRI[t,4] <- ad.test(triangle, ptri, min=-1, max=1, mode=0)$p.value
  TRI[t,5] <- ifelse(TRI[t,4] < 0.05,0,1)
  histo <- hist(triangle, breaks=c(quantile(triangle, probs = seq(0,1,0.25))), plot = FALSE) 
  null.probs <- rollapply(ptri(histo$breaks,-1, 1, 0), 2, function(x) x[2]-x[1])
  TRI[t,6] <- chisq.test(histo$counts, p=null.probs, rescale.p = T)$p.value
  TRI[t,7] <- ifelse(TRI[t,6] < 0.05,0,1)
  TRI[t,8] <- runs.test(triangle)
  
  ###########################################
  ### Transformed Density Rejection (TDR) ###
  ###########################################
  
  ### Norm(0,1) ###
  f <- function (x) { (1/sqrt(2*pi))*exp(-0.5*x^2) }
  method<-Sys.time()
  norm <- ur(tdr.new(pdf=f, dpdf=NULL, lb=-Inf, ub=Inf, islog=FALSE),N)
  TDR[t,1] <- as.numeric(Sys.time()-method)

  ### Exp(1) ###
  e <- function (x) { exp(-x) }
  method<-Sys.time()
  exp <- ur(tdr.new(pdf=e, dpdf=NULL, lb=0, ub=Inf, islog=FALSE),N)
  TDR[t,9] <- as.numeric(Sys.time()-method)

  ### Gamma(2,1) ###
  g <- function (x) { (1/gamma(2))*x*exp(-x) }
  method<-Sys.time()
  gamma <- ur(tdr.new(pdf=g, dpdf=NULL,lb = 0, ub = Inf, islog=FALSE),N)
  TDR[t,17] <- as.numeric(Sys.time()-method)

  # KS,AD,CHI,RUNS
  TDR <- all.tests(TDR,t)
  
  ##############################################
  ### Simple Ratio-Of-Uniforms Method (SROU) ###
  ##############################################
  
  ### Norm(0,1) ###
  f <- function (x) { (1/sqrt(2*pi))*exp(-0.5*x^2) }
  method<-Sys.time()
  norm <- ur(srou.new(pdf=f, lb=-Inf, ub=Inf, islog=FALSE, mode=0, area=integrate(f,-Inf,Inf)$value),N)
  ROU[t,1] <- as.numeric(Sys.time()-method)
  
  ### Exp(1) ###
  e <- function (x) { exp(-x) }
  method<-Sys.time()
  exp <- ur(srou.new(pdf=e, lb=0, ub=Inf, islog=FALSE, mode=0, area=integrate(e,0,Inf)$value),N)
  ROU[t,9] <- as.numeric(Sys.time()-method)
  
  ### Gamma(2,1) ###
  g <- function (x) { (1/gamma(2))*x*exp(-x) }
  method<-Sys.time()
  gamma <- ur(srou.new(pdf=g, lb=0, ub=Inf, islog=FALSE, mode=1, area=integrate(g,0,Inf)$value),N)
  ROU[t,17] <- as.numeric(Sys.time()-method)

  # KS,AD,CHI,RUNS
  ROU <- all.tests(ROU,t)
  
  
  ####################   Discrete Distributions   ######################################################################
  
  ######################################
  ### Transformation and Properties  ###
  ######################################
  
  # Did not make it to the Poster as it is no universal method but rather uses transformations to get
  # closed form representation using polar coordinates, inversion, convolution and further properties.
  
    ### Bin(40,0.5) ###
    binom <- vector(length = N)
    method<-Sys.time()
    for (i in 1:N) {
      u <- runif(40)
      binom[i] = sum(u > 1-0.5) 
    }
    TNS.d[t,1] <- as.numeric(Sys.time()-method) 
    TNS.d[t,2] <- chi.test(binom,pbinom,40,0.5)$p.value
    TNS.d[t,3] <- ifelse(TNS.d[t,2] < 0.05,0,1)
    TNS.d[t,4] <- runs.test(binom)
   
    ### Pois(10) ###
    method<-Sys.time()
    pois <- sapply(1:N, function(i){
      u<-runif(100)
      x<--log(1-u)/10
      y=cumsum(x)
      length(which(y<=1))
    })
    TNS.d[t,5] <- as.numeric(Sys.time()-method) 
    histo <- hist(pois, breaks=c(quantile(pois, probs = seq(0,1,0.25))), plot = FALSE)
    null.probs <- rollapply(ppois(histo$breaks,10), 2, function(x) x[2]-x[1])
    TNS.d[t,6] <- chisq.test(histo$counts, p=null.probs, rescale.p = T)$p.value
    TNS.d[t,7] <- ifelse(TNS.d[t,6] < 0.05,0,1)
    TNS.d[t,8] <- runs.test(pois)
  
  #########################
  ### R Default Methods ###
  #########################
  
  ### Bin(40,0.5) ###
  method<-Sys.time()
  binom <- rbinom(N,40,0.5)
  R.d[t,1] <- as.numeric(Sys.time()-method)
  
  R.d[t,2] <-  chi.test(binom,pbinom,40,0.5)$p.value
  R.d[t,3] <- ifelse(R.d[t,2] < 0.05,0,1)
  R.d[t,4] <-  runs.test(binom)
  
  ### Pois(10) ###
  method<-Sys.time()
  pois <- rpois(N, 10)
  R.d[t,5] <- as.numeric(Sys.time()-method)
  
  histo <- hist(pois, breaks=c(quantile(pois, probs = seq(0,1,0.25))), plot = FALSE)
  null.probs <- rollapply(ppois(histo$breaks,10), 2, function(x) x[2]-x[1])
  R.d[t,6] <- chisq.test(histo$counts, p=null.probs, rescale.p = T)$p.value
  R.d[t,7] <- ifelse(R.d[t,6] < 0.05,0,1)
  R.d[t,8] <- runs.test(pois)
  
  ####################
  ### Alias Method ###
  ####################
  
  ### Bin(40,0.5) ###
  method<-Sys.time()
  binom <- ur(daud.new(udbinom(size=40, prob = 0.5)),N)
  ALS[t,1] <- as.numeric(Sys.time()-method)
  
  ALS[t,2] <- chi.test(binom,pbinom,40,0.5)$p.value
  ALS[t,3] <- ifelse(ALS[t,2] < 0.05,0,1)
  ALS[t,4] <- runs.test(binom)
  
  ### Pois(10) ###
  method<-Sys.time()
  pois <- ur(dau.new(pv=dpois(0:N, 10), from=1),N)
  ALS[t,5] <- as.numeric(Sys.time()-method) 

  histo <- hist(pois, breaks=c(quantile(pois, probs = seq(0,1,0.25))), plot = FALSE)
  null.probs <- rollapply(ppois(histo$breaks,10), 2, function(x) x[2]-x[1])
  ALS[t,6] <- chisq.test(histo$counts, p=null.probs, rescale.p = T)$p.value
  ALS[t,7] <- ifelse(ALS[t,6] < 0.05,0,1)
  ALS[t,8] <- runs.test(pois)
  
}
print(Sys.time()-strt)



##############################################################################
##############################  Results Tables  ##############################
##############################################################################

### Result tables (continuous)
Results <- matrix(nrow=3, ncol=40)
for (i in 1:8){
  Results[1,i] <- colMeans(NINV)[i]
  Results[1,i+8] <- colMeans(TNS.c)[i]
  Results[1,i+16] <- colMeans(R.c)[i]
  Results[1,i+24] <- colMeans(TDR)[i]
  Results[1,i+32] <- colMeans(ROU)[i]
  Results[2,i] <- colMeans(NINV)[i+8]
  Results[2,i+8] <- colMeans(TNS.c)[i+8]
  Results[2,i+16] <- colMeans(R.c)[i+8]
  Results[2,i+24] <- colMeans(TDR)[i+8]
  Results[2,i+32] <- colMeans(ROU)[i+8]
  Results[3,i] <- colMeans(NINV)[i+16]
  Results[3,i+8] <- colMeans(TNS.c)[i+16]
  Results[3,i+16] <- colMeans(R.c)[i+16]
  Results[3,i+24] <- colMeans(TDR)[i+16]
  Results[3,i+32] <- colMeans(ROU)[i+16]
}
Results <- round(Results, digits = 5) #round
rownames(Results) <- c("N(0,1)","Exp(1)","Gamma(2,1)")
colnames(Results) <- c("NINV.Time","NINV.ks.E[p]","NINV.ks.C","NINV.ad.E[p]","NINV.ad.C","NINV.chi.E[p]","NINV.chi.C","NINV.Run.Test",   #.C corresonds to non rejection rates
                       "TNS.Time","TNS.ks.E[p]","TNS.ks.C","TNS.ad.E[p]","TNS.ad.C","TNS.chi.E[p]","TNS.chi.C","TNS.Run.Test",
                       "R.Time","R.ks.E[p]","R.ks.C","R.ad.E[p]","R.ad.C","R.chi.E[p]","R.chi.C","R.Run.Test",
                       "TDR.Time","TDR.ks.E[p]","TDR.ks.C","TDR.ad.E[p]","TDR.ad.C","TDR.chi.E[p]","TDR.chi.C","TDR.Run.Test",
                       "ROU.Time","ROU.ks.E[p]","ROU.ks.C","ROU.ad.E[p]","ROU.ad.C","ROU.chi.E[p]","ROU.chi.C","ROU.Run.Test")
Results.Time <- Results[,c(1,9,17,25,33)]
Results.ks <- Results[,c(2,10,18,26,34,3,11,19,27,35)]
Results.ad <- Results[,c(4,12,20,28,36,5,13,21,29,37)]
Results.chi <- Results[,c(6,14,22,30,38,7,15,23,31,39)]
Results.run <- Results[,c(8,16,24,32,40)]
print(Results.Time)               #Runtimes
print(Results.ks)                 #KS-Tests
print(Results.ad)                 #AD-Tests
print(Results.chi)                #Chi^2-Tests
print(Results.run)                #Runs Tests

### Result table TRI 
Results.TRI <- rbind(colMeans(TRI)[1:8],colMeans(TRI)[9:16])
rownames(Results.TRI) <- c("NINV","Comp")
colnames(Results.TRI) <- c("Time","ks.E[p]","ks.C","ad.E[p]","ad.C","chi.E[p]","chi.C","Run.Test")   #.C corresonds to non rejection rates
print(Results.TRI)

### Result table (discrete)
Results.d <- matrix(nrow=2, ncol=12)
for (i in 1:4){
  Results.d[1,i] <- colMeans(TNS.d)[i]
  Results.d[1,i+4] <- colMeans(R.d)[i]
  Results.d[1,i+8] <- colMeans(ALS)[i]
  Results.d[2,i] <- colMeans(TNS.d)[i+4]
  Results.d[2,i+4] <- colMeans(R.d)[i+4]
  Results.d[2,i+8] <- colMeans(ALS)[i+4]
}
Results.d <- round(Results.d, digits = 5) #round
rownames(Results.d) <- c("Bin","Pois")
colnames(Results.d) <- c("TNS.Time","TNS.chi.E[p]","TNS.chi.C","TNS.Run.Test",    #.C corresonds to non rejection rates
                       "R.Time","R.chi.E[p]","R.chi.Cover","R.Run.Test",
                       "ALS.Time","ALS.chi.E[p]","ALS.chi.C","ALS.Run.Test")
print(Results.d)                 #All results (discrete)



#######################################################################
############################## Graphics  ##############################
#######################################################################

############
### Data ###
############

### Runtimes: Extract N(0,1), Bin(40,0.5) und Triangle(-1,1,0) ###
TimeN1 <- c(R.c[,1],NINV[,1], TDR[,1], ROU[,1])
TimeB1 <- c(R.d[,1], ALS[,1])
TimeT1 <- c(TRI[,1], TRI[,9])
time.all1 <- c(TimeN1, TimeB1, TimeT1)

### Extract Run-Tests ###
runn1 <- data.frame(cbind(R.c[,8],NINV[,8], TDR[,8], ROU[,8]))
runb1 <- data.frame(cbind(R.d[,4], ALS[,4]))
runt1 <- data.frame(cbind(TRI[,8], TRI[,16]))


#############
### Plots ###
#############

# NOTE: Graphics may only be visible in the created PNG file for which they are scaled and may look strange in the R window

### Runtimes ###
y1 <-c(rep(1, 1000), rep(2, 1000), rep(3,1000), rep(4,1000),
       rep(5, 1000), rep(6, 1000), rep(7,1000), rep(8,1000))

data.time.n <- data.frame(cbind(time.all1, y1))

png(filename= "Time.png", width = 800, height= 800)
plot(data.time.n, col= c(rep("red", 1000), rep("blue", 1000), rep("seagreen", 1000), rep("orange", 1000),
                         rep("red", 1000), rep("sienna", 1000), rep("blue", 1000), rep("darkcyan", 1000)),
     pch=4, xlab= "run time in sec", ylab= "RVG", cex=3, cex.axis=2.5, cex.lab= 2)
abline(h=4.5, lwd = 5)
abline(h=6.5, lwd = 5)
legend("topright", 
       legend = c("Composition", "Alias", "Ratio-of-Uniforms", "Transformed Density Rejection", "num.Inversion", "implemented in R"),
       col= c("darkcyan", "sienna", "orange", "seagreen", "blue", "red"),
       pch=4, cex = 2.5)
dev.off()


### Run-Tests ###
dat.run1 <- data.frame(cbind(runn1, runb1, runt1))
colnames(dat.run1) <- c("R", "NINV", "TDR", "ROU", "R", "Alias", "NINV", "Comp")

png(filename= "Run.png",width = 900, height= 800)
boxplot(dat.run1, col= c("red","blue", "seagreen", "orange", "red", "sienna", "blue", "darkcyan"),
        ylab= "p-values of run test", xlab= "RVG", cex.axis=2.3, cex.lab= 2,
        mar= c(5,4,4,2)+0.5)
abline(v= 4.5, lwd=5)
abline(v= 6.5, lwd=5)
abline(h=0.05, col= "darkgrey", lwd=4, lty=4)
dev.off()


############################
### Graphics Theory Part ###
############################

### Composition Method ###
x <- seq(0,4, length=1000) 
hx <- dtri(x, min=1, max=3, mode=2)

png(filename="Comp.png" ,width = 800, height= 800)
plot(x, hx, type="l", lty=1, lwd=3 ,xlab = "x",
     ylab = "y",cex.axis = 2, cex.lab=2) 
abline(v=2, col = "seagreen", lwd = 5, lty = 2) #split pdf to describe it as sum
dev.off()

### Alias Method ###
al <- c(8/28, 3/28, 6/28, 11/28)
al2 <- c(7/28,3/28, 6/28, 7/28)
al3 <- c(0,7/28 ,7/28,0)
# visualization of theoretical background
png(filename = "Alias.png", width = 800, height= 800)
par(mfrow=c(1,2))
barplot(al, beside = T, col= c("seagreen", "blue", "red", "orange"), 
        ylim=c(0,0.4), cex.axis = 2.5, xlab="K", ylab = "q(K)", cex.lab=2)
abline(h=0.25, col= "darkgrey", lty= 4, lwd=4)
barplot(al3, col= c("blue", "orange", "seagreen", "red"), 
        ylim = c(0,0.4), cex.axis = 2.5, xlab="K",ylab = "q(K)", cex.lab=2)
barplot(al2, beside = T, col= c("seagreen", "blue", "red", "orange"),
        add=T, ylim = c(0,0.4), cex.axis = 2.5)
abline(h=0.25, col= "darkgrey", lty= 4, lwd=4)
dev.off()

### Ratio-of-Uniform ###
u <- runif(10000) 
v1 <- 2*u*(-log(u))^(1/2)   #v and v1 by solving u â‰¤ sqrt(f(v/u)) for v
v <- (-2)*u*(-log(u))^(1/2)
# boundaries of Cf and enclosing rectangle for N(0,1)
png(filename= "ROU.png",width = 800, height= 800) 
plot(u,v, ylim = c(-1,1), col= "seagreen", cex.axis= 2.5, cex.lab=2)
points(u,v1, col= "seagreen")
abline(h=0)
abline(h=(max(v1)), col= "blue", lwd=4)
abline(h=(min(v)), col= "blue", lwd=4)
abline(v=1, col="blue", lwd=4)
abline(v=0, col="blue", lwd=4)
dev.off()



