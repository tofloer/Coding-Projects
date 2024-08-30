##########################
## MI for PS Methods    ##
## Date: March 4, 2021  ##
##########################

##################################
## Packages and stored objects: ##
##################################

require(mi)
require(dplyr)
require(MASS)
require(boot)
require(ggplot2)
require(ggpubr)


##############################
## Data generating process: ##
##############################

## Create experimental design matrix:

rho <- c(0.3, 0.6)
RR <- c(1, 2)
gamma_y <- c(0, -0.4)
design.matrix <- expand.grid(rho, RR, gamma_y)
colnames(design.matrix) <- c("rho", "RR", "gamma")
rownames(design.matrix) <- paste0("case-", 1:8)

## Treatment assignment process as function of covariates:

logit.treatmentprob <- function(x1, x2, x3){
  logit = -1.15 + 0.7*x1 + 0.6*x2 + 0.6*x3
  return(logit)
}

## Logit of (binary) outcome probability from covariates and treatment indicator:

logit.outcome <- function(x1, x2, x3, z, RR, rho){
  
  if(rho == 0.3 & RR == 1){theta_c = 0}
  if(rho == 0.3 & RR == 2){theta_c = 1.221}
  if(rho == 0.6 & RR == 1){theta_c = 0}
  if(rho == 0.6 & RR == 2){theta_c = 1.289}
  
  logit <- -1.5 + 0.5*x1 + 0.5*x2 + 0.3*x3 + theta_c*z
  return(logit)
}

## Logit of missingness probability for covariates X1 and X3:

logit.missing.covariate <- function(gamma_0, Z, X2, gamma_y, Y){
  logit.missingness <- gamma_0 + Z + X2 + gamma_y*Y
  return(logit.missingness)
}

## Create data matrix as function of simulation parameters:

gen.dat.list <- function(design.matrix = design.matrix, 
                         N.datasets = 5,
                         n.obs = 2000,
                         condition.index){
  
  ## Get parameters from experimental design matrix:
  
  pars <- design.matrix[condition.index,]
  
  ## Mean vector for X1, X2, X3:
  
  mu <- rep(0, 3)
  
  ## Setup covariance matrix for X1, X2, X3:
  
  Sigma <- matrix(pars$rho, nrow = 3, ncol = 3)
  diag(Sigma) = 1
  colnames(Sigma) <- paste0("X", 1:3)
  rownames(Sigma) <- paste0("X", 1:3)
  
  ## Prepare dat.list to collect all the resulting datasets for
  ## the experimental condition:
  
  dat.list = list(); length(dat.list) = N.datasets
  names(dat.list) <- paste0("dataset-", 1:N.datasets)
  
  ## Start loop to create a dataset with each iteration:
  
  for(j in 1:N.datasets){
    
    if(j %% 50 == 0){print(paste0("dataset ", j, " out of ", N.datasets))}
    
    ## Draw covariates from multivariate normal:
    
    covariates.raw <- mvrnorm(n = n.obs, 
                              mu = mu, 
                              Sigma = Sigma) %>%
      as.data.frame
    
    ## Dichotomize X3:
    
    X3.pos.loc <- which(covariates.raw$X3 >= 0)
    covariates <- covariates.raw
    covariates$X3 <- 0
    covariates$X3[X3.pos.loc] <- 1
    
    ## Generate logit of treatment probability:
    
    logit.treatprobs <- apply(covariates, 1, 
                              function(x){logit.treatmentprob(x[1], x[2], x[3])})
    
    ## Generate treatment probability:
    
    treatprobs <- inv.logit(logit.treatprobs)
    
    ## Generate treatment assignment:
    
    treat <- vapply(treatprobs, 
                    function(x){rbinom(1, 1, x)},
                    FUN.VALUE = numeric(1))
    
    ## Join covariate and treatment data:
    
    dat <- data.frame(covariates, 
                      "Z" = treat)
    
    ## Generate outcome probability:
    
    logit.Y <- apply(dat, 1, FUN = function(x){
      logit.outcome(x1 = x[1], 
                    x2 = x[2], 
                    x3 = x[3], 
                    z = x[4],
                    rho = pars$rho,
                    RR = pars$RR
                    )}
      )
    
    prob.Y <- inv.logit(logit.Y)
    
    dat$outcome <- vapply(prob.Y, FUN = function(x){
      rbinom(1, 1, x)
    }, FUN.VALUE = numeric(1))
    
    ## Set value of gamma_0 that is used for missingness indicators:
    
    if(pars$gamma == 0){gamma_0 <- -1.5} else if(
      pars$gamma == -0.4){gamma_0 <- -1.3} else{
        gamma_0 <- NA
        stop("gamma_y unbekannt")
      }
    
    ## Generate logit of missingness probabilities for X1 and X3:
    
    logit.miss.covariate <- apply(dat, 1, function(x){
      logit.missing.covariate(gamma_0 = gamma_0, 
                              Z = x[4], 
                              X2 = x[2],
                              gamma_y = pars$gamma, 
                              Y = x[5])
      })
    
    ## Generate probability of missingness for X1 and X3:
    
    prob.miss.covariate <- inv.logit(logit.miss.covariate)
    
    ## Generate missingness indicators for X1 and X3:
    
    dat$miss.x1 <- vapply(prob.miss.covariate, 
                          FUN = function(x){rbinom(1, 1, x)}, 
                          FUN.VALUE = numeric(1))
    
    dat$miss.x3 <- vapply(prob.miss.covariate, 
                          FUN = function(x){rbinom(1, 1, x)}, 
                          FUN.VALUE = numeric(1))
    
    ## Duplicate Columns for X1 and X3 to keep their values 
    ## after missingness:
    
    dat$X1.full <- dat$X1
    dat$X3.full <- dat$X3
    
    ## Delete data for X1 and X3 where the missingness indicator is 1:
    
    dat$X1[dat$miss.x1 == 1] <- NA
    dat$X3[dat$miss.x3 == 1] <- NA
    
    ## Save current dataset in list of all datasets in this 
    ## experimental condition:
    
    dat.list[[j]] <- dat
    
  }
  return(dat.list)
}

## Create a list of datasets for each experimental condition:
set.seed(123)

all.dat.lists <- list(); length(all.dat.lists) <- nrow(design.matrix)
names(all.dat.lists) <- paste0("experimental-condition-",
                               1:nrow(design.matrix))

for(i in 1:nrow(design.matrix)){

   print(paste0("starting condition ", i))
 
   all.dat.lists[[i]] <- gen.dat.list(design.matrix = design.matrix,
                                      N.datasets = 1000,
                                      n.obs = 1000,
                                      condition.index = i)
}



#######################################################
## Variance and Standardized-Differences Functions:  ##
#######################################################

## Variance following Williamson et al. equation (4) and result matrices
Var.Williamson <- function(dat,ps.model,vcov,mu1,mu0) {
  
  n <- nrow(dat)
  
  w0 <- mean((1-dat$Z)/(1-dat$PS))
  w1 <- mean(dat$Z/dat$PS)
  
  K0 <- 1/mu0
  K1 <- 1/mu1
  
  V_un <- (K1/w1)**2*mean((dat$outcome-mu1)**2*dat$Z/dat$PS**2) +
          (K0/w0)**2*mean((dat$outcome-mu0)**2*(1-dat$Z)/(1-dat$PS)**2)
  
  
  mat <- as.matrix(t(cbind(rep(1,n),ps.model$model[2:4])))
  vpart1<-NULL
  vpart2<-NULL
  
  for (i in seq(1:n)){
    x <- mat[,i]
    a <- ((dat$outcome[i]-mu1)*dat$Z[i]*(1-dat$PS[i]))/dat$PS[i]
    b <- ((dat$outcome[i]-mu0)*(1-dat$Z[i])*dat$PS[i])/(1-dat$PS[i])
    ax <- a*x
    bx <- b*x
    vpart1 <- cbind(vpart1,ax)
    vpart2 <- cbind(vpart2,bx)
  }
  part1 <- rowMeans(vpart1)
  part2 <- rowMeans(vpart2)
  v <- (K1/w1)*part1 + (K0/w0)*part2
  
  nVar <- V_un - t(v)%*% (vcov*n) %*%v
  
  var_corr <- nVar/n
 return(var_corr) 
}

Var.Full <- Var.CC <- Var.MIte <- Var.MIps <- Var.MIpar <- matrix(nrow = 1000, ncol = 8)


## SDiff and result matrices
SDiff <- function(vec1,vec0) {
  SD <- (100*abs(mean(vec1) - mean(vec0)))/sqrt((var(vec1) + var(vec0))/2)
  return(SD) 
}

SD.Crude <- SD.Full <- SD.CC <- SD.MIte <- SD.MIps <- SD.MIpar <- matrix(nrow = 1000, ncol = 8) 

####################################
## Method 0: Before Data Deletion ##
####################################

if(!"RR.list" %in% ls(envir = globalenv())){
  RR.list <- list()
}

if(!"Full" %in% names(RR.list)){
  
  Full <- list(); length(Full) = length(all.dat.lists)
  names(Full) = names(all.dat.lists)
  
  for(i in 1:length(all.dat.lists)){
    dat.list.iter <- all.dat.lists[[i]]
    
    print(paste0("Experimental Condition ", i))
    
    for(j in 1:length(dat.list.iter)){
      dat.iter <- dat.list.iter[[j]]

      ps.model <- glm(Z ~ X1.full + X2 + X3.full,
                      data = dat.iter,
                      family = binomial(link = "logit"),
                      maxit = 100)
      
      dat.iter$PS <- inv.logit(predict(ps.model))
      
      # Calculate IPTW Treatment Effect:
      
      dat.treated <- dat.iter[which(dat.iter$Z == 1),]
      num <- sum(dat.treated$outcome/dat.treated$PS)
      denom <- sum(1/dat.treated$PS)
      mu1.iptw <- num/denom
      
      dat.untreated <- dat.iter[which(dat.iter$Z == 0),]
      num <- sum(dat.untreated$outcome/(1 - dat.untreated$PS))
      denom <- sum(1/(1 - dat.untreated$PS))
      mu0.iptw <- num/denom
      
      Full[[i]] <- c(Full[[i]], mu1.iptw/mu0.iptw)
      
      # Calculate Standardized Differences:
      
      SD.Crude[j,i] <- SDiff(dat.treated$X1.full, dat.untreated$X1.full)
      
      SD.Full[j,i] <- SDiff(dat.treated$X1.full/dat.treated$PS, dat.untreated$X1.full/(1-dat.untreated$PS))
      
      # Calculate Variance:
      
      Var.Full[j,i] <- Var.Williamson(dat.iter, ps.model, vcov(ps.model), mu1.iptw, mu0.iptw)
      
      
      if(j %% 50 == 0){print(paste0("dataset ", j, " out of ", 
                                    length(dat.list.iter)))}
    }
  }
  RR.list$Full <- Full
}

RR.Full.mean <- lapply(RR.list$Full, mean) %>%
  unlist %>% 
  round(2)


##############################
## Method 1: Complete Cases ##
##############################

if(!"CC" %in% names(RR.list)){
  
  CC <- list(); length(CC) = length(all.dat.lists)
  names(CC) = names(all.dat.lists)
 
  for(i in 1:length(all.dat.lists)){
    dat.list.iter <- all.dat.lists[[i]]
    
    print(paste0("Experimental Condition ", i))
    
    for(j in 1:length(dat.list.iter)){
      dat.iter <- dat.list.iter[[j]]
      dat.cc <- na.omit(dat.iter)
      
      ps.model <- glm(Z ~ X1 + X2 + X3, 
                      data = dat.cc,
                      family = binomial(link = "logit"),
                      maxit = 100)
      
      dat.cc$PS <- inv.logit(predict(ps.model))
      
      # Calculate IPTW Treatment Effect:
      
      dat.treated <- dat.cc[which(dat.cc$Z == 1),]
      num <- sum(dat.treated$outcome/dat.treated$PS)
      denom <- sum(1/dat.treated$PS)
      mu1.iptw <- num/denom
      
      dat.untreated <- dat.cc[which(dat.cc$Z == 0),]
      num <- sum(dat.untreated$outcome/(1 - dat.untreated$PS))
      denom <- sum(1/(1 - dat.untreated$PS))
      mu0.iptw <- num/denom
      
      CC[[i]] <- c(CC[[i]], mu1.iptw/mu0.iptw)
      
      # Calculate Standardized Differences:

      SD.CC[j,i] <- SDiff(dat.treated$X1.full/dat.treated$PS, dat.untreated$X1.full/(1-dat.untreated$PS))
      
      # Calculate Variance:
      
      Var.CC[j,i] <- Var.Williamson(dat.cc, ps.model, vcov(ps.model), mu1.iptw, mu0.iptw)
      
      if(j %% 50 == 0){print(paste0("dataset ", j, " out of ", 
                                    length(dat.list.iter)))}
    }
  }
  RR.list$CC <- CC
}

RR.CC.mean <- lapply(RR.list$CC, mean) %>%
  unlist %>% 
  round(2)


##############################
## Methods: MIte,MIps,MIpar ##
##############################

## Number of Imputations:
set.seed(123)
M = 10

strt<-Sys.time()
if(!"MIte" %in% names(RR.list)){
  
  MIte <- list(); length(MIte) = length(all.dat.lists)
  names(MIte) = names(all.dat.lists)
  
  MIps <- list(); length(MIps) = length(all.dat.lists)
  names(MIps) = names(all.dat.lists)
  
  MIpar <- list(); length(MIpar) = length(all.dat.lists)
  names(MIpar) = names(all.dat.lists)
  
  for(i in 1:length(all.dat.lists)){
    dat.list.iter <- all.dat.lists[[i]]
    
    print(paste0("Experimental Condition ", i))
    
    for(j in 1:length(dat.list.iter)){
      dat.iter <- dat.list.iter[[j]][c("X1","X2","X3","Z","outcome")]
      
      mi.imp <- mi(subset(dat.iter, select = -c(4)) , n.chains = M, n.iter = 15, parallel = TRUE, seed = 123)
      dat.imp <- complete(mi.imp) # checked covergence for n.iter with Rhats(mi.imp, statistic = "moments")
      
      mu0.iptw.MIte <- mu1.iptw.MIte <- SD.MI <- Var.MI <- SD.MI.full <- vector(length = M)
      ps.MIps <- matrix(nrow = nrow(dat.iter), ncol = M)
      par.MIpar <- matrix(nrow = M, ncol = 4)
      X1.MIpar <- X2.MIpar <- X3.MIpar <- matrix(nrow = nrow(dat.iter), ncol = M) #rep(0L, nrow(dat.iter))
      
      for(h in 1:M){
        dat.MI.raw <- cbind(subset(dat.iter, select = c(4)), dat.imp[[h]])
        
        loc.fac <- which(sapply(dat.MI.raw, is.factor))
        loc.log <- which(sapply(dat.MI.raw, is.logical))
        
        loc.recode <- c(loc.fac, loc.log)
        loc.ok <- setdiff(c(1:ncol(dat.MI.raw)), loc.recode)
        
        dat.MI.ok <- dat.MI.raw[, loc.ok]
        dat.MI.fac <- dat.MI.raw[, loc.fac]
        dat.MI.log <- dat.MI.raw[, loc.log]
        
        dat.MI.fac.rec <- apply(dat.MI.fac, 
                                2, 
                                function(x) as.numeric(as.vector(x))) %>%
          data.frame
        
        dat.MI.log.rec <- apply(dat.MI.log, 
                                2, 
                                function(x) as.numeric(x)) %>%
          data.frame
        
        dat.MI <- data.frame(matrix(NA, 
                                    ncol = ncol(dat.MI.raw), 
                                    nrow = nrow(dat.MI.raw)))
        
        dat.MI[,loc.ok] <- dat.MI.ok
        dat.MI[,loc.fac] <- dat.MI.fac.rec
        dat.MI[,loc.log] <- dat.MI.log.rec
        colnames(dat.MI) <- colnames(dat.MI.raw)

        
        ps.model <- glm(Z ~ X1 + X2 + X3, 
                        data = dat.MI,
                        family = binomial(link = "logit"),
                        maxit = 100)
        
        dat.MI$PS <- inv.logit(predict(ps.model))
        
        # Calculate MIte IPTW Treatment Effect:
        dat.treated <- dat.MI[which(dat.MI$Z == 1),]
        num <- sum(dat.treated$outcome/dat.treated$PS)
        denom <- sum(1/dat.treated$PS)
        mu1.iptw.MIte[h] <- num/denom
        
        dat.untreated <- dat.MI[which(dat.MI$Z == 0),]
        num <- sum(dat.untreated$outcome/(1 - dat.untreated$PS))
        denom <- sum(1/(1 - dat.untreated$PS))
        mu0.iptw.MIte[h] <- num/denom
        
        # Calculate Standardized Differences MIte:
        SD.MI[h] <- SDiff(dat.treated$X1/dat.treated$PS, dat.untreated$X1/(1-dat.untreated$PS))
        
        # Calculate Variance MIte:
        Var.MI[h] <- Var.Williamson(dat.MI, ps.model, vcov(ps.model), mu1.iptw.MIte[h], mu0.iptw.MIte[h])
        
        # Collect necessary data for MIps and MIpar from M imputed sets
        ps.MIps[,h] <- inv.logit(predict(ps.model))
        par.MIpar[h,] <-  ps.model$coefficients
        X1.MIpar[,h] <- dat.imp[[h]]$X1
        X2.MIpar[,h] <- dat.imp[[h]]$X2
        X3.MIpar[,h] <- as.numeric(dat.imp[[h]]$X3)-1
      }
      
        # Combine MIte calculations
        SD.MIte[j,i] <- mean(SD.MI)
        Var.MIte[j,i] <- mean(Var.MI)
        
        
        # Bind MIpar data
        dat.MIpar <- cbind(1,rowMeans(X1.MIpar, dims = 1), rowMeans(X2.MIpar, dims = 1), rowMeans(X3.MIpar, dims = 1))

        
        # Calculate MIps IPTW Treatment Effect:
        PS <- rowMeans(ps.MIps, dims = 1)
        dat.MIps <- cbind(subset(dat.iter, select = c(4,5)), PS, dat.MIpar)
        
        dat.treated <- dat.MIps[which(dat.MIps$Z == 1),]
        num <- sum(dat.treated$outcome/dat.treated$PS)
        denom <- sum(1/dat.treated$PS)
        mu1.iptw.MIps <- num/denom
        
        dat.untreated <- dat.MIps[which(dat.MIps$Z == 0),]
        num <- sum(dat.untreated$outcome/(1 - dat.untreated$PS))
        denom <- sum(1/(1 - dat.untreated$PS))
        mu0.iptw.MIps <- num/denom
        
        dat.vcov <- as.data.frame(cbind(dat.MI$Z,dat.MIpar[,2:4]))
        ps.model <- glm(V1 ~ V2 + V3 + V4, 
                        data = dat.vcov,
                        family = binomial(link = "logit"),
                        maxit = 100)
        
        # Calculate Standardized Differences MIps:
        SD.MIps[j,i] <- SDiff(unlist(dat.treated[5], use.names=FALSE)/dat.treated$PS, unlist(dat.untreated[5], use.names=FALSE)/(1-dat.untreated$PS))
        
        # Calculate Variance MIps:
        Var.MIps[j,i] <- Var.Williamson(dat.MIps, ps.model, vcov(ps.model), mu1.iptw.MIps, mu0.iptw.MIps)
        
        
        # Calculate MIpar IPTW Treatment Effect:
        par.MI.par <- colMeans(par.MIpar)
        
        PS <- exp(dat.MIpar%*%par.MI.par)/(1+exp(dat.MIpar%*%par.MI.par))

        dat.MIpar.ps <- cbind(subset(dat.iter, select = c(4,5)), PS, dat.MIpar)
        
        
        dat.treated <- dat.MIpar.ps[which(dat.MIpar.ps$Z == 1),]
        num <- sum(dat.treated$outcome/dat.treated$PS)
        denom <- sum(1/dat.treated$PS)
        mu1.iptw.MIpar <- num/denom
        
        dat.untreated <- dat.MIpar.ps[which(dat.MIpar.ps$Z == 0),]
        num <- sum(dat.untreated$outcome/(1 - dat.untreated$PS))
        denom <- sum(1/(1 - dat.untreated$PS))
        mu0.iptw.MIpar <- num/denom
        
        # Calculate Standardized Differences MIpar:
        SD.MIpar[j,i] <- SDiff(unlist(dat.treated[5], use.names=FALSE)/dat.treated$PS, unlist(dat.untreated[5], use.names=FALSE)/(1-dat.untreated$PS))
        
        # Calculate Variance MIpar:
        Var.MIpar[j,i] <- Var.Williamson(dat.MIpar.ps, ps.model, vcov(ps.model)-cov(par.MIpar), mu1.iptw.MIps, mu0.iptw.MIpar)   # Cov(alpha_bar) = (1+1/M)B
        
      
            
      
      MIte[[i]] <- c(MIte[[i]], mean(mu1.iptw.MIte/mu0.iptw.MIte))
      MIps[[i]] <- c(MIps[[i]], mu1.iptw.MIps/mu0.iptw.MIps)
      MIpar[[i]] <- c(MIpar[[i]], mu1.iptw.MIpar/mu0.iptw.MIpar)
      
      
      if(j %% 25 == 0){print(paste0("dataset ", j, " out of ", 
                                    length(dat.list.iter)))}
    }
  }
  RR.list$MIte <- MIte
  RR.list$MIps <- MIps
  RR.list$MIpar <- MIpar
}
print(Sys.time()-strt)

RR.MIte.mean <- lapply(RR.list$MIte, mean) %>%
  unlist %>% 
  round(2)

RR.MIps.mean <- lapply(RR.list$MIps, mean) %>%
  unlist %>% 
  round(2)

RR.MIpar.mean <- lapply(RR.list$MIpar, mean) %>%
  unlist %>% 
  round(2)



#################################################
## Analysis of Bias, Coverage Rates and SDiff: ##
#################################################

## Design matrix
rho <- c(0.3, 0.6)
RR <- c(1, 2)
gamma_y <- c(0, -0.4)
design.matrix <- expand.grid(rho, RR, gamma_y)
colnames(design.matrix) <- c("rho", "RR", "gamma")
rownames(design.matrix) <- paste0("case-", 1:8)


## Estimated variance coverage rates with log transformation for CI
tval <- c(1,1,2,2,1,1,2,2)

bias <- cover <- wil.var <- matrix(nrow=5,ncol=8)
wil.var <- rbind(colMeans(Var.BD),colMeans(Var.CC),colMeans(Var.MIte),colMeans(Var.MIps),colMeans(Var.MIpar))

for (k in 1:8) {
  for (j in 1:5) {
    Q <- unlist(RR.list[[j]][k], use.names=FALSE)
    
    bias[j,k] <- mean(Q) - tval[k] # bias
    
    T <- wil.var[j,k] # variance
    CIlow <- log(Q)-qnorm(0.975)*sqrt(T) #confidence intervals
    CIup <- log(Q)+qnorm(0.975)*sqrt(T)
    c <- rep(0L, 1000)
    for (i in 1:1000){
      c[i] <- ifelse(exp(CIlow[i]) <= tval[k] && exp(CIup[i]) >= tval[k],1,0) #coverage
    }
    cover[j,k] <- mean(c)
  }
}
rownames(bias) <- rownames(cover) <- rownames(wil.var) <- c( "Full","CC","MIte", "MIps", "MIpar")
colnames(bias) <- colnames(cover) <- colnames(wil.var) <- paste0("case-", 1:8)


## Empirical coverage rates with log transformation for CI
emp.cover <- emp.var <- matrix(nrow=5,ncol=8)
for (k in 1:8) {
  for (j in 1:5) {
    Q <- unlist(RR.list[[j]][k], use.names=FALSE)
    
    T <- emp.var[j,k] <- var(log(Q)) # variance
    CIlow <- log(Q)-qnorm(0.975)*sqrt(T) #confidence intervals
    CIup <- log(Q)+qnorm(0.975)*sqrt(T)
    c <- rep(0L, 1000)
    for (i in 1:1000){
      c[i] <- ifelse(exp(CIlow[i]) <= tval[k] && exp(CIup[i]) >= tval[k],1,0) #coverage
      }
    emp.cover[j,k] <- mean(c)
  }
}
rownames(emp.cover) <- rownames(emp.var) <- c( "Full","CC","MIte", "MIps", "MIpar")
colnames(emp.cover) <- colnames(emp.var) <- paste0("case-", 1:8)


## Standardized Differences
SD <- matrix(nrow=6,ncol=8)

SD <- rbind(colMeans(SD.crude),colMeans(SD.full),colMeans(SD.CC),colMeans(SD.MIte),colMeans(SD.MIps),colMeans(SD.MIpar))
rownames(SD) <- c( "Crude","Full","CC","MIte", "MIps", "MIpar")
colnames(SD) <- paste0("case-", 1:8)


## Results
print(design.matrix)
print(bias)
print(cover)
print(emp.cover)
print(wil.var)
print(emp.var)
print(SD)



############
## Plots: ##
############

## Missingness associated with Y
Bias<-bias

D1<-data.frame(Methods=c("Full", "CC", "MIte", "MIps", "MIpar"), Bias=c(as.vector(Bias[,6])))

Plot1<-ggplot(data=D1, aes(x=Methods, y=Bias, fill=Methods)) + ylim(0,0.35) + geom_bar(stat="identity", show.legend=FALSE)+theme(plot.margin = margin(1.5, 0.5, 0.5, 1, "cm"),
                                                                                                                 axis.line = element_line(colour = "black"),
                                                                                                                 panel.grid.major = element_blank(),
                                                                                                                 panel.grid.minor = element_blank(),
                                                                                                                 panel.border = element_blank(),
                                                                                                                 axis.title.x=element_blank(),
                                                                                                                 panel.background = element_blank(),
                                                                                                                 plot.title= element_text(face="bold", size=8, hjust=0.5)) + ggtitle("Missingness associated with Y")
Plot1


## Missingness independent of Y
D2<-data.frame(Methods=c("Full", "CC", "MIte", "MIps", "MIpar"), Bias=c(as.vector(Bias[,2])) )

Plot2<-ggplot(data=D2, aes(x=Methods, y=Bias, fill=Methods)) + ylim(0,0.35) + geom_bar(stat="identity", show.legend=FALSE)+theme(plot.margin = margin(1.5, 0.5, 0.5, 1, "cm"),
                                                                                                                axis.line = element_line(colour = "black"),
                                                                                                                panel.grid.major = element_blank(),
                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                panel.border = element_blank(),
                                                                                                                axis.title.x=element_blank(),
                                                                                                                panel.background = element_blank(),
                                                                                                                plot.title= element_text(face="bold", size=8, hjust=0.5)) + ggtitle("Missingness independent of Y")
Plot2


## Missingness associated with Y
D3<-data.frame(Methods=c("Full", "CC", "MIte", "MIps", "MIpar"), Bias= c(as.vector(Bias[,8])))

Plot3<-ggplot(data=D3, aes(x=Methods, y=Bias, fill=Methods)) + ylim(0,0.35) + geom_bar(stat="identity", show.legend=FALSE)+theme(plot.margin = margin(1.5, 0.5, 0.5, 1, "cm"),
                                                                                                                axis.line = element_line(colour = "black"),
                                                                                                                panel.grid.major = element_blank(),
                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                panel.border = element_blank(),
                                                                                                                axis.title.x=element_blank(),
                                                                                                                panel.background = element_blank(),
                                                                                                                plot.title= element_text(face="bold", size=8, hjust=0.5)) + ggtitle("Missingness associated with Y")
Plot3


## Missingness independent of Y
D4<-data.frame(Methods=c("Full", "CC", "MIte", "MIps", "MIpar"), Bias= c(as.vector(Bias[,4])))

Plot4<-ggplot(data=D4, aes(x=Methods, y=Bias,fill=Methods)) + ylim(0,0.35) + geom_bar(stat="identity", show.legend=FALSE) +theme(plot.margin = margin(1.5, 0.5, 0.5, 1, "cm"),
                                                                                                                axis.line = element_line(colour = "black"),
                                                                                                                panel.grid.major = element_blank(),
                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                panel.border = element_blank(),
                                                                                                                axis.title.x=element_blank(),
                                                                                                                panel.background = element_blank(),
                                                                                                                plot.title= element_text(face="bold", size=8, hjust=0.5)) + ggtitle("Missingness independent of Y")
Plot4


## Arrangement of the plots
Plotfull <- ggarrange(Plot1,Plot2, Plot3, Plot4,ncol = 2, nrow = 2, labels=c("RR=1", "RR=2"), label.x=c(0.85,-0.15), label.y=c(1,0.117), vjust=4, font.label = list(size = 10, color = "black", face = "bold", family = NULL))
Plotfull


Cover<-cover
S<-data.frame(Methods=c(replicate(8,"Full"), replicate(8,"CC"), replicate(8,"MIte"), replicate(8,"MIps"), replicate(8, "MIpar")), CI=c(Cover[1,],Cover[2,],Cover[3,],Cover[4,],Cover[5,]))
Figure3<-ggplot(S, aes(x=Methods, y= CI, fill=Methods)) + ylim(0.8,1) + geom_boxplot(show.legend = FALSE)+theme(plot.margin = margin(1, 1, 1, 1, "cm"),axis.line = element_line(colour = "black"),
                                                                                                panel.grid.major = element_blank(),
                                                                                                panel.grid.minor = element_blank(),
                                                                                                panel.border = element_blank(),
                                                                                                axis.text.x = element_text(face="bold", size=15, color="black"),
                                                                                                axis.title.x= element_blank(),
                                                                                                axis.title.y=element_text(face="bold", size=12),
                                                                                                panel.background = element_blank())  + ylab("Coverage Rate")+ geom_hline(yintercept=0.95, linetype="dashed", color = "black")
Figure3  


