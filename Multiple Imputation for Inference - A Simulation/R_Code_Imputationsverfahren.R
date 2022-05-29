# R-Code zur Fu:Stat Seminararbeit "Inferenz mit unvollständigen Daten - ein Vergleich von singulären und multiplen Imputationsverfahren" von Tobias Flörchinger (608227, Master Statistik)

##################################
########## Pakete Laden ##########
##################################

### Benutzte Pakete
packages =  c("mice"       #Imputationen und Rubins Rules zu beta2, beta3
             ,"miceadds"   #Combining Rules für Cor/R^2 (für Ergebnisse der Arbeit nicht unmittelbar relevant!)
             ,"rpart"      #CART
             ,"VIM"        #Plots
             ,"plot3D"     #Plots
             ,"rgl"        #Plots
             ,"rpart.plot" #Plots
              )

### Installation der Pakete, falls nicht vorhanden
inst <- packages %in% installed.packages()    
if(length(packages[!inst]) > 0) install.packages(packages[!inst])

### Pakete in die Umgebung einbinden
loading <- lapply(packages, require, character.only=TRUE)



######################################
##########   Vorbereitung   ##########
######################################

##########   Wahre Werte (mit WLLN)   ##################################################
set.seed(123)

K <- 100000
n <- 1000
true_vals <- matrix(nrow=K,ncol=3)

### Iterative Erstellung der Zufallsstichproben mit den Schätzungen zu beta2, beta3 und mu4
strt<-Sys.time()
for(i in 1:K) {
if(i %% 5000 == 0) cat("Iteration", i, "of", K,"\n")
  Y1 <- rnorm(n,3,1)
  Y2 <- runif(n,min=2,max=4)
  Y3 <- 0.75*Y1 + 0.1*Y2 + rnorm(n,0,1)
  Y4 <- 3*Y1 + 0.75*Y2^2 + 0.5*Y3 + (rchisq(n,3)-3)
  data <- as.data.frame(cbind(Y1,Y2,Y3,Y4)) 
  fit <- lm(Y4~Y1+Y2+Y3, data=data)
  true_vals[i,1] <- mean(Y4)
  true_vals[i,2:3] <- coef(fit)[3:4]
}
print(Sys.time()-strt)

### Gemittelte Schätzungen als Annäherung der wahren Werte
true_vals <- colMeans(true_vals)


##########   Funktionen für die Simulation   ##################################################

### Konfidenzintervalle und deren Breite zu beta2, beta3 und mu4
CI_Vals <- function(dat,model) {
  ci_low_b2 <- confint(model, 'Y2', level=0.95)[1]
  ci_up_b2 <- confint(model, 'Y2', level=0.95)[2]
  ci_width_b2 <- ci_up_b2 - ci_low_b2
  ci_low_b3 <- confint(model, 'Y3', level=0.95)[1]
  ci_up_b3 <- confint(model, 'Y3', level=0.95)[2]
  ci_width_b3 <- ci_up_b3 - ci_low_b3
  ci_low_Y4 <- t.test(dat$Y4)$conf.int[1]
  ci_up_Y4 <- t.test(dat$Y4)$conf.int[2]
  ci_width_Y4 <- ci_up_Y4 - ci_low_Y4
  return(cbind(ci_low_Y4,ci_up_Y4,ci_width_Y4,ci_low_b2,ci_up_b2,ci_width_b2,ci_low_b3,ci_up_b3,ci_width_b3))
}

### Deckungsrate
Coverage <- function(value, CI_low, CI_up) {
  ifelse(CI_low <= value && CI_up >= value,1,0)
}

### Rubins Rules für mu4
MI_comb <- function(dat,M){
  Qhat <- aggregate(Y4 ~.imp, data=dat, mean)$Y4  # Schätzer der einzelnen ergänten Datensätze h
  Uhat <- aggregate(Y4 ~.imp, data=dat, var)$Y4/n # Varianz der Schätzer für h
  Qbar <- sum(Qhat)/M # Kombinierter Schätzer
  Ubar <- sum(Uhat)/M # Within-Varianz Ubar
  B <- (1/(M-1))*sum((Qhat-Qbar)^2) # Between-Varianz B
  T <- Ubar+B+B/M # total variance
  df <- (M-1)*(1+(M/(M+1))*Ubar/B)^2 # Freiheitsgrade
  # Konfidenzintervall und dessen Breite
  CIlow <- Qbar-qt(0.975,df)*sqrt(T) 
  CIupper <- Qbar+qt(0.975,df)*sqrt(T)
  CIwidth <- CIupper - CIlow
  return(cbind(Qbar,CIlow,CIupper,CIwidth,T))
}

### Liste aller Vergleichskriterien zu mu4, beta2, beta3 für singuläre Imputationen (gegenüber der Seminararbeit um Cor(X3,X4) und R^2 erweitert)
criteria_s <- function(dat) {
  l <- c(c((mean(dat$Y4)-true_vals[1])/true_vals[1], Coverage(true_vals[1],ci_vals[1],ci_vals[2]), ci_vals[3], var(dat$Y4)/n, cor(dat$Y3,dat$Y4)), # zu mu4
            c((coef(fit)[3]-true_vals[2])/true_vals[2], Coverage(true_vals[2],ci_vals[4],ci_vals[5]), ci_vals[6], sqrt(diag(vcov(fit)))[3], # zu beta2
              (coef(fit)[4]-true_vals[3])/true_vals[3], Coverage(true_vals[3],ci_vals[7],ci_vals[8]), ci_vals[9], sqrt(diag(vcov(fit)))[4], # zu beta3
              summary(fit)$r.squared)) 
  return(l)
}

### Liste aller Vergleichskriterien zu mu4, beta2, beta3 für multiple Imputationen (gegenüber der Seminararbeit mit miceadds um Cor(X3,X4) und R^2 erweitert)
criteria_m <- function(MI_ana) {
  l <- c(c((MI_ana[1]-true_vals[1])/true_vals[1], Coverage(true_vals[1], MI_ana[2], MI_ana[3]),MI_ana[4], MI_ana[5], cor[1,3]), # zu mu4
            c((sum_reg[3,2]-true_vals[2])/true_vals[2], Coverage(true_vals[2],sum_reg[3,7],sum_reg[3,8]), sum_reg[3,8]-sum_reg[3,7], sum_reg[3,3], # zu beta2
              (sum_reg[4,2]-true_vals[3])/true_vals[3], Coverage(true_vals[3],sum_reg[4,7],sum_reg[4,8]), sum_reg[4,8]-sum_reg[4,7], sum_reg[4,3],# zu beta3
              R2))
  return(l)
}



##########################################
##########   Simulation Study   ##########
##########################################

set.seed(123)
Z <- 1000
runs_T <- 5 # MICE Iterationen
m <- 5 # Anzahl der Imputationen

### Tabellen der Schätzungen, die mit jeder Iteration z befüllte werden
BD <- CC <- SI_mean <- SI_hot <- SI_reg <- SI_norm <- SI_pmm <- SI_cart <- MI_norm <- MI_pmm <- MI_cart <- matrix(nrow=Z,ncol=5)
BD_reg <- CC_reg <- SI_mean_reg <- SI_hot_reg <- SI_reg_reg <- SI_norm_reg <- SI_pmm_reg <- SI_cart_reg <- MI_norm_reg <- MI_pmm_reg <- MI_cart_reg <- matrix(nrow=Z,ncol=9)
PD <- matrix(nrow=Z,ncol=2)


strt<-Sys.time()
for(z in 1:Z) { if(z %% 50 == 0) cat("Iteration", z, "of", Z,"...\n")
  
  ### Datengenerierender Prozess
  Y1 <- rnorm(n,3,1)
  Y2 <- runif(n,min=2,max=4)
  Y3 <- 0.75*Y1 + 0.1*Y2 + rnorm(n,0,1)
  Y4 <- 3*Y1 + 0.75*Y2^2 + 0.5*Y3 + (rchisq(n,3)-3)
  data <- as.data.frame(cbind(Y1,Y2,Y3,Y4)) #Datensatz in z
  
  ### Datenausfälle 
  data_miss <- data
  logit_mar <- Y1 + Y2 - 7 #für fehlende Werte in Y3
  prob <- exp(logit_mar)/(1+exp(logit_mar))
  miss_mar <- rbinom(n, 1, prob)
  data_miss$Y3[as.logical(miss_mar)] <- NA
  logit_mar <- 0.725*Y2 - 1*Y1 #für fehlende Werte in Y4
  prob <- exp(logit_mar)/(1+exp(logit_mar))
  miss_mar <- rbinom(n, 1, prob)
  data_miss$Y4[as.logical(miss_mar)] <- NA
  PD[z,1] <- sum(is.na(data_miss$Y3)/length(data_miss$Y3)) #Anteil der fehlenden Werte
  PD[z,2] <- sum(is.na(data_miss$Y4)/length(data_miss$Y4))
  
  ### Before Deletion (BD) 
  fit <- lm(Y4~Y1+Y2+Y3, data=data) #Schätzungen
  ci_vals <- CI_Vals(data,fit)
  BD[z,] <- criteria_s(data)[1:5] #Vergleichskriterien
  BD_reg[z,] <- criteria_s(data)[6:14]
  
  ### Complete Cases (CC) (In der Arbeit nicht genutzt, aber mögliche Alternative zur Imputation für valider Inferenz zu beta2 und beta3)
  data_cc <- data_miss[complete.cases(data_miss),]
  fit <- lm(Y4~Y1+Y2+Y3, data=data_cc) #Schätzungen
  ci_vals <- CI_Vals(data_cc,fit)
  CC[z,] <- criteria_s(data_cc)[1:5] #Vergleichskriterien
  CC_reg[z,] <- criteria_s(data_cc)[6:14]
  
  
  ##########   Sinuläre Imputation   ##################################################
  
  ### Mittelwertergänzung
  data_meanImp <- data_miss #Imputation
  data_meanImp$Y4[is.na(data_meanImp$Y4)] <- mean(data_meanImp$Y4, na.rm=TRUE)
  data_meanImp$Y3[is.na(data_meanImp$Y3)] <- mean(data_meanImp$Y3, na.rm=TRUE)
  fit <- lm(Y4~Y1+Y2+Y3, data=data_meanImp) #Schätzungen
  ci_vals <- CI_Vals(data_meanImp,fit)
  SI_mean[z,] <- criteria_s(data_meanImp)[1:5] #Vergleichskriterien
  SI_mean_reg[z,] <- criteria_s(data_meanImp)[6:14]
  
  ### Hot-deck Imputation
  data_hotImp <- data_miss #Imputation
  data_hotImp$Y3[is.na(data_hotImp$Y3)] <- sample(Y3[which(!is.na(data_hotImp$Y3))],length(data_hotImp$Y3[is.na(data_hotImp$Y3)]),replace=TRUE)
  data_hotImp$Y4[is.na(data_hotImp$Y4)] <- sample(Y4[which(!is.na(data_hotImp$Y4))],length(data_hotImp$Y4[is.na(data_hotImp$Y4)]),replace=TRUE)
  fit <- lm(Y4~Y1+Y2+Y3, data=data_hotImp) #Schätzungen
  ci_vals <- CI_Vals(data_hotImp,fit)
  SI_hot[z,] <- criteria_s(data_hotImp)[1:5] #Vergleichskriterien
  SI_hot_reg[z,] <- criteria_s(data_hotImp)[6:14]
  
  ### MICE Prediction Matrix initialisieren
  ini <- mice(data_miss,m=1,maxit=0)
  pred <- ini$pred
  
  ### Lineare Regressions Imputation
  imp_reg <- mice(data_miss, m=1, maxit=runs_T, method=c('','','norm.predict','norm.predict'), predictorMatrix = pred, print=F) #Imputation
  mice_reg <- complete(imp_reg, include=FALSE)
  fit <- lm(Y4~Y1+Y2+Y3, data=mice_reg) #Schätzungen
  ci_vals <- CI_Vals(mice_reg,fit)
  SI_reg[z,] <- criteria_s(mice_reg)[1:5] #Vergleichskriterien
  SI_reg_reg[z,] <- criteria_s(mice_reg)[6:14]
  
  ### Predictiv Mean Matching SI
  imp_pmm <- mice(data_miss, m=1, maxit=runs_T, method=c('','','pmm','pmm'), predictorMatrix = pred, print=F) #Imputation
  mice_pmm <- complete(imp_pmm, include=FALSE)
  fit <- lm(Y4~Y1+Y2+Y3, data=mice_pmm) #Schätzungen
  ci_vals <- CI_Vals(mice_pmm,fit)
  SI_pmm[z,] <- criteria_s(mice_pmm)[1:5] #Vergleichskriterien
  SI_pmm_reg[z,] <- criteria_s(mice_pmm)[6:14]
  
  ### Bayesian Linear Regression SI 
  imp_norm <- mice(data_miss, m=1, maxit=runs_T, method=c('','','norm','norm'), predictorMatrix = pred, print=F) #Imputation
  mice_norm <- complete(imp_norm, include=FALSE)
  fit <- lm(Y4~Y1+Y2+Y3, data=mice_norm) #Schätzungen
  ci_vals <- CI_Vals(mice_norm,fit)
  SI_norm[z,] <- criteria_s(mice_norm)[1:5] #Vergleichskriterien
  SI_norm_reg[z,] <- criteria_s(mice_norm)[6:14]
  
  ### CART Regression Tree SI 
  imp_cart <- mice(data_miss, m=1, maxit=runs_T, method=c('','','cart','cart'), predictorMatrix = pred, print=F) #Imputation
  mice_cart <- complete(imp_cart, include=FALSE)
  fit <- lm(Y4~Y1+Y2+Y3, data=mice_cart) #Schätzungen
  ci_vals <- CI_Vals(mice_cart,fit)
  SI_cart[z,] <- criteria_s(mice_cart)[1:5] #Vergleichskriterien
  SI_cart_reg[z,] <- criteria_s(mice_cart)[6:14]
  
  
  ##########   Multiple Imputation   ##################################################
  
  ### Predictiv Mean Matching MI
  imp_pmm_MI <- mice(data_miss, m=m, maxit=runs_T, method=c('','','pmm','pmm'), predictorMatrix = pred, print=F) #Imputation
  mice_pmm_MI <- complete(imp_pmm_MI, action="long", include=FALSE)
  ana_pmm <- MI_comb(mice_pmm_MI,m) #Schätzungen
  sum_reg <- summary(pool(with(imp_pmm_MI, lm(Y4~Y1+Y2+Y3))), conf.int = TRUE)
  cor <- miceadds::micombine.cor(imp_pmm_MI, variables=c(4,3))
  R2 <- pool.r.squared(pool(with(imp_pmm_MI, lm(Y4~Y1+Y2+Y3))), adjusted = FALSE)[1]
  MI_pmm[z,] <- criteria_m(ana_pmm)[1:5] #Vergleichskriterien
  MI_pmm_reg[z,] <- criteria_m(ana_pmm)[6:14]
  
  ### Bayesian Linear Regression MI 
  imp_norm_MI <- mice(data_miss, m=m, maxit=runs_T, method=c('','','norm','norm'), predictorMatrix = pred, print=F) #Imputation
  mice_norm_MI <- complete(imp_norm_MI, action="long", include=FALSE)
  ana_norm <- MI_comb(mice_norm_MI,m) #Schätzungen
  sum_reg <- summary(pool(with(imp_norm_MI, lm(Y4~Y1+Y2+Y3))), conf.int = TRUE)
  cor <- miceadds::micombine.cor(imp_norm_MI, variables=c(4,3))
  R2 <- pool.r.squared(pool(with(imp_norm_MI, lm(Y4~Y1+Y2+Y3))), adjusted = FALSE)[1]
  MI_norm[z,] <- criteria_m(ana_norm)[1:5] #Vergleichskriterien
  MI_norm_reg[z,] <- criteria_m(ana_norm)[6:14]
  
  ### CART Regression Tree MI
  imp_cart_MI <- mice(data_miss, m=m, maxit=runs_T, method=c('','','cart','cart'), predictorMatrix = pred, print=F) #Imputation
  mice_cart_MI <- complete(imp_cart_MI, action="long", include=FALSE)
  ana_cart <- MI_comb(mice_cart_MI,m) #Schätzungen
  sum_reg <- summary(pool(with(imp_cart_MI, lm(Y4~Y1+Y2+Y3))), conf.int = TRUE)
  cor <- miceadds::micombine.cor(imp_cart_MI, variables=c(4,3))
  R2 <- pool.r.squared(pool(with(imp_cart_MI, lm(Y4~Y1+Y2+Y3))), adjusted = FALSE)[1]
  MI_cart[z,] <- criteria_m(ana_cart)[1:5] #Vergleichskriterien
  MI_cart_reg[z,] <- criteria_m(ana_cart)[6:14]
  
}
print(Sys.time()-strt)



##########   Ergebnisse durch aggregation über Z   ##################################################

### Ergebnisse zu mu4
set.seed(123)
Results <- matrix(nrow=5, ncol=11)
for (i in 1:5){
  Results[i,1] <- colMeans(BD)[i]
  Results[i,2] <- colMeans(CC)[i]
  Results[i,3] <- colMeans(SI_mean)[i]
  Results[i,4] <- colMeans(SI_hot)[i]
  Results[i,5] <- colMeans(SI_reg)[i]
  Results[i,6] <- colMeans(SI_norm)[i]
  Results[i,7] <- colMeans(SI_pmm)[i]
  Results[i,8] <- colMeans(SI_cart)[i]
  Results[i,9] <- colMeans(MI_norm)[i]
  Results[i,10] <- colMeans(MI_pmm)[i]
  Results[i,11] <- colMeans(MI_cart)[i]
}
Results <- round(Results, digits = 3)
rownames(Results) <- c("Bias","Coverage","CI_Width","Var[Mean(Y4)]","Cor(Y4,Y3)")
colnames(Results) <- c("BD","CC","SI_mean","SI_hot","SI_reg","SI_norm","SI_pmm","SI_cart","MI_norm","MI_pmm","MI_cart")
print(Results)

### Ergebnisse zu beta2,beta3
set.seed(123)
Results_Reg <- matrix(nrow=9, ncol=11)
for (i in 1:9){
  Results_Reg[i,1] <- colMeans(BD_reg)[i]
  Results_Reg[i,2] <- colMeans(CC_reg)[i]
  Results_Reg[i,3] <- colMeans(SI_mean_reg)[i]
  Results_Reg[i,4] <- colMeans(SI_hot_reg)[i]
  Results_Reg[i,5] <- colMeans(SI_reg_reg)[i]
  Results_Reg[i,6] <- colMeans(SI_norm_reg)[i]
  Results_Reg[i,7] <- colMeans(SI_pmm_reg)[i]
  Results_Reg[i,8] <- colMeans(SI_cart_reg)[i]
  Results_Reg[i,9] <- colMeans(MI_norm_reg)[i]
  Results_Reg[i,10] <- colMeans(MI_pmm_reg)[i]
  Results_Reg[i,11] <- colMeans(MI_cart_reg)[i]
}
Results_Reg <- round(Results_Reg, digits = 3)
rownames(Results_Reg) <- c("Bias_b2", "Coverage_b2", "CIwidth_b2", "SE_b2", "Bias_b3", "Coverage_b3", "CIwidth_b3", "SE_b3","R^2")
colnames(Results_Reg) <- c("BD", "CC", "SI_mean", "SI_hot", "SI_reg", "SI_norm", "SI_pmm", "SI_cart", "MI_norm", "MI_pmm", "MI_cart")
print(Results_Reg)

### Aggregierter Anteil der Ausfälle
colMeans(PD)[1]
colMeans(PD)[2]



##################################
##########   Grafiken   ##########
##################################

### Exemplarisches Datenausfallmuster
aggr(data_miss,numbers=TRUE,prop=TRUE,oma = c(10,5,5,3),ylabs =c("",""), combined =TRUE)
md.pattern(data_miss,plot=T,rotate.names = TRUE)

### Datenausfall in Y3,Y4
par(mfrow=c(1,2))
barMiss(x=data_miss, pos=3,miss.labels="NA")
barMiss(x=data_miss, pos=4,miss.labels="NA")
par(mfrow=c(1,1))

### MICE Konvergenz und Plot der Verteilung
set.seed(123)
Imp0 <- mice(data_miss, maxit = 20, printFlag = TRUE)
plot(Imp0)
densityplot(Imp0)

### Exemplarischer CART Plot
cart <- rpart(formula = Y4 ~ .,
              data    = data,
              method  = "anova")
rpart.plot(cart, type = 2)

### Plot der Mittelwertergänzungen
plot(Y4~Y3,col=rgb(0,0,1,0.3),pch=16, xlim = c(-1.5,6), ylim = c(1,40))
par(new=TRUE)
plot(Y4~Y3,data=data_meanImp,col=rgb(1,0,0,0.3),pch=16, xlim = c(-1.5,6),
     ylim = c(1,40), xlab ="", ylab= "")
legend(x=3.4,y=10,legend = c("tatsächliche Werte","nach Ergänzung"),fill=c("blue","red"))

### Plot der LR Ergänzung
LROrg <- lm(Y4~Y3+Y2+Y1)
LRImp <- lm(Y4~Y3+Y2+Y1, data = mice_reg)
summary(LROrg)
summary(LRImp)
# 3D scatterplot
Y2LRImp <- mice_reg$Y2
Y3LRImp <- mice_reg$Y3
Y4LRImp <- mice_reg$Y4
plot3d(Y4~Y3+Y2, col=rgb(0,0,1,0.1), xlim = c(-2,6),
       ylim=c(-2,6), zlim=c(0,21))
plot3d(Y2LRImp, Y3LRImp, Y4LRImp, col=rgb(1,0,0,0.2),add=TRUE)


