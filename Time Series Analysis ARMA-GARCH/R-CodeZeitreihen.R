#########Datensatz##########
#packages#
install.packages("tseries")
install.packages("urca")
install.packages("fGarch")
install.packages("zoo")
install.packages("dplyr")
install.packages("xts")
rm(list=ls())
library(tseries)
library(fGarch)
library(zoo)
library(dplyr)
library(xts)
library(urca)


#Daten einlesen#
data<-read.csv("S&P500daten.csv",header=TRUE)
head(data)
data<-data %>% select(Date,Adj.Close)

#Preis Plot#
par(mfrow=c(1,1))
data[,1]<-as.Date(data[,1])
data.plot<-xts(data[-1],order.by=data$Date)
plot.zoo(data.plot,ylab="Preis",xlab="Zeit",main="S&P Preis von Juli 2002 bis Ende Juni 2017")

#Renditen berechnen und plotten
data<- data %>% mutate(Return.Daily=log(Adj.Close/lag(Adj.Close)))
data<- data %>% select(Date,Return.Daily)
data<-na.omit(data)
data.plot2<-xts(data[-1],order.by=data$Date)
plot.zoo(data.plot2,ylab="Rendite",xlab="Zeit",main="S&P Rendite Juli 2002 bis Ende Juni 2017")
#quadrierte Renditen:
data.plot2<-xts(data[-1]^2,order.by=data$Date)
plot.zoo(data.plot2,ylab="quad. Rendite",xlab="Zeit",main="S&P quadrierte Rendite Juli 2002 bis Ende Juni 2017")
#kein Trends zu erkennen -> DF Tests ohne Trend

#Autokorrelationsfunktionen#
par(mfrow=c(2,2))
acf(data$Return.Daily,lag.max=200,main="Autokorrelationen der Returns")
acf(data$Return.Daily^2,lag.max=200,main="Autokorrelationen der R²")
pacf(data$Return.Daily,lag.max=200,main="PACF der Returns")
pacf(data$Return.Daily^2,lag.max=200,main="PACF der R²")

#Test und Training Period
data.train<-filter(data, Date < "2012-07-01")
data.test<-filter(data, Date >= "2012-07-01")
#########Modellierung und Evaluation##########

#Dickey Fuller Unit Root Test
summary(ur.df(data.train$Return.Daily,type = "none",selectlags = "AIC"))
summary(ur.df(data.train$Return.Daily^2,type = "none",selectlags = "AIC"))
#Weder in ARMAs für Renditen noch für Varianz liegen Einheitswurzeln vor

#AR#
AICs<- c(1:10)
  for (i in 1:10) {
    model<-arma(data.train$Return.Daily, order = c(i, 0))
    AICs[i] <- summary(model)$aic
  }
AICs
which.min(AICs)
#geringstes AIC: AR(7); sehr viele Koeffizienten insignifk, und  Konfidenzintervalle ohnehinn zu klein...
AR7<-arma(data.train$Return.Daily, order = c(7, 0))
AR1<-arma(data.train$Return.Daily, order = c(1, 0))
summary(AR7)
summary(AR1) #minimal größeres AIC aber signifikantere Koeff.

########################################################################GARCH###
#LM Test für GARCH/ARCH Effekte#
udach <- AR7$residuals
LMtest <- lm(udach ~ lag(udach,1) + lag(udach,2) + lag(udach,3) + lag(udach,4) + lag(udach,5) + lag(udach,6) + lag(udach,7) + lag(udach,8) + lag(udach,9) + lag(udach,10) +lag(udach,11) + lag(udach,12) + lag(udach,13) + lag(udach,14) + lag(udach,15) + lag(udach,16))
summary(LMtest)
udach <- AR1$residuals
LMtest <- lm(udach ~ lag(udach,1) + lag(udach,2) + lag(udach,3) + lag(udach,4) + lag(udach,5) + lag(udach,6) + lag(udach,7) + lag(udach,8) + lag(udach,9) + lag(udach,10) +lag(udach,11) + lag(udach,12) + lag(udach,13) + lag(udach,14) + lag(udach,15) + lag(udach,16))
summary(LMtest)

#GARCH Modelle für verschiedene AR#
#Kein AR Modell "AR(0)"
results0<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~garch(p,q), list(p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results0[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results0==min(results0),arr.ind = TRUE)
#GARCH(2,1) am besten nach AIC (erste Spalte: GARCH(i,0), zweite: GARCH(i,1))
#AR1
results1<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(1,0)+garch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results1[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results1==min(results1),arr.ind = TRUE)
#GARCH(2,1) am besten: 
#AR2
results2<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(2,0)+garch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results2[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results2==min(results2),arr.ind = TRUE)
#GARCH(2,1) am besten:
#AR3
results3<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(3,0)+garch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results3[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results3==min(results3),arr.ind = TRUE)
#GARCH(2,1) am besten:
#AR4
results4<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(4,0)+garch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results4[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results4==min(results4),arr.ind = TRUE)
#GARCH(4,1) am besten:
#AR5
results5<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(5,0)+garch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results5[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results5==min(results5),arr.ind = TRUE)
#GARCH(2,1) am besten:
#AR6
results6<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(6,0)+garch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results6[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results6==min(results6),arr.ind = TRUE)
results6
#GARCH(2,1) am besten;
#AR7
results7<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
for (i in 1:5) {
model<-garchFit(substitute(~arma(7,0)+garch(p,q), list(p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
results7[i,j+1]<-as.numeric(model@fit$ics[1])
}}
which(results7==min(results7),arr.ind = TRUE)
#GARCH(2,1) am besten nach AIC 
##################AR0-AR7 Vergleich mit jeweils besten GARCH
garchAIC <- vector("integer",8)
for (i in 0:7) {
  if(i==5) {
    model <-garchFit(~arma(5,0)+garch(4,1),data=data.train$Return.Daily, trace = FALSE)
    garchAIC[i+1]<-as.numeric(model@fit$ics[1])
   }
  if(i==0) {
    model <-garchFit(~garch(2,1),data=data.train$Return.Daily, trace = FALSE)
    garchAIC[i+1]<-as.numeric(model@fit$ics[1])
   } 
  if(i!=0&i!=5) {
    model <-garchFit(substitute(~arma(p,0)+garch(2,1), list(p=i)),data=data.train$Return.Daily, trace = FALSE)
    garchAIC[i+1]<-as.numeric(model@fit$ics[1])
    }
}
which.min(garchAIC)
#AR(6)+GARCH(2,1) nach AIC bestes Modell
garch <- garchFit(~arma(6,0)+garch(2,1),data=data.train$Return.Daily, trace = FALSE)
#############################################################APARCH######

#APARCH Modelle für verschiedene AR#
#Kein AR Modell "AR(0)"
results0<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~aparch(p,q), list(p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results0[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results0==min(results0),arr.ind = TRUE)
#aparch(2,1) am besten: 
#AR1
results1<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(1,0)+aparch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results1[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results1==min(results1),arr.ind = TRUE)
#aparch(2,1) am besten: 
#AR2
results2<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(2,0)+aparch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results2[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results2==min(results2),arr.ind = TRUE)
#aparch(2,1) am besten:
#AR3
results3<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(3,0)+aparch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results3[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results3==min(results3),arr.ind = TRUE)
#aparch(2,1) am besten:
#AR4
results4<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(4,0)+aparch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results4[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results4==min(results4),arr.ind = TRUE)
#aparch(2,1) am besten:
#AR5
results5<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(5,0)+aparch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results5[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results5==min(results5),arr.ind = TRUE)
#aparch(2,2) am besten:
#AR6
results6<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(6,0)+aparch(p,q), list( p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results6[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results6==min(results6),arr.ind = TRUE)
results6
#aparch(2,2) am besten;
#AR7
results7<-matrix(0,nrow=5,ncol=6)
for (j in 0:5){
  for (i in 1:5) {
    model<-garchFit(substitute(~arma(7,0)+aparch(p,q), list(p=i, q=j)),data=data.train$Return.Daily, trace = FALSE)
    results7[i,j+1]<-as.numeric(model@fit$ics[1])
  }}
which(results7==min(results7),arr.ind = TRUE)
#aparch(2,2) am besten nach AIC
##################AR0-AR7 Vergleich mit jeweils besten ARCH
aparchAIC <- vector("integer",8)
for (i in 0:7) {
  if(i>0&i<5) {
    model <-garchFit(~arma(5,0)+aparch(2,1),data=data.train$Return.Daily, trace = FALSE)
    aparchAIC[i+1]<-as.numeric(model@fit$ics[1])
  }
  if(i==0) {
    model <-garchFit(~aparch(2,1),data=data.train$Return.Daily, trace = FALSE)
    aparchAIC[i+1]<-as.numeric(model@fit$ics[1])
  } 
  if(i>4) {
    model <-garchFit(substitute(~arma(p,0)+aparch(2,2), list(p=i)),data=data.train$Return.Daily, trace = FALSE)
    aparchAIC[i+1]<-as.numeric(model@fit$ics[1])
  }
}
which.min(aparchAIC)
#bestes APARCH nach AIC: AR(6)+APARCH(2,1)
aparch <-garchFit(~arma(6,0)+aparch(2,1),data=data.train$Return.Daily, trace = FALSE)
################Vergleich der besten Modelle#########################################
summary(garch)
plot(garch@residuals)
plot(aparch@fitted)
plot(AR7$residuals)
summary(aparch)
summary(AR7)
#geringstes AIC:
min(AICs)
min(garchAIC)
min(aparchAIC)
#größtes R^2:
RQuad <- c(0,0,0)
RQuad[1] <- 1 - var(AR7$residuals, na.rm=TRUE)/var(data.train$Return.Daily)
RQuad[2] <- 1 - var(garch@residuals)/var(data.train$Return.Daily)
RQuad[3] <- 1 - var(aparch@residuals)/var(data.train$Return.Daily)
RQuad
#Ljung Box Fehler
Box.test(AR7$residuals, lag=15, type="Ljung-Box", fitdf=7)
Box.test(AR1$residuals, lag=15, type="Ljung-Box", fitdf=1)
Box.test(AR7$residuals, lag=20, type="Ljung-Box", fitdf=7)
Box.test(AR1$residuals, lag=20, type="Ljung-Box", fitdf=1)
#Jarque Bera Normalverteilung der Fehler
jarque.bera.test(removeNA(AR1$residuals))
jarque.bera.test(removeNA(AR7$residuals))
jarque.bera.test(garch@residuals/sqrt(garch@h.t))
jarque.bera.test(aparch@residuals/sqrt(aparch@h.t))
#Histogramme
par(mfrow=c(2,2))
hist(removeNA(AR1$residuals))
hist(removeNA(AR7$residuals))
hist(garch@residuals/sqrt(garch@h.t))
hist(aparch@residuals/sqrt(aparch@h.t))
head(AR1$residuals)
plot(data.train$Date[-1],removeNA(AR1$residuals), ylab="Residuen",xlab="Zeit", main="AR(1)")
plot(data.train$Date[8:length(data.train$Date)],removeNA(AR7$residuals), ylab="Residuen",xlab="Zeit", main="AR(7)")
plot(data.train$Date,garch@residuals/sqrt(garch@h.t), ylab="standardisierte Residuen",xlab="Zeit", main="AR(6)+GARCH(2,1)")
plot(data.train$Date,aparch@residuals/sqrt(aparch@h.t), ylab="standardisierte Residuen",xlab="Zeit",main="AR(6)+APARCH(2,1)")
which.min(garch@residuals/sqrt(garch@h.t))
jarque.bera.test(garch@residuals[-1172]/sqrt(garch@h.t[-1172]))
jarque.bera.test(aparch@residuals[-1172]/sqrt(aparch@h.t[-1172]))
plot(data.train$Date,rnorm(2519,0,mean((garch@residuals/sqrt(garch@h.t))^2)),ylab="simulierte Residuen",xlab="Zeit",main="Simulierte, normalverteilte Residuen")

###########################################################MSEs##################################
#AR1
estimate0 <-0
for(i in 2:(length(data.test$Return.Daily)-1))
{
  estimate0[i+1]<-AR1$coef[2]+AR1$coef[1]*data.test$Return.Daily[i]
}
MSEAR1 <- sum((data.test$Return.Daily[3:length(data.test$Return.Daily)]-estimate0[3:length(estimate0)])^2)*1/(length(data.test$Return.Daily)-3)
#AR7
estimate1 <-0
for(i in 8:(length(data.test$Return.Daily)-1))
{
  estimate1[i+1]<-AR7$coef[8]+AR7$coef[1]*data.test$Return.Daily[i]+AR7$coef[2]*data.test$Return.Daily[i-1]+AR7$coef[3]*data.test$Return.Daily[i-2]+AR7$coef[4]*data.test$Return.Daily[i-3]+AR7$coef[5]*data.test$Return.Daily[i-4]+AR7$coef[6]*data.test$Return.Daily[i-5]+AR7$coef[7]*data.test$Return.Daily[i-6]
}
MSEAR7<-sum((data.test$Return.Daily[9:length(estimate1)]-estimate1[9:length(estimate1)])^2)*1/(length(estimate1)-9)
#AR6+GARCH
garch@fit$coef[7]
estimate2 <-0
for(i in 7:(length(data.test$Return.Daily)-1))
{
  estimate2[i+1]<-garch@fit$coef[1]+garch@fit$coef[2]*data.test$Return.Daily[i]+garch@fit$coef[3]*data.test$Return.Daily[i-2]+garch@fit$coef[4]*data.test$Return.Daily[i-3]+garch@fit$coef[5]*data.test$Return.Daily[i-4]+garch@fit$coef[6]*data.test$Return.Daily[i-5]+garch@fit$coef[7]*data.test$Return.Daily[i-6]
}
MSEAR6<-sum((data.test$Return.Daily[8:length(estimate2)]-estimate2[8:length(estimate2)])^2)*1/(length(estimate2)-8)
#MSEs:
MSEAR6
MSEAR1
MSEAR7

