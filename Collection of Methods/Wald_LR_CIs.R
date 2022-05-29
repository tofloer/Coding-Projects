####################################################################################
##3. HEALTH INSURANCE --> YES/NO?
####################################################################################
#VietnamI- Expenses
#Source: http://microdata.worldbank.org/index.php/catalog/2694/study-description


#Load Data & Packages

install.packages("mosaic") #needed
install.packages("proportion") #needed
install.packages("Ecdat") 
install.packages("xtable")#needed
install.packages("emplik")

library("emplik")
library("Ecdat")
library("mosaic")
library("xtable")
library("proportion")
data(VietNamI)
View(VietNamI)

alpha = 0.01
n = nrow(VietNamI)

#Estimate of the relative frequency of health insured persons in the Vietnamese population
p = mean(VietNamI$insurance)  #estimate of the relative frequency of insured persons

#######################################################################################
## 3.1 Wald Confidence interval
#######################################################################################

CIinfB   = p - qnorm(1-alpha/2)*sqrt(p*(1-p)/n) 
CIsupB   = p + qnorm(1-alpha/2)*sqrt(p*(1-p)/n)

#exact CI:
#binom.test(p*n,n,alternative = "two.sided", conf.level=1-alpha)

#####################
#likelihoodratio
#https://cran.r-project.org/web/packages/proportion/proportion.pdf

#ciAAllx(n*p, n, 0.01, 0)
#####################################################################################
##3.2 Likelihood ratio Confidence interval
####################################################################################

###Programming function
# see http://www.ms.uky.edu/~mai/sta635/Two%20Examples.pdf


#Calculating -2*log-likelihood ratio ( which is Chi2 distributed)
LogLikRatio = function(n,x,p){
  phat = x/n
  temp= -2*(x*log(p/phat)+(n-x)*log((1-p)/(1-phat)))
  list("-2LLR"=temp)
}

#Solving the equation: -2*Log-Likelihood-ratio = chi2_quantile(1-alpha)
CI= findUL( step=0.00000001, fun = LogLikRatio, MLE=p, n=n, x=round(p*n), level= qchisq(1-alpha, df=1))


#Export results in table

CIBinTab = data.frame("p"=p,"CI-Inf"= CIinfB,"CI-Sup" = CIsupB )

CIBinTab = rbind(CIBinTab, c(p, CI$Low, CI$Up))

row.names(CIBinTab) = c("Wald", "Likelihood-Ratio")

print(xtable(CIBinTab,
             caption="Estimation and confidence interval for the relative frequency of insured persons in the Vietnamese population in 1997", 
             type = "latex", file = "CBinomial.tex", digits= 5),
      include.rownames=TRUE)