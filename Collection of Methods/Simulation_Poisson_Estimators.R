#Simulation of poisson estimators Lambda_1, Lambda_2, Lambda_3, Lambda_4

#install.packages("xtable")
library(xtable)
#init
Lambda = 6
Sample_Length = 200
nsamples= 10000
Samples = list()

lambda_1_vec = rep(NA, length = nsamples)
lambda_2_vec = rep(NA, length = nsamples)
lambda_3_vec = rep(NA, length = nsamples)
lambda_4_vec = rep(NA, length = nsamples)
Ones = rep(1, Sample_Length )


#generate "nsamples" with length "Sample_Length"
for (i in 1:nsamples) {  
  set.seed(i)   #taking a new random sequence
  Samples[[i]] = rpois(Sample_Length,Lambda)
  
  #estimating lambdas for each sample
  lambda_1_vec[i] = mean(Samples[[i]])
  lambda_2_vec[i] = 0.5*(min(Samples[[i]])+max(Samples[[i]]))
  lambda_3_vec[i] = 0.5*(Samples[[i]][1]+Samples[[i]][Sample_Length])
  lambda_4_vec[i] = (prod(Samples[[i]]+Ones))^(1/Sample_Length)
}

#computing the mean over all samples
lambda_1 = mean(lambda_1_vec)
lambda_2 = mean(lambda_2_vec)
lambda_3 = mean(lambda_3_vec)
lambda_4 = mean(lambda_4_vec)

#computing the SD of each estimator 
sd1 = sd(lambda_1_vec)
sd2 = sd(lambda_2_vec)
sd3 = sd(lambda_3_vec)
sd4 = sd(lambda_4_vec)



#Export results in table
Results = data.frame("lambda1"=lambda_1,
                     "lambda2"=lambda_2,
                     "lambda3"=lambda_3, 
                     "lambda4"=lambda_4 )

Results = rbind(Results, c(sd1, sd2, sd3, sd4))

row.names(Results)= c("Mean", "Standard <deviation")

print(xtable(Results,
             caption="Simulation of Mean-Estimators for a Poisson distribution with a mean of 6 (based on 1000 random samples with length 200)", 
             type = "latex", file = "PoissonEstimators.tex", digits= 5))