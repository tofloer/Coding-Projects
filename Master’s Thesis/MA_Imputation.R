# Tobias Floerchinger (608227), floercht@hu-berlin.de


# R-Code of Master Thesis: 'Predicting Interdependent Diseases in Individuals 
#                           Using Survey Data and Machine Learning Algorithms'


# In order to reproduce the MA results, this script requires:
#   - R version 4.1.2 from (https://cran.r-project.org/bin/windows/base/old/)
#   - RStudio version 2021.9.2.382 from (https://www.rstudio.com/products/rstudio/older-versions/#2021092)


################################################################################
############################   Setup Environment   #############################
################################################################################

### Packages ###
packages =  c('rstudioapi',
              'haven',
              'dplyr',
              'mice')

### Install package if not available ###
inst = packages %in% installed.packages()    
if(length(packages[!inst]) > 0) install.packages(packages[!inst])

### Load packages in environment ###
lapply(packages, require, character.only=TRUE)

### Set WD to the R script directory ###
path.wd = dirname(getSourceEditorContext()$path)
setwd(path.wd)



################################################################################
#########################  Imputation for Values MAR   #########################
################################################################################

### Read data from directory of R script ###
X = read.csv('X.csv')

### Replace NULL values with NA ###
X$X = NULL
X[X == ''] = NA

### Change type of character variables to factor ###
tofactor <- colnames(X[, sapply(X, class) == 'character'])
X[tofactor] <- lapply(X[tofactor], factor)

### Initialize MICE prediction matrix ###
ini <- mice(X,m=1, maxit=0)
pred <- ini$pred

### Imputation via CART and 10 iterations ###
set.seed(123)
imp <- mice(X, m=1, maxit=10, predictorMatrix = pred, print=F, method = "cart")
X_imp <- complete(imp, include=FALSE)

### Check for missing values and log events ###
sapply(X_imp, function(x) round((sum(is.na(x))/nrow(X_imp))*100,2) )
imp$loggedEvents$out

### Write imputed data to .csv file ###
write.csv(x=X_imp, file='X_imp.csv')

