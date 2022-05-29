library("Ecdat")
data(VietNamI)

### divide data by insurance status, preliminaries

VietNamI0 = VietNamI[VietNamI$insurance==FALSE,]
n0 = length(VietNamI0$lnhhexp)
MEANExp0 = mean(VietNamI0$lnhhexp)
SDExp0 = sd(VietNamI0$lnhhexp)

VietNamI1 = VietNamI[VietNamI$insurance==TRUE,]
n1 = length(VietNamI1$lnhhexp)
MEANExp1 = mean(VietNamI1$lnhhexp)
SDExp1 = sd(VietNamI1$lnhhexp)

##### 1D Delta method for the median ###

## median
exp(MEANExp0)
## sd, asymptotic sd
sqrt( (exp(SDExp0^2 / n0) - 1) * exp(2 * MEANExp0 + SDExp0^2 / n0) )
sqrt( exp(2 * MEANExp0) * SDExp0^2 / n0 )

## median
exp(MEANExp1)
## sd, asymptotic sd
sqrt( (exp(SDExp1^2 / n1) - 1) * exp(2 * MEANExp1 + SDExp1^2 / n1) )
sqrt( exp(2 * MEANExp1) * SDExp1^2 / n1 )


##### 2D Delta method for the mean ###

## mean
exp(MEANExp0 + SDExp0^2 /2)
## asymptotic sd
sqrt( exp(2 * MEANExp0 + SDExp0^2) * (SDExp0^2 / n0 + SDExp0^2 / (2 * n0) ) )


## mean
exp(MEANExp1 + SDExp1^2 /2)
## asymptotic sd
sqrt( exp(2 * MEANExp1 + SDExp1^2) * (SDExp1^2 / n1 + SDExp1^2 / (2 * n1) ) )