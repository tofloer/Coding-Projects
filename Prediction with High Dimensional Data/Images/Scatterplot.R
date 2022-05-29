if (!require(MASS)) install.packages('MASS')
library(MASS)

data("Boston")
Boston <- Boston[,-c(4,9)]

#Scatterplot
dev.new(width=6, height=5)
par(mfrow=c(4,3))

for (i in 1:11) {
  x <- Boston[,i]
  plot(Boston$medv, x, xlab="medv", ylab=colnames(Boston)[i], pch=20, cex.lab=1.5)
  #abline(lm(Boston$medv~x), col="red") # regression line (y~x)
  lines(lowess(Boston$medv, x),lty="solid" , col="red", lwd=2) # lowess line (x,y)
}