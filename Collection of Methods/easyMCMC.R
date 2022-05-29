# https://stephens999.github.io/fiveMinuteStats/MH-examples1.html # November 2020

# target distribution
target = function(x){
  if(x<0){
    return(0)
    } else {
    return( exp(-x))
  }
}

easyMCMC = function(niter, startval, proposalsd){
  x = rep(0,niter)
  x[1] = startval     
  for(i in 2:niter){
    currentx = x[i-1]
    proposedx = rnorm(1,mean=currentx,sd=proposalsd) 
    A = target(proposedx)/target(currentx)
    if(runif(1)<A){
      x[i] = proposedx       # accept move with probabily min(1,A)
    } else {
      x[i] = currentx        # otherwise "reject" move, and stay where we are
    }
  }
  return(x)
}

z1=easyMCMC(100000,30,1)
z2=easyMCMC(100000,30,2)
z3=easyMCMC(100000,30,3)

plot(z1[1:300],type="l")
lines(z2[1:300],col=2)
lines(z3[1:300],col=3)

par(mfrow=c(2,1))
plot(x=seq(0,8,by=0.01),y=dexp(seq(0,8,by=0.01),rate=1),col="brown",lwd=3)
lines(density(z1,from=0),col=1,lwd=2)
lines(density(z2,from=0),col=2,lwd=2)
lines(density(z3,from=0),col=3,lwd=2)

plot(x=seq(0,8,by=0.01),y=dexp(seq(0,8,by=0.01),rate=1),col="brown",lwd=3)
lines(density(z1[-c(1:5000)],from=0),col=1,lwd=2)
lines(density(z2[-c(1:5000)],from=0),col=2,lwd=2)
lines(density(z3[-c(1:5000)],from=0),col=3,lwd=2)


# problem with density at the left can be handeled by using logspline(z1, lbound = 0)

# See also:
# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
# https://bookdown.org/rdpeng/advstatcomp/metropolis-hastings.html
# http://www.mas.ncl.ac.uk/~ndjw1/teaching/sim/metrop/

