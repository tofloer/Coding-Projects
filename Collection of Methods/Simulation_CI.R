# select the parameters for the simulation, set-up
n = 100 # sample size
p = 0.1 # success probability
alpha = 0.05 # 1 - confidence level
nsim = 100 # number of simulations
z <- qnorm(1 - alpha / 2)
contained <- rep(NA, nsim) # indicator vector, for which simulations the confidence interval contains the true value pi

# simulate nsim samples of size n i.i.d. Bernoulli(p), compute estimate and confidence interval for each simulation
x_bar <- replicate(nsim, mean(rbinom(n, 1, p))) # \hat\pi
lower_bound <- x_bar - z * sqrt((x_bar * (1 - x_bar))/n) # lower
upper_bound <- x_bar + z * sqrt((x_bar * (1 - x_bar))/n) # and upper bounds of the confidence intervals
contained <- (lower_bound < p & p < upper_bound)
mean(contained) # empirical coverage in the nsim simulations, should be close to the nominal 1-alpha

# plot for each simulation the estimate and the confidence interval. Mark CIs that do not cover pi in red.
plot(1:nsim, x_bar, ylim = c(min(lower_bound), max(upper_bound)), xlab = "iterations", ylab = "mean and CI", pch=19, cex=0.5, col=4)
abline(h = p, col = 3)
arrows(x0 = 1:nsim, y0 = lower_bound, y1 = upper_bound, col = 1, angle = 90, code = 3, length  = 0.05)
arrows(x0 = (1:nsim)[!contained], y0 = lower_bound[!contained], y1 = upper_bound[!contained], col = 2, angle = 90, code = 3, length  = 0.05)