###############################################################################
##                            Bootstrap: Normal                              ##
###############################################################################

### Setting ###

# Number of simulations
B <- 2000 
# Alpha level
alpha <- .05 
# Sample size
n <- 50 
# mu
mu <- 0


### Drawing one data set, estimating mu once ###

# Raw data
raw <- rnorm(n,mu, 1) 
# mu.hat
mu_hat <- mean(raw)


### Bootstrap sample ###

# Drawing the samples
boot_samples <- matrix(sample(raw, size = B * n, replace = TRUE), B, n)
# one example bootstrap sample
sort(boot_samples[1,])
# Estimated mean in each bootstrap sample
boot_mean <- apply(boot_samples, 1, mean)


### Estimating bias and variance ###

# Calculating standard deviations
sd_boot <- sd(boot_mean)
sd_raw <- sd(raw)

# The bias of the sample mean as an estimator is zero
# The bootstrap sample will on average replicate the empirical bias in the original sample due to finite n
bias_raw <- mu_hat - mu
bias_boot <- mean(boot_mean) - mu


### Confidence intervals  ###

# normal quantile 
z_1alpha2 <- qnorm(1 - alpha / 2)

# Normal confidence interval from original sample
ci_raw <- c(mu_hat - z_1alpha2 * sd_raw / sqrt(n),
            mu_hat + z_1alpha2 * sd_raw / sqrt(n) )

# Bootstrap percentile confidence interval
ci_percentile <- quantile(boot_mean, probs = c(alpha/2, 1 - alpha/2))

# Asymptotic bootstrap confidence interval
ci_asymptotic <- c(mu_hat - z_1alpha2 * sd_boot, 
                   mu_hat + z_1alpha2 * sd_boot) 

# Bootstrap t confidence interval
zb_star <- boot_mean - mu_hat
zb_star <- zb_star / sd(zb_star)
z_low  <- quantile(zb_star, probs = 1 - alpha/2)
z_up  <- quantile(zb_star, probs = alpha/2)
ci_t <- c(mu_hat - z_low * sd_boot,
          mu_hat - z_up  * sd_boot)


### Plotting the results ###
#par(mfrow = c(1, 2))

# Sample
hist(raw, nclass = 50, col = 3, xlab = "observations",
     main = "Histogram of the True Sample Values")
abline(v = ci_raw, col = "turquoise")
abline(v = 0, col = "black", lty = 2)
abline(v = mu_hat, col = "pink", lty = 2)

# Bootstrap
hist(boot_mean, nclass = 50, col = 3, xlab = "Bootstrap means",
     main = "Histogram of Bootstrap Means")
abline(v = ci_raw, col = "turquoise")
abline(v = ci_percentile, col = "red")
abline(v = ci_asymptotic, col = "blue")
abline(v = ci_t, col = "green")
abline(v = 0, col = "black", lty = 2)
abline(v = mu_hat, col = "pink", lty = 2)
legend("topright", 
       legend = c("CI original sample", "Bootstrap percentile", "Bootstrap asymptotic", "Bootstrap t", "true mean", "estimate in sample"),
       col = c("turquoise", "red", "blue", "green", "black", "pink"), lty = c(1,1,1,1,2,2))
