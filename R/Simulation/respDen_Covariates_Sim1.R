### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
##       Simulation of response densities with Covariates   ##
##               Author: Johannes                           ##
##                      3.12.2020                           ##
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###


# Data Simulation ---------------------------------------------------------

## bayes-space linear model with beta-distribution-intercept
# reparametrized distribution function for beta distribution such that shape1-1 = a and shape2-1 = b
dbeta2 <- function(x, a, b, ncp = 0, log = FALSE) dbeta(x, a+1, b+1, ncp, log)
rbeta2 <- function(x, a, b, ncp = 0) rbeta(x, a+1, b+1, ncp)
# intercept parameters
a0 <- 2
b0 <- 2
# time grid
t <- seq(0,1, length.out = 50)
# plot intercept distribution
plot(t, dbeta2(t, a0, b0), type = "l")
lines(t, dbeta(t, 3,3), type = "l", col= "red")
## specify slope
slope1 <- 1
slope2 <- 2
slope3 <- 1
# sample size
N <- 60


# Covariates --------------------------------------------------------------
## sample covariates >-1/(slope*min(a0,b0)) from exponential distribution substracted by the same constant
set.seed(493822)
x1 <- rexp(N)-1/(slope1*min(a0, b0))
range(x1)
# we want x to linearly decrease the variance, without shifting the expectation,
# i.e. x increases density in center and decreases the tails
x2 <- rep(seq(0,2, by = 1), N/3)
# we want x2 to be a factor that shifts the density to the right through increasing 
# paramter a
x3 <- runif(N,1,3)
# x3 is a continuous covariate which shifts the distribution to the left


# Random Effects ----------------------------------------------------------
# Number of levels
G <- 10
b_0 <- rnorm(G, 0.5)
group <- rep(b_0, N/G)

# Linear Effects ----------------------------------------------------------


a <- slope1*x1*a0 + slope2*x2 + group
range(a)
hist(a)
b <- slope1*x1*b0 + slope3*x3 + group
range(b)
hist(b)
## recode factors
x2 <- as.factor(x2)
group <- as.factor(group)

## center continous covarites
x1 <- as.vector(scale(x1, center = T, scale = F))
x3 <- as.vector(scale(x3, center = T, scale = F))

# Sample from response distributions --------------------------------------
# plot true distribution functions
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l")
# there are a few strange distributions but they should not worry us

## sample response curves in a realistic procedure by sampling from the corresponding distributions
# number of samples per density function
n <- 500
nsamples <- apply(cbind(a,b), 1, function(p) rbeta2(n, p[1], p[2]))
# estimate densities from sample
density_esimates <- structure(apply(nsamples, 2, density, bw = .05, n = length(t), from = 0, to = 1))
y <- sapply(density_esimates, function(x) x$y)
# plot estimated densities
matplot(t, y, type = "l")

## apply the clr-transform of Talska et al. 2017 to get standard L^2 functions
yclr <- scale( log(y), center = TRUE, scale = FALSE )

# check clr-constraint 
# each column is one density
stopifnot(all(colMeans(yclr)<1e-10))
matplot(t, yclr, type = "l")



# Data Format Array and long -------------------------------------------------------

## Prepare data for FDboost and refund ## 
## Array format
dat_ar <- list(yclr = t(yclr), t = t, x1 = x1, x2 = x2, x3 = x3, group = group)

## Long format
dat_lon <- data.frame(yclr = as.vector(yclr), 
                     t = rep(t,N), x1 = rep(x1, each = length(t)), 
                     x2 = rep(x2, each = length(t)), x3 = rep(x3, each = length(t)),
                     group = rep(group, each = length(t)),
                     id = rep(factor(1:N), each = length(t)))


# FDboost -----------------------------------------------------------------

# 1. Model only with Almonds Covariate
fd_mod1 <- FDboost(yclr ~ 1 + bolsc(x1, df = 2), timeformula = ~bbsc(t, df = 2), data = dat_ar)
plot(fd_mod1)
# offset seems to extreme (should translate into beta(3,3))

# 2. Model with all Covariates
fd_mod2 <- FDboost(yclr ~ 1 + bolsc(x1, df = 2) + bolsc(x2, df = 2) + bolsc(x3, df = 2), timeformula = ~bbsc(t, df = 2), data = dat_ar)
plot(fd_mod2)


# Predictions -------------------------------------------------------------

yclrpred1 <- t(predict(fd_mod2))
# we check iff  all clr-densities sum up to zero!
stopifnot(all(colMeans(yclrpred1)<1e-10))
# apply the reverse clr-transformation to obtain density data
ypred1 <- sweep(exp(yclrpred1), 2, colMeans(exp(yclrpred1)), "/")
# check if densities sum up to one
stopifnot(all(abs(colMeans(ypred1)-1)<1e-10))

## compare true curves, data curves, and predicted curves
opar <- par(mfrow = c(1,3))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l", main = "true underlying beta-densities")
matplot(t, ypred1, type = "l", main = "predicted densities in linear model")
matplot(t, y, type = "l", main = "data consisting of kernel-density-estimates")
par(opar)

# predict with covariate effects on a0 and b0 set to zero
newdat <- dat_ar
# add one row
newdat$x1[1] <- 1/slope1
newdat$x2[1] <- "0"
newdat$x3[1] <- 0
a0b0estimateclr1 <- predict(fd_mod2, newdata = newdat)[1,]
# transform back into bayes space
a0b0estimate1 <- exp(a0b0estimateclr1)/mean(exp(a0b0estimateclr1))

# in fact a beta(3,3) distribution
plot(t, dbeta2(t, a0, b0), t = "l", ylim =  lty = 1, main = "True basis density of beta(a0,b0)-distribution \n and predicted density by linear model")
lines(t, a0b0estimate1, col = "cornflowerblue", lty = 2)
legend("bottomright", legend = c("true", "predicted"), col = c("black", "cornflowerblue"), lty = 1:2)
# that could be related to the group factor

## for Almonds model
a0b0estimateclr2 <- predict(fd_mod1, newdata = newdat)[1,]
# transform back into bayes space
a0b0estimate2 <- exp(a0b0estimateclr2)/mean(exp(a0b0estimateclr2))

# in fact a beta(3,3) distribution
plot(t, dbeta2(t, a0, b0), t = "l", ylim = c(0,2),  lty = 1, main = "True basis density of beta(a0,b0)-distribution \n and predicted density by linear model")
lines(t, a0b0estimate2, col = "cornflowerblue", lty = 2)
legend("bottomright", legend = c("true", "predicted"), col = c("black", "cornflowerblue"), lty = 1:2)
# better but still disturbed, which is a hint to the group variable

## is it better with smooth effects?
## fit flexible model with smooth base learners
fd_mod3 <- FDboost(yclr ~ 1+ bbsc(x1, df = 2), timeformula = ~bbsc(t, df = 2), data = dat_ar) # set mstop to 500 for letting the fitting algo run longer
fd_mod4 <- FDboost(yclr ~ 1+ bbsc(x1, df = 2) + bolsc(x2, df = 2) + bbsc(x3, df = 2), timeformula = ~bbsc(t, df = 2), data = dat_ar) # set mstop to 500 for letting the fitting algo run longer
plot(fd_mod4)
# looks more smooth (haha)

## compare estimated intercept with true undelying intercept 
# (or rather the curve with the parameters specified as the intercept above)
# produce new data with one observation at x = 1/slope corresponding to a0 and b0

a0b0estimateclr3 <- predict(fd_mod4, newdata = newdat)[1,]
a0b0estimate3 <- exp(a0b0estimateclr3)/mean(exp(a0b0estimateclr3))

plot(t, dbeta2(t, a0, b0), t = "l", lty = 1, main = "True basis density of beta(a0,b0)-distribution \n and predicted densities")
lines(t, a0b0estimate1, col = "cornflowerblue", lty = 2)
lines(t, a0b0estimate3, col = "darkseagreen", lty = 3)
legend("bottomright", legend = c("true", "linear model", "flexible model"), col = c("black", "cornflowerblue", "darkseagreen"), lty = 1:3)
## indeed smooth model is better in dealing with those group effects


# refund/mgcv -------------------------------------------------------------


library(refund)
## Not clear how to implement the sum-to-zero constraint also for the development over time in pffr
re_mod1 <- refund::pffr( yclr ~ 1 + s(x1, bs = "ps"), yind = t, data = dat_ar, bs.yindex = list(mc = c(TRUE,TRUE)))
plot(re_mod1)
# smooth intercept is equal to FDboost, x1 effect seems off 

## Let's do a full model with random effects
re_mod2 <- refund::pffr( yclr ~ 0 + c(0) + group + c(group) + s(x1, bs = "ps") + s(x2, bs = "re") + s(x3, bs = "ps") , yind = t, data = dat_ar)
# this model gives us a random intercept for each group; problem is that the general random intercept is not corectly specified
# but now x1 looks much better
# x3 is fuzzy but meets expectation
# x2?? -> check coefficient vector and estimate at hand -> what we should do anyway
re_mod2b <- refund::pffr( yclr ~ 1 + group + c(group) + s(x1, bs = "ps") + s(x2, bs = "re") + s(x3, bs = "ps") , yind = t, data = dat_ar)
re_mod2c <- refund::pffr( yclr ~ 1 + c(group) + s(x1, bs = "ps") + s(x2, bs = "re") + s(x3, bs = "ps") , yind = t, data = dat_ar)

re_mod2$formula

yclrpred3 <- t(predict(re_mod2))
colMeans(yclrpred3)
# apply the reverse clr-transformation to obtain density data
ypred3 <- sweep(exp(yclrpred3), 2, colMeans(exp(yclrpred3)), "/")
stopifnot(all(abs(colMeans(ypred3)-1)<1e-10))
# densities are still fine??
## compare estimated with actual curves
opar <- par(mfrow = c(1,3))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l", main = "true underlying beta-densities")
matplot(t, ypred3, type = "l", main = "predicted densities in pffr with varying school effects")
matplot(t, y, type = "l", main = "data consisting of kernel-density-estimates")
par(opar)
## in fact it looks far to "wiggely" -> we could reduce smoothness

## Test
test <- t(predict(re_mod2))
stopifnot(max(abs(colSums(test))) < 1e-12)

## with scalar random intercept for school
re_mod3 <- refund::pffr( yclr ~ c(s(group, bs = "re")) + s(x1, bs = "ps") + s(x2, bs = "re") + s(x3, bs = "ps") , yind = t, data = dat_ar)
re_mod3b <- refund::pffr( yclr ~ 1 + c(s(group, bs = "re")) + s(x1, bs = "ps") + s(x2, bs = "re") + s(x3, bs = "ps") , yind = t, data = dat_ar)

plot(re_mod3)
re_mod3$pffr$termmap
# gives an "ugly" random intercept
# checks gaussian assumption of group level
# gives a correct plot for x1 and presumably for x3
## x2?
yclrpred4 <- t(predict(re_mod3))
colMeans(yclrpred4)
# apply the reverse clr-transformation to obtain density data
ypred4 <- sweep(exp(yclrpred4), 2, colMeans(exp(yclrpred4)), "/")
stopifnot(all(abs(colMeans(ypred4)-1)<1e-10))
# densities are still fine??
## compare estimated with actual curves
opar <- par(mfrow = c(1,3))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l", main = "true underlying beta-densities")
matplot(t, ypred4, type = "l", main = "predicted densities in pffr with scalar random intercept")
matplot(t, y, type = "l", main = "data consisting of kernel-density-estimates")
par(opar)
## that looks better even though "flatter" curves are still pretty wiggely

## Combine the two
# re_mod6 <- refund::pffr( yclr ~ 1 + x2 + c(x2) + c(s(group, bs = "re")) + s(x1, bs = "ps") + s(x3, bs = "ps") , yind = t, data = dat_ar)

re_mod6 <- refund::pffr( yclr ~ 1 + x2 + c(s(group, bs = "re")) + s(x1, bs = "ps") + s(x3, bs = "ps") , yind = t, data = dat_ar)
re_mod6$formula


yclrpred6 <- t(predict(re_mod6))
colMeans(yclrpred6)
# apply the reverse clr-transformation to obtain density data
ypred6 <- sweep(exp(yclrpred6), 2, colMeans(exp(yclrpred6)), "/")
stopifnot(all(abs(colMeans(ypred4)-1)<1e-10))
# densities are still fine??
## compare estimated with actual curves
opar <- par(mfrow = c(1,3))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l", main = "true underlying beta-densities")
matplot(t, ypred6, type = "l", main = "predicted densities in pffr Model6")
matplot(t, y, type = "l", main = "data consisting of kernel-density-estimates")
par(opar)


# Big question: what is the difference between the two?? ------------------


# MGCV --------------------------------------------------------------------

library(mgcv)

gam_mod1 <- gam(yclr ~ 0 +
                ti(t, bs = "ps", mc = TRUE, k = 8) + # functional intercept
                ti(x1, t, # covariate and time
                   bs = c("ps", "ps"), # use B-splines for both
                   mc = c(TRUE, TRUE), # sum-to-zero constraints for both x and t
                   k = c(8, 8)), # 8 basis functions for both x and t (before constraint)
              data = dat_lon) 
plot(gam_mod1)
# random intercept is closer to FDboost
# x1 is fuzzy
# check constraint
dat_lon$yclr_pred <- predict(gam_mod1)
yclrpred_gam <- array(dat_lon$yclr_pred, dim = dim(yclr))
stopifnot(max(abs(colSums(yclrpred_gam))) < 1e-12)

## re-gain and plot predicted densities
ypred_gam <- sweep(exp(yclrpred_gam), 2, colMeans(exp(yclrpred_gam)), "/")
matplot(ypred_gam, t = "l")
# Hmm, seems the model wants to center the data quite a lot, all deviant densities are just pushed
# into the middle. Since we only account for x1 and an intercept that is not completly surprising


## now let's try to implement the other
gam_mod2 <- gam(yclr ~ 0 +
                  ti(t, bs = "ps", mc = TRUE, k = 8) + # functional intercept
                  s(group, bs = "re") +
                  ti(x1, t, # covariate and time
                     bs = c("ps", "ps"), # use B-splines for both
                     mc = c(TRUE, TRUE), # sum-to-zero constraints for both x and t
                     k = c(8, 8)), # 8 basis functions for both x and t (before constraint)
                data = dat_lon) 

plot(gam_mod2)
# random intercept is closer to FDboost
# x1 is fuzzy
# check constraint
dat_lon$yclr_pred <- predict(gam_mod2)
yclrpred_gam <- array(dat_lon$yclr_pred, dim = dim(yclr))
stopifnot(max(abs(colSums(yclrpred_gam))) < 1e-12)

## re-gain and plot predicted densities
ypred_gam <- sweep(exp(yclrpred_gam), 2, colMeans(exp(yclrpred_gam)), "/")
matplot(ypred_gam, t = "l")

# plus covariates
gam_mod3 <- gam(yclr ~ 0 +
                  ti(t, bs = "ps", mc = TRUE, k = 8) + # functional intercept
                  ti(group, bs = "re", mc = TRUE) + # random effect
                  s(t, by = x2, bs = "ps", k = 5, m = c(2,1)) + # factor covariate
                  ti(x1, t, # covariate and time
                     bs = c("ps", "ps"), # use B-splines for both
                     mc = c(TRUE, TRUE), # sum-to-zero constraints for both x and t
                     k = c(8, 8)) + # 8 basis functions for both x and t (before constraint)
                  ti(x3, t, # covariate and time
                     bs = c("ps", "ps"), # use B-splines for both
                     mc = c(TRUE, TRUE), # sum-to-zero constraints for both x and t
                     k = c(8, 8)),            
                  data = dat_lon) 

plot(gam_mod3)

dat_lon$yclr_pred <- predict(gam_mod3)
yclrpred_gam <- array(dat_lon$yclr_pred, dim = dim(yclr))
stopifnot(max(abs(colSums(yclrpred_gam))) < 1e-12)

## re-gain and plot predicted densities
ypred_gam <- sweep(exp(yclrpred_gam), 2, colMeans(exp(yclrpred_gam)), "/")
matplot(ypred_gam, t = "l")


## compare true curves, data curves, and predicted curves from flexible FDboost and flexible gam model
ylim <- range(ypred1, y, ypred_gam)

opar <- par(mfrow = c(1,3))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), ylim = ylim, type = "l", main = "true underlying beta-densities")
matplot(t, y, type = "l", ylim = ylim, main = "data consisting of kernel-density-estimates")
matplot(t, ypred_gam, type = "l", ylim = ylim, main = "predicted densities in flexible gam(/refund) model")
par(opar)


# Tensor Smooths ----------------------------------------------------------


