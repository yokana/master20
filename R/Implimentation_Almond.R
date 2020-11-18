### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
##       Sketch of regression with a density as response    ##
##               Author: Almond Stoecker                    ##
##                      6.12.2017                           ##
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

## Note:
# This is only a brief sketch of possible implementations. 
# Please, read/use carefully as parts may be not optimal.
# e.g., there are no functional random intercepts applied to account for in-curve dependencies
# and the exact constrains for the effects are worth a closer look.

#### Simulate data ####

## bayes-space linear model with beta-distribution-intercept
# reparametrized distribution function for beta distribution such that shape1-1 = a and shape2-1 = b
dbeta2 <- function(x, a, b, ncp = 0, log = FALSE) dbeta(x, a+1, b+1, ncp, log)
rbeta2 <- function(x, a, b, ncp = 0) rbeta(x, a+1, b+1, ncp)
# intercept parameters
a0 <- 3
b0 <- 3
# time grid
t <- seq(0,1, length.out = 50)
# plot intercept distribution
plot(t, dbeta2(t, a0, b0), type = "l")
# plot(t, dbeta(t, 4, 4), type = "l") 

## specify slope
slope <- 1

# sample size
N <- 50

## sample covariates >-1/(slope*min(a0,b0)) from exponential distribution substracted by the same constant
set.seed(493822)
x <- rexp(50)-1/(slope*min(a0, b0))
range(x)

## determine vectors of true underlying parameters
a <- slope*x*a0
b <- slope*x*b0
# Slope only changes the "thickness* of the density; 
# COvariate X conributes to the density as well: higher X values -> higher density in the center
# We can actually see the three distributions with very high x values


# plot true distribution functions
# use matplot to plot matrices
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l")
paraMatrix <- cbind(a,b)
plot(dbeta(t, paraMatrix[1,1], paraMatrix[1,2]), type = "l")

## sample response curves in a realistic procedure by sampling from the corresponding distributions
# number of samples per density function
n <- 1000
nsamples <- apply(cbind(a,b), 1, function(p) rbeta2(n, p[1], p[2]))
# use structure to put density objects in a list
density_esimates <- structure(apply(nsamples, 2, density, bw = .05, n = length(t), from = 0, to = 1))
# subset each density object 
# each column is a density; row represents evaluation points
y <- sapply(density_esimates, function(x) x$y)
# How far are those from the true functions
par(mfrow = c(1,2))
plot(dbeta(t, paraMatrix[1,1], paraMatrix[1,2]), type = "l")
plot(density_esimates[[1]])
# there is a bit of noise, but not much => try with smaller n

# for all
par(mfrow = c(1,2))
matplot(t, y, type = "l", main = "Sample", ylim = c(0,5))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l", main = "Truth", ylim = c(0,5))
# even though we sample 1000 estimates from each function we still get a lot of noise
# what if we only observe around 20 points?


## apply the clr-transform of Talska et al. 2017 to get standard L^2 functions
yclr <- base::scale( log(y), center = TRUE, scale = FALSE )
# is this equivalent to clr()?
# clr takes a 
test1 <- as.matrix(t(y))
test <- clr(test1)
# almost identical
sum(yclr[,1])
sum(test[1,])
# both some up to zero though

# not at all; it seems clr estimation (geometric mean instead of mean returned by scale function) doesnt spread as much
plot(yclr[,1], type = "l", col = "blue", ylim = c(-5,2))
lines(unclass(test[1,]), type = "l", col = "red")
# identical -> clr produces is class type object, which might be more suitable for Boogharts approach

stopifnot(all(colMeans(yclr)<1e-10))
matplot(t, yclr, type = "l")

## Prepare data for FDboost and refund ## 
dat <- list(yclr = t(yclr), t = t, x = x)
# data structure is the same for both packages

#### Fit regression model with FDboost ####

library(FDboost)

## fit linear model (as original model is linear in Bayes-space) ##
# remember x is continuous -> use bolsc which is equivalent to zero-sum constraint => Do we need this for categorical variables??
model1 <- FDboost(yclr ~ 1 + bolsc(x, df = 2), timeformula = ~bbsc(t, df = 2), data = dat) # the base-learners bolsc and bbsc correspond to a linear effect bols and a smooth spline effect bbs, but implement a sum-to-zero constraint over the values of x and t, respectively. See ?bolsc, ?bbsc.
plot(model1, which = 2)

# Smooth Intercept Model does resemble clr()-beta(4,4) -> transform backwards to check
# smooth effect of x corresponds to our expectation

# predict respondent densities
# predict returns a matrix with each row representing a density
yclrpred1 <- t(predict(model1))
clrTest <- predict(model1)
plot(clrTest[1,], type = "l")
plot(yclr[,1], type = "l")
# check for sum-to-zero constraint
stopifnot(all(colMeans(yclrpred1)<1e-10))
# apply the reverse clr-transformation to obtain density data from each column
# divide exp(y_i) by geometric mean of exp(y_i)
ypred1 <- sweep(exp(yclrpred1), 2, colMeans(exp(yclrpred1)), "/")
test <- clrInv(clrTest)
opar <- par(mfrow = c(1,2))
matplot(t,ypred1[,1],type = "l")
matplot(t,test[1,],type = "l")
# the clr inverse doesnt look like a density even though its shape is identical
par(opar)


# check iff column mean is 1
stopifnot(all(abs(colMeans(ypred1)-1)<1e-10))

## compare true curves, data curves, and predicted curves
opar <- par(mfrow = c(1,3))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l", main = "true underlying beta-densities")
matplot(t, ypred1, type = "l", main = "predicted densities in linear model")
matplot(t, y, type = "l", main = "data consisting of kernel-density-estimates")
par(opar)

## compare estimated intercept with true underlying intercept 
# (or rather the curve with the parameters specified as the intercept above)
# produce new data with one observation at x = 1/slope corresponding to a0 and b0
newdat <- dat
# replace first value of X with the inverse of the slope to cancel out the effect of the slope and get b0 and a0
newdat$x[1] <- 1/slope
# predict the first distribution given that x equals 1/slope
a0b0estimateclr1 <- predict(model1, newdata = newdat)[1,]
# backtransform
a0b0estimate1 <- exp(a0b0estimateclr1)/mean(exp(a0b0estimateclr1))

plot(t, dbeta2(t, a0, b0), t = "l", lty = 1, main = "True basis density of beta(a0,b0)-distribution \n and predicted density by linear model")
lines(t, a0b0estimate1, col = "cornflowerblue", lty = 2)
legend("bottomright", legend = c("true", "predicted"), col = c("black", "cornflowerblue"), lty = 1:2)


plot(model1, which = 1) # is that the same
plot(t, a0b0estimateclr1, t = "l", lty = 1, main = "True basis density of beta(a0,b0)-distribution \n and predicted density by linear model")
# it is corresponding to the smooth Intercept model

## fit flexible model ##
# that is with a smooth effect for continous x
model2 <- FDboost(yclr ~ 1+ bbsc(x, df = 2), timeformula = ~bbsc(t, df = 2), data = dat) # set mstop to 500 for letting the fitting algo run longer
plot(model2, which = 2)
# offset remains the same
# effects are more distinct but overall ppicture remains the same


yclrpred2 <- t(predict(model2))
stopifnot(all(colMeans(yclrpred2)<1e-10))
# apply the reverse clr-transformation to obtain density data
ypred2 <- sweep(exp(yclrpred2), 2, colMeans(exp(yclrpred2)), "/")
stopifnot(all(abs(colMeans(ypred2)-1)<1e-10))

## compare true curves, data curves, and predicted curves
opar <- par(mfrow = c(2,2))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), type = "l", main = "true underlying beta-densities")
matplot(t, ypred1, type = "l", main = "predicted densities in linear model")
matplot(t, y, type = "l", main = "data consisting of kernel-density-estimates")
matplot(t, ypred2, type = "l", main = "predicted densities in flexible model")
par(opar)

## compare estimated intercept with true undelying intercept 
# (or rather the curve with the parameters specified as the intercept above)
# produce new data with one observation at x = 1/slope corresponding to a0 and b0
newdat <- dat
newdat$x[1] <- 1/slope
a0b0estimateclr2 <- predict(model2, newdata = newdat)[1,]
a0b0estimate2 <- exp(a0b0estimateclr2)/mean(exp(a0b0estimateclr2))

plot(t, dbeta2(t, a0, b0), t = "l", lty = 1, main = "True basis density of beta(a0,b0)-distribution \n and predicted densities")
lines(t, a0b0estimate1, col = "cornflowerblue", lty = 2)
lines(t, a0b0estimate2, col = "darkseagreen", lty = 3)
legend("bottomright", legend = c("true", "linear model", "flexible model"), col = c("black", "cornflowerblue", "darkseagreen"), lty = 1:3)


#### Fit regression model with refund ####

library(refund)
## Not clear how to implement the sum-to-zero constraint also for the development over time in pffr
model3 <- pffr( yclr ~ 1 + s(x, bs = "ps"), yind = t, data = dat, bs.yindex = list(mc = c(TRUE,TRUE)))
# i.e. we would need mc = c(TRUE,TRUE) in the mgcv formula:
model3$formula
## -> so manually specify the function in mgcv::gam (pffr is a wrapper for gam)
# 
coef(model3)
library(mgcv)
## built the data for gam (just a loong data.frame with y in one column as for scalar data)
datgam <- data.frame(yclr = as.vector(yclr), 
                     t = rep(t,N), x = rep(x, each = length(t)), 
                     id = rep(factor(1:N), each = length(t)))
model4 <- gam(yclr ~ 0 +
                ti(t, bs = "ps", mc = TRUE, k = 8) + # functional intercept
                ti(x, t, # covariate and time
                   bs = c("ps", "ps"), # use B-splines for both
                   mc = c(TRUE, TRUE), # sum-to-zero constraints for both x and t
                   k = c(8, 8)), # 8 basis functions for both x and t (before constraint)
              data = datgam) 

## check zero summation
datgam$yclr_pred <- predict(model4)
yclrpred_gam <- array(datgam$yclr_pred, dim = dim(yclr))
stopifnot(max(abs(colSums(yclrpred_gam))) < 1e-12)

## re-gain and plot predicted densities
ypred_gam <- sweep(exp(yclrpred_gam), 2, colMeans(exp(yclrpred_gam)), "/")
matplot(ypred_gam, t = "l")

## compare true curves, data curves, and predicted curves from flexible FDboost and flexible gam model
ylim <- range(ypred1, y, ypred_gam)

opar <- par(mfrow = c(2,2))
matplot(t, apply(cbind(a,b), 1, function(p) dbeta2(t, p[1], p[2])), ylim = ylim, type = "l", main = "true underlying beta-densities")
matplot(t, ypred1, type = "l", ylim = ylim, main = "predicted densities in flexible FDboost model")
matplot(t, y, type = "l", ylim = ylim, main = "data consisting of kernel-density-estimates")
matplot(t, ypred_gam, type = "l", ylim = ylim, main = "predicted densities in flexible gam(/refund) model")
par(opar)

# Comparing the two different approaches, the FDboost regression model seems to perform slightly worse at the borders of the domain.
# However, this might be due to the unspecified stopping iteration in the boosting algorithm defaulting to mstop = 100. 
# This is of special relevance, as in the present case, 
# where a relatively simple model is fitted and perfectly specified, the optimal mstop is expected to be much higher. 