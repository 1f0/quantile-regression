#-----data generation

#  super parameters
nr <- 200 # train
ns <- 50  # test
n <- nr + ns
set.seed(496)

method <- "LCQR-MDM"

#  samples
X <- runif(n)


#  error config:
#  normal error
#  - N(0, 1)
#  - t(3)-distribution
#  - 0.9N(0,1)+0.1N(0,10^2)
#  - Cauchy(0, 1)
#  non-sym error
#  - Exp(1)
#  - LN(0, 1)
#  - chi^2(1)
#  - Pa(1/3, 1, 0)
#  - Pa(1/2, 1, 0)
#  - Pa(1, 1, 0)
eps <- rnorm(n)


generate_data <- function(X, eps) {
  return(X * sin(2 * pi * X) + (2 + cos(2 * pi * X)) / 20 * eps)
}


Y <- generate_data(X, eps)
dat <- data.frame(X, Y)

# Linear model
if (method == "OLS") {
  mod <- lm(Y ~ X)
  print(summary(mod))
  plot(Y ~ X)
  abline(mod, col="red")
}

if (method == "LL-IS") {
  alpha <- 0.2 # proportion of neighbor

  mod <- loess(Y ~ X, degree = 1, span = alpha)
  Z <- predict(mod)
  print(mod)

  plot(Y ~ X)

  ## plot method 1, subtlely different from method 2
  # j <- order(X)
  # lines(X[j], Z[j], col="red")

  # plot method 2, since scatter.smooth cannot set line color
  lines(loess.smooth(X, Y, degree = 1, span = alpha), col = "red")

  #method 3: not work yet
  # lattice::xyplot(Y ~ X, data=dat, type=c('p'), col.line='red')
}


# LCQR-MDM
bandwidth <- 0.2
r <- 5 # or 9, 19 # should be odd num so that median exists

if (method == "LCQR-MDM") {
  # quantile position
  tau_r <- seq(from = 1/(r+1), to = 1 - 1/(r+1), length.out = r)
  
  #LQR
  # mod <- McSpatial::qreglwr(Y~X, taumat=tau_r, kern="epan",target="alldata")
  # j <- order(X)
  # plot(Y~X)
  # for (i in 1:r) {
  #   lines(X[j], mod$yhat[j,i], col="red")
  # }

  library("quantreg")
  if (!exists("cqreglwr", mode="function")) source("cqreglwr.R")
  mod <- cqreglwr(Y~X, taumat=tau_r, kern="epan",target="alldata")
  print(mod)
}


# LCER-MDM





