#-----data generation

#  super parameters
nr <- 200 # train
ns <- 50  # test
n <- nr + ns
r <- 5 # or 9, 19

method = "LR"
#  location
X <- runif(n)
# kern == max(0, epan)
kern <- function(psi) {
  q <- 1 - psi^2
  q[q < 0] <- 0
  return(0.75 * q)
}


bandwidth <- 0

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

Y = generate_data(X, eps)

# LL-IS

# Linear Regression
if (method == "LR") {
  dat = data.frame(X, Y)
  mod = lm(Y ~ X, data = dat)
  print(summary(mod))
  plot(Y ~ X, data = dat, cex = 0.5)
  abline(mod, col="red")
}

if (method == "LL-IS") {

}

# LCQR-MDM



# LCER-MDM





