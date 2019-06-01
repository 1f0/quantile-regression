# data generation

#  super parameters
n <- 200
r <- 5 # 9, 19

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
eps <- rnorm(0, 1)


generate_data <- function(X, eps) {
  return(X * sin(2 * pi * X) + (2 + cos(2 * pi * X)) / 20 * eps)
}

Y = generate_data(X, eps)
hist(Y)


# LCQR-MDM


# LCER-MDM


# LL-IS


