# Haibo Jin
# haibo.nick.jin@gmail.com
# 014343698
##################################################################
# generate data
alpha <- 7; beta <- 9; nY <- 50
set.seed( 100 )
y <- rgamma( nY, alpha, beta)

# unnormalized posterior distribution
target <- function( a, b, lambda, y) {
  - lambda * ( a + b ) + length(y) * a * log( b ) - 
  length( y ) * lgamma( a ) + ( a - 1 ) * sum( log( y ) ) - 
  b * sum( y )  
}

# use grid method to find approximate value of alpha and beta
grid <- seq(0.1, 10, length.out = 100)
x <- matrix(0, 100,100)
for ( ii in seq(1,100)) {
  for (jj in seq(1,100)) {
    x[ii, jj] <- exp( target(grid[ii], grid[jj], 0.001, y) )
  }
}
index <- which(x == max(x), arr.ind = TRUE)
alpha_est <- grid[index[1]]
beta_est <- grid[index[2]]

# finally we get the approximate mode of the posterior 
# alpha = 6.9, beta = 9.4

# calculate hessian at the mode for normal approximation
hessian <- matrix( c( - length( y ) * trigamma( alpha_est ),
                      length( y ) / beta_est,
                      length( y ) / beta_est,
                      - length( y ) * alpha_est / ( beta_est ^ 2 ) ), 
                   2, 2)
# the covariance of normal approximation is the inverse of the negative hessian
sigma_est <- solve( - hessian )
# so here is the parameter for normal approximation
# mu = (6.9, 9.4) sigma = [1.82, 2.48, 2.48, 3.63]

# sampling loops
a <- 5; Sigma <- a * sigma_est;
nSamples <- 10000; lambda <- 0.001; nAccepted <- 0;
x <- matrix( 0, nSamples, 2 ); x[1, ] <- c(3,3);
for ( ii in seq( 2, nSamples ) ) {
  x[ ii, ] <- mvtnorm::rmvnorm( 1, x[ ii - 1, ], Sigma)
  if ( x[ii, 1] <= 0 || x[ii, 2] <= 0 ) {
    x[ii, ] <- x[ii - 1, ]
    next
  }
  r <- exp( target( x[ii, 1], x[ii, 2], lambda, y) - 
                target( x[ii - 1, 1], x[ii - 1, 2], lambda, y) )
  if ( runif( 1 ) > r ) {
    x[ii, ] <- x[ii - 1, ]
  } else {
    nAccepted <- nAccepted + 1
  }
}
# acceptance rate
nAccepted / nSamples
# trace plot
plot( seq(1, nSamples), x[ , 1], "l", main = "trace plot of alpha")
plot( seq(1, nSamples), x[ , 2], "l", main = "trace plot of beta" )
# summary statistics
summary(x)
# standard error of alpha
sd( x[,1] ) / sqrt( nSamples )
# standard error of beta
sd( x[,2] ) / sqrt( nSamples )
# autocorrelation
acf(x[,1], main = "autocorrelation of alpha")
acf(x[,2], main = "autocorrelation of beta")
# cumulative average plot of alpha
cumsum_a <- cumsum(x[,1])
cumave_a <- cumsum_a / seq(1, nSamples)
plot(seq(1, nSamples), cumave_a, "l", main = "cumulative average of alpha")
# cumulative average plot of beta
cumsum_b <- cumsum(x[,2])
cumave_b <- cumsum_b / seq(1, nSamples)
plot(seq(1, nSamples), cumave_b, "l", main = "cumulative average of beta")
