# Haibo Jin
# haibo.nick.jin@gmail.com
# 014343698
##################################################################
# generate data
alpha <- 7; beta <- 9; nY <- 50
set.seed( 100 )
y <- rgamma( nY, alpha, beta)

# unnormalized posterior distribution
target <- function( phi, psi, lambda, y) {
  - lambda * ( exp( phi ) + exp( phi - psi ) ) +
    length(y) * exp( phi ) * ( phi - psi ) - 
    length( y ) * lgamma( exp( phi ) ) + 
    ( exp( phi ) - 1 ) * sum( log( y ) ) - 
    exp( phi - psi) * sum( y ) + 2 * phi - psi  
}

# sampling loops
a <- 0.1; Sigma <- a * matrix( c( 1, 0, 0, 1), 2, 2 );
nSamples <- 10000; lambda <- 0.001; nAccepted <- 0;
x <- matrix( 0, nSamples, 2 ); x[1, ] <- c(1,1);
for ( ii in seq( 2, nSamples ) ) {
  x[ ii, ] <- mvtnorm::rmvnorm( 1, x[ ii - 1, ], Sigma)
  r <- exp( target( x[ii, 1], x[ii, 2], lambda, y) - 
              target( x[ii - 1, 1], x[ii - 1, 2], lambda, y) )
  if ( runif( 1 ) > r ) {
    x[ii, ] <- x[ii - 1, ]
  } else {
    nAccepted <- nAccepted + 1
  }
}
# transform phi and psi to alpha and beta
xx <- matrix( 0, nSamples, 2 )
xx[, 1] <- exp( x[, 1] ); xx[, 2] <- exp( x[, 1] - x[, 2] );
# acceptance rate
nAccepted / nSamples
# trace plot
plot( seq(1, nSamples), xx[ , 1], "l", main = "trace plot of alpha")
plot( seq(1, nSamples), xx[ , 2], "l", main = "trace plot of beta" )
# summary statistics
summary(xx)
# standard error of alpha
sd( xx[,1] ) / sqrt( nSamples )
# standard error of beta
sd( xx[,2] ) / sqrt( nSamples )
# autocorrelation
acf(xx[,1], main = "autocorrelation of alpha")
acf(xx[,2], main = "autocorrelation of beta")
# cumulative average plot of alpha
cumsum_a <- cumsum(xx[,1])
cumave_a <- cumsum_a / seq(1, nSamples)
plot(seq(1, nSamples), cumave_a, "l", main = "cumulative average of alpha")
# cumulative average plot of beta
cumsum_b <- cumsum(xx[,2])
cumave_b <- cumsum_b / seq(1, nSamples)
plot(seq(1, nSamples), cumave_b, "l", main = "cumulative average of beta")
