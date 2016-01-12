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

# proposal distribution
proposal <- function( a, b, lambda, y) {
  log( dgamma(b, length( y ) * a + 1, lambda + sum( y ) ) )
}

# sampling loops
a <- 1; Sigma <- a * 1;
nSamples <- 10000; lambda <- 0.001; nAccepted <- 0;
x <- matrix( 0, nSamples, 2 ); x[1, ] <- c(1,1);
for ( ii in seq( 2, nSamples ) ) {
  x[ii, 1] <- rnorm( 1, x[ ii - 1, 1], Sigma) 
  x[ii, 2] <- rgamma( 1, length(y)*exp(x[ii,1])+1, lambda+sum(y))
  r <- exp( target( exp(x[ii, 1]), x[ii, 2], lambda, y) - 
              target( exp(x[ii - 1, 1]), x[ii - 1, 2], lambda, y) +
              proposal( exp(x[ii - 1, 1]), x[ii - 1, 2], lambda, y ) - 
              proposal( exp(x[ii, 1]), x[ii, 2], lambda, y ) )
  if ( runif( 1 ) > r ) {
    x[ii, ] <- x[ii - 1, ]
  } else {
    nAccepted <- nAccepted + 1
  }
}
# transform phi to alpha
x[,1] <- exp(x[,1])
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
