
rm( list = ls( ) )
graphics.off( )

require( dplyr ); require( data.table )

n_samp <- 100000
I <- 14

data_capture <- 
  matrix(
    data = c( 2, 4, 2, 4, 4, 2, 3, 4, 2, 3, 5, 10, 2, 5,
              0, 0, 2, 4, 4, 2, 3, 2, 2, 3, 5, 6, 2, 5 ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(  c( 'ni', 'mi' ), c( 1:14 ) ) )

r <- (data_capture[ 'ni',] - data_capture[ 'mi',]) %>% sum
lambda <- 50
set.seed( 123456 )

# set N vector and N(0)
f_N <- numeric( n_samp )
# set matrix of pis and pi(0)
f_p <- matrix( nrow = n_samp, ncol = I )

f_N[ 1 ] <- 
  rpois( 1, lambda )

f_p[ 1, ] <- runif( 14,min = 0, max = 1 )

for( n in 2:n_samp ){

  f_N[ n ] <- 
    rpois( 1, 
           lambda = ( lambda * prod( ( 1 - f_p[ n - 1, ] ) ) ) ) + r
  
  for( j in 1:14 ){
    f_p[ n, j ] <-
      rbeta( 1,
             shape1 = ( data_capture[ 'ni', j ] + 1 ), 
             shape2 = ( f_N[ n ] - data_capture[ 'ni', j ] + 1 ) )  
    
  }
    
}

quantile(f_N, probs = seq(0,1,0.1))

plot( density(f_N) , type = 'l')
lines( density(rpois(10000,lambda)), col = 'red')
barplot()
line( )