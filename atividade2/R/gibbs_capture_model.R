### Codes for Exercise 2 - MI602-2s2020
### Gibbs Sampling applied to capture-recapture models
### Author: josehcms
### Last update: 2020-11-18

### Housekeeping #----------------------
rm( list = ls( ) )
graphics.off( )
require( dplyr ); require( data.table ); 
require( ggplo2 )
#########################################

### Set initial variables #--------------

# set data matrix with capture-recapture data
data_capture <- 
  matrix(
    data = c( 2, 4, 2, 4, 4, 2, 3, 4, 2, 3, 5, 10, 2, 5,
              0, 0, 2, 4, 4, 2, 3, 2, 2, 3, 5, 6, 2, 5 ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(  c( 'ni', 'mi' ), c( 1:14 ) ) )

# set values for lambda, r, sample size n_samp and I 
r <- ( data_capture[ 'ni',] - data_capture[ 'mi',] ) %>% sum
lambda <- 50
n_samp <- 100000
I <- 14
##########################################

### Run Gibbs sampling #------------------
set.seed( 123456 )

# set N vector 
f_N <- numeric( n_samp )
# set matrix of pis
f_p <- matrix( nrow = n_samp, ncol = I )

# set priors for N ~ Pois( lambda = 50 )
f_N[ 1 ] <- 
  rpois( 1, lambda )

# set priors for pi ~ uniform( 0, 1 )
f_p[ 1, ] <- 
  runif( 14, min = 0, max = 1 )

# run sample size n_samp
for( n in 2:n_samp ){

  # generate f(Ni) from f(pi-1)
  f_N[ n ] <- 
    rpois( 1, 
           lambda = ( lambda * prod( ( 1 - f_p[ n - 1, ] ) ) ) ) + r
  
  # generate f(pi) from f(Ni)
  for( j in 1:14 ){
    f_p[ n, j ] <-
      rbeta( 1,
             shape1 = ( data_capture[ 'ni', j ] + 1 ), 
             shape2 = ( f_N[ n ] - data_capture[ 'ni', j ] + 1 ) )  
    
  }
    
}

# results in data.table
res_gibbs <- 
  cbind( f_N, f_p, 1:100000 ) %>%
  as.data.table %>%
  setnames( c( 'N', paste0( 'p',1:14 ), 't' ) ) %>%
  melt(
    id.vars       = c( 't' ),
    measure.vars  = c( 'N', paste0( 'p',1:14 ) ),
    value.name    = 'est',
    variable.name = 'type'
  )
##########################################

### Verify convergence and autocorrelation #------

# plot 200 first samples to verify convergence
res_gibbs[ t %in% 0:200 ] %>%
  ggplot() +
  geom_line( aes( x = t, 
                  y = est ),
             size = 0.75 ) +
  labs( x = 'Tempo (t)',
        y = 'Valores gerados\n(Sequência de Gibbs)' ) +
  facet_wrap( ~ type, scales = 'free_y', nrow = 3 ) +
  theme_light() +
  theme( 
    strip.text = element_text( color = 'black', size = 12 ),
    axis.title = element_text( color = 'black', size = 13 ),
    axis.text  = element_text( color = 'black', size = 11 )   
  )

# autocorrelation data by variable
acf_dat <- data.table()
for( i in unique( res_gibbs$type ) ){
  aux <- 
    res_gibbs[ type == i ] 
  
  acf_dat <- 
    rbind(
      acf_dat,
      data.table(
        type   = i,
        acf    = acf( aux$est, plot = F )$acf,
        lag    = acf( aux$est, plot = F )$lag
      )
    )
}

# 95% CI
ciline = qnorm( (1 - 0.95 ) / 2 ) / sqrt(  length( f_N ) )

# plot acf results
acf_dat[ lag <= 20 & lag > 0 ] %>%
  ggplot() +
  geom_segment( aes( x = lag, xend = lag, 
                     y = 0,   yend = acf ),
             size = 0.75 ) +
  geom_point( aes( x = lag,
                   y = acf ),
              size = 1 ) +
  geom_hline( yintercept = ciline,
              linetype = 'longdash',
              size = 0.25,
              color = 'red' ) +
  geom_hline( yintercept = -ciline,
              linetype = 'longdash',
              size = 0.25,
              color = 'red' ) +
  scale_x_continuous( limits = c( 0, 20 ) ) +
  labs( x = 'Lag',
        y = 'Autocorrelação' ) +
  facet_wrap( ~ type, nrow = 3 ) +
  theme_light() +
  theme( 
    strip.text = element_text( color = 'black', size = 12 ),
    axis.title = element_text( color = 'black', size = 13 ),
    axis.text  = element_text( color = 'black', size = 11 )   
  )


##########################################


### Compute descriptive statistics of fidings #------
res_gibbs_final <- 
  res_gibbs[ t %in% seq( 100, n_samp, 10 ) ]

res_gibbs_final[ ,
                 .( t025 = round( quantile( est, probs = 0.025 ), 3 ),
                    t975 = round( quantile( est, probs = 0.975 ), 3 ),
                    mean = round( mean( est ), 3 ),
                    med  = round( median( est ), 3 ) ),
                 .( type ) ]

##########################################

### The End