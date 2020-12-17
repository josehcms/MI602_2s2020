
rm( list = ls( ) )

require( dplyr ); require( data.table ); require( sf )
require( ggplot2 ); require( rgdal ); require( rgeos )
require( matlib )

### Read adjacencies matrix #--------------------
rj_adj <- 
  fread( 'atividade3/data/rio_adj_at3.csv' ) %>%
  .[ , - c( 'Bangu', 'Gericinó', 'Jabour', 'Vila Kennedy' ) ] %>%
  .[ ! ( V1 %in% c( 'Bangu', 'Gericinó', 'Jabour', 'Vila Kennedy' ) ) ]

# create labels for names
labels <-  paste0( 'B', 1 : nrow( rj_adj ) )
names( labels ) <- rj_adj$V1

A <- 
  rj_adj[,-1] %>%
  as.matrix

row.names( A ) <- labels
colnames( A ) <- labels
###################################

### Read dependent variabel #------
Y <- 
  fread( 'atividade3/data/rio_atv3.csv' ) %>%
  .[ , cod_neighb := factor( labels[ as.character( neighborhood ) ],
                             levels = labels ) ] %>%
  setorder( cod_neighb ) %>%
  .[ , tx := 100000 * count / Inhabitants ] %>%
  pull( count )
  
grad <- function( rho, mu, v, A, Y ){
    
    n = nrow( A )
    I = diag( n )
    one <- matrix( 1, ncol = 1, nrow = n ) 
    
    dldv <-  
      -n / v + 1 / (v^2) * t( ( Y - mu * 1 ) ) %*% 
      matrix.power( ( I - rho * A ), 2 )  %*% ( Y - mu * 1 )  
    
    dldmu  <-  
      1 / ( 2 * v ) * t( Y ) %*% ( I - rho * A  ) %*% one +
      1 / ( 2 * v ) * t( one ) %*% matrix.power( ( I - rho * A ), 2 ) %*% Y +
      mu / ( v ^ 2 ) * t( one ) %*% ( I - rho * A ) %*% one
    
    dldrho <- 
      tr( ( A %*% solve(I-rho*A,I) ) )  + 
      ( t( Y - mu * one ) %*% t( I - rho * A ) %*% A %*% ( Y - mu * one ) ) / v
    
    return( matrix( c( dldmu, dldv, dldrho ),
                    ncol = 1,
                    byrow = 1 ) )
  }
  
hess <- function( rho, mu, v, A, Y ){
    
    n = nrow( A )
    I = diag( n )
    one <- matrix( 1, ncol = 1, nrow = n ) 
    
    dl2dmu2 <- 
      - ( 1 / v ) * t( one ) %*% matrix.power( ( I - rho * A ), 2 ) %*% one
    
    dl2dmudv <- 
      - 1 / ( 2 * ( v ^ 2 ) ) * t( Y ) %*% matrix.power( ( I - rho * A ), 2 ) %*% one -
      1 / ( 2 * ( v ^ 2 ) ) * t( one ) %*% matrix.power( ( I - rho * A ), 2 ) %*% Y +
      mu / ( v ^ 2 ) * t( one ) %*% matrix.power( ( I - rho * A ), 2 ) %*% one 
  
    dl2dmudrho <- 
      - 1 / ( v ) * 
      ( 
        t( Y ) %*% t( I - rho * A ) %*% A %*% one - 
        t( one ) %*% t( I - rho * A ) %*% A %*% Y +
          2 * mu * t( one ) %*% t( I - rho * A ) %*% A %*% one
      )
      
    dl2dv2 <- 
      n / ( 2 * ( v ^ 2 ) ) - 
      1 / ( v ^ 3 ) * t ( Y - mu * one ) %*% matrix.power( ( I - rho * A ), 2 ) %*% 
      ( Y - mu * one )
    
    dl2dvdmu = dl2dmudv
    
    dl2dvdrho <- 
      - 1 / ( v ^ 2 ) * t ( Y - mu * one ) %*% t( I - rho * A ) %*% 
      A %*% ( Y - mu * one )
    
    dl2drho2 <- 
      - tr( A %*% solve(I-rho*A,I) %*% A %*% solve(I-rho*A,I) ) +
      ( t( Y - mu * one ) %*% t( A ) %*% A %*% ( Y - mu * one  ) ) / v
    
    dl2drhodv <- dl2dvdrho
    
    dl2drhodmu <- dl2dmudrho
    
    hess_out <- 
      matrix( c( dl2dmu2, dl2dmudv, dl2dmudrho,
                 dl2dvdmu, dl2dv2, dl2dvdrho,
                 dl2drhodmu, dl2drhodv, dl2drho2 ),
              ncol = 3,
              byrow = T )
  
    return( hess_out )
  }
 
x1 <- 
  c( pars[1], pars[2]^2, pars[3] )

for( i in 1:4 ){
  x1 = 
    x1 - 
    solve( hess( rho = x1[3], mu = x1[1], v = x1[2], A, Y  ),
           grad( rho = x1[3], mu = x1[1], v = x1[2], A, Y  ) )
    
}


mloglik <- 
  function( pars, Y, A ){
    require( expm )
    require( matrixcalc )
    mu = pars[1]
    sdev = pars[2]
    rho = pars[3]
    n = nrow( A )
    I = diag( n )
    one <- matrix( 1, ncol = 1, nrow = n ) 
    
    # Likfunc <-
    #   (2*pi)^(-3/2) * sdev^(-3)*det(I-rho*A)*
    #   exp((-1/(2*(sdev^2) ))*t(Y-mu*one)%*%matrix.power(I-rho*A,2)%*%(Y-mu*one))
    # 
    # logLik <- log( Likfunc )
    # 
    logLik <-
      - ( n / 2 ) * log( 2 * pi ) -
     n * log( sdev ) +
      log( det( I - rho * A ) ) -
      ( 1 / ( 2 * (sdev^2) ) ) * 
      t( Y - mu * one ) %*% ( matrix.power( ( I - rho * A ), 2 ) ) %*%
      ( Y - mu * one )

    return( - logLik )
  }

pars = c( mean(Y), sd(Y), 0.72 )


mle <- 
  nlm( mloglik,
       pars,
       hessian = TRUE,
       Y = Y,
       A = A
  )

mle$estimate[1]
mle$estimate[2] 
mle$estimate[3]

require(MASS)
require(matrixcalc)
install.packages('mvtnorm')
require(mvtnorm)

dmvnorm( mean = ( mu * one ), 
         sigma = ( v * matrix.power( ( I - rho * A ), -2 ) ) )
