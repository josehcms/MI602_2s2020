### Codes for Exercise 3 - MI602-2s2020
### Newton-Raphson optimization for SAR model 
### parameter estimation
### Author: josehcms
### Last update: 2020-12-16

### Housekeeping #----------------------
rm( list = ls( ) )
graphics.off( )

req_packs <- 
  c( 'dplyr', 'data.table', 'sf',
     'ggplot2', 'rgdal', 'rgeos',
     'matlib', 'matrixcalc' )

lapply( req_packs, require, character.only = TRUE )

#########################################

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
atv3_dt <- 
  fread( 'atividade3/data/rio_atv3.csv' ) %>%
  .[ , cod_neighb := factor( labels[ as.character( neighborhood ) ],
                             levels = labels ) ] %>%
  setorder( cod_neighb ) %>%
  .[ , tx := 100000 * count / Inhabitants ]  %>%
  .[ , countClass := cut( count, 
                          breaks = c( 0, 5, 20, 50, Inf ),
                          labels = c( '0-4', '5-19', '20-49', '50+' ),
                          right  = FALSE ) ] 

Y <- 
  atv3_dt %>%
  copy %>%
  pull( count )
###################################

### Log Likelihood Minimization Function #--------

mloglik <- 
  function( pars, Y, A ){
    
    mu = pars[1]
    sdev = pars[2]
    rho = pars[3]
    n = nrow( A )
    I = diag( n )
    one <- matrix( 1, ncol = 1, nrow = n ) 
    
    logLik <-
      - ( n / 2 ) * log( 2 * pi ) -
      n * log( sdev ) +
      log( det( I - rho * A ) ) -
      ( 1 / ( 2 * (sdev^2) ) ) *
      t( Y - mu * one ) %*% ( matrix.power( ( I - rho * A ), 2 ) ) %*%
      ( Y - mu * one )
    
    return( - logLik )
  }

###################################

### Define list of different initial values #------

# First list: fix mu and sigma2 in mean(Y) and sd(Y)
pars_list1 <- 
  lapply( seq( 0.01, 0.95, 0.01 ),
          function( x ){
            c( mean(Y), sd(Y), x )
          } )

saved_estpars1 <- data.table()

for( i in 1:length( pars_list1 ) ){
  
  aux_pars <- pars_list1[[i]]
  
  NRest <- 
    nlm( mloglik,
         aux_pars,
         hessian = TRUE,
         Y = Y,
         A = A
    )
  
  p0_char <- 
    paste0( 
      round( aux_pars, 2 ),
      collapse = ', ' )
  
  saved_estpars1 <- 
    rbind(
      saved_estpars1,
      data.table(
        p0       = paste0( '(', p0_char, ')' ),
        mu_0     = round( aux_pars[1], 2 ),
        sigma2_0 = round( aux_pars[2], 2 ),
        rho_0    = round( aux_pars[3], 3 ),
        mu       = round( NRest$estimate[ 1 ], 2 ),
        sigma2   = round( NRest$estimate[ 2 ], 2 ), 
        rho      = round( NRest$estimate[ 3 ], 3 ),
        k        = NRest$iterations
      )
    )

}

# Second list: fix mu and sigma in 25 and 25
pars_list2 <- 
  lapply( seq( 0.01, 0.95, 0.01 ),
          function( x ){
            c( 25, 25, x )
          } )

saved_estpars2 <- data.table()

for( i in 1:length( pars_list2 ) ){
  
  aux_pars <- pars_list2[[i]]
  
  NRest <- 
    nlm( mloglik,
         aux_pars,
         hessian = TRUE,
         Y = Y,
         A = A
    )
  
  p0_char <- 
    paste0( 
      round( aux_pars, 2 ),
      collapse = ', ' )
  
  saved_estpars2 <- 
    rbind(
      saved_estpars2,
      data.table(
        p0       = paste0( '(', p0_char, ')' ),
        mu_0     = round( aux_pars[1], 2 ),
        sigma2_0 = round( aux_pars[2], 2 ),
        rho_0    = round( aux_pars[3], 3 ),
        mu       = round( NRest$estimate[ 1 ], 2 ),
        sigma2   = round( NRest$estimate[ 2 ], 2 ), 
        rho      = round( NRest$estimate[ 3 ], 3 ),
        k        = NRest$iterations
      )
    )
  
}

# Third list: fix mu and sigma in 50 and 2
pars_list3 <- 
  lapply( seq( 0.01, 0.95, 0.01 ),
          function( x ){
            c( 50, 2, x )
          } )

saved_estpars3 <- data.table()

for( i in 1:length( pars_list3 ) ){
  
  aux_pars <- pars_list3[[i]]
  
  NRest <- 
    nlm( mloglik,
         aux_pars,
         hessian = TRUE,
         Y = Y,
         A = A
    )
  
  p0_char <- 
    paste0( 
      round( aux_pars, 2 ),
      collapse = ', ' )
  
  saved_estpars3 <- 
    rbind(
      saved_estpars3,
      data.table(
        p0       = paste0( '(', p0_char, ')' ),
        mu_0     = round( aux_pars[1], 2 ),
        sigma2_0 = round( aux_pars[2], 2 ),
        rho_0    = round( aux_pars[3], 3 ),
        mu       = round( NRest$estimate[ 1 ], 2 ),
        sigma2   = round( NRest$estimate[ 2 ], 2 ), 
        rho      = round( NRest$estimate[ 3 ], 3 ),
        k        = NRest$iterations
      )
    )
  
}

# Fourth list: fix rho in 0.01, sigma to sd(Y) and moves mu
pars_list4 <- 
  lapply( seq( 5, 100, 5 ),
          function( x ){
            c( x, 100, 0.01 )
          } )

saved_estpars4 <- data.table()

for( i in 1:length( pars_list4 ) ){
  
  aux_pars <- pars_list4[[i]]
  
  NRest <- 
    nlm( mloglik,
         aux_pars,
         hessian = TRUE,
         Y = Y,
         A = A
    )
  
  p0_char <- 
    paste0( 
      round( aux_pars, 2 ),
      collapse = ', ' )
  
  saved_estpars4 <- 
    rbind(
      saved_estpars4,
      data.table(
        p0       = paste0( '(', p0_char, ')' ),
        mu_0     = round( aux_pars[1], 2 ),
        sigma2_0 = round( aux_pars[2], 2 ),
        rho_0    = round( aux_pars[3], 3 ),
        mu       = round( NRest$estimate[ 1 ], 2 ),
        sigma2   = round( NRest$estimate[ 2 ], 2 ), 
        rho      = round( NRest$estimate[ 3 ], 3 ),
        k        = NRest$iterations
      )
    )
  
}


# Fifth list: fix rho in 0.01, mu to 10 and moves sigma
pars_list5 <- 
  lapply( seq( 5, 200, 10 ),
          function( x ){
            c( 10, x, 0.01 )
          } )

saved_estpars5 <- data.table()

for( i in 1:length( pars_list5 ) ){
  
  aux_pars <- pars_list5[[i]]
  
  NRest <- 
    nlm( mloglik,
         aux_pars,
         hessian = TRUE,
         Y = Y,
         A = A
    )
  
  p0_char <- 
    paste0( 
      round( aux_pars, 2 ),
      collapse = ', ' )
  
  saved_estpars5 <- 
    rbind(
      saved_estpars5,
      data.table(
        p0       = paste0( '(', p0_char, ')' ),
        mu_0     = round( aux_pars[1], 2 ),
        sigma2_0 = round( aux_pars[2], 2 ),
        rho_0    = round( aux_pars[3], 3 ),
        mu       = round( NRest$estimate[ 1 ], 2 ),
        sigma2   = round( NRest$estimate[ 2 ], 2 ), 
        rho      = round( NRest$estimate[ 3 ], 3 ),
        k        = NRest$iterations
      )
    )
  
}

saved_estpars <- 
  rbind( saved_estpars1, saved_estpars2,
         saved_estpars3, saved_estpars4,
         saved_estpars5 )

write.table( saved_estpars,
             row.names = F,
             file = 'atividade3/rmd/estimates_tab.csv' )
###################################

### Plot Map 1 #-------------------
## read shape files
rj_sf <- 
  st_read( 'atividade3/data/Limite_Bairro-shp/Limite_Bairro.shp' )

## merge with count information
map1_dt <-
  rj_sf %>%
  merge(
    atv3_dt[ , .( neighborhood, countClass ) ],
    by.x = 'NOME',
    by.y = 'neighborhood',
    all.x = T
  ) %>%
  mutate( countClass = ifelse( is.na( countClass ),
                               'Sem Dados',
                               paste0( countClass ) ) ) %>%
  mutate( countClass = factor( countClass,
                               levels = c( '0-4', '5-19', '20-49',
                                           '50+', 'Sem Dados' ) ) )
## create plot
map1_desc <- 
  ggplot() +
  geom_sf( data = map1_dt,
           aes( fill = countClass ),
           color = 'gray90',
           size = 0.001 ) +
  scale_fill_manual( values = c( 'Sem Dados' = "gray25", 
                                 '0-4'       = "steelblue", 
                                 '5-19'      = "orange",
                                 '20-49'     = "tomato1", 
                                 '50+'       = "firebrick3" ),
                     name = 'Tiroteios\nRegistrados' ) +
  labs( title = '',
        subtitle = '',
        x = '',
        y = '' ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    strip.background = element_rect( color = 'black',
                                     fill = 'gray90' ),
    strip.text = element_text( size = 15, color = 'black' ),
    legend.text = element_text( size = 12, color = 'black' ),
    legend.title = element_text( size = 13, color = 'black' ),
    axis.ticks = element_blank(),
    axis.text = element_text( size = 10, color = 'black' ),
    plot.title = element_text( size = 17, color = 'black', 
                               face = 'bold' ),
    plot.subtitle = element_text( size = 16, color = 'black' ),
    plot.caption = element_text( size = 13, color = 'black', 
                                 hjust = 0 )
  )

### save descriptive map1
pdf( file = 'atividade3/rmd/map1_desc.pdf',
     height = 5, width = 7 )
print( map1_desc )
dev.off()
###################################

### Construction of map with estimated pars #------

# Construct smoothed estimates with parameters 

n = nrow( A )
I = diag( n )
one <- matrix( 1, ncol = 1, nrow = n ) 
mu = 9.36
sigma = 14.96
rho = 0.012
Yhat <- numeric( n )

for( i in 1:n ){
  S <- sigma^2 * matrix.power( ( I - rho * A ), -2 )
  Yhat[ i ] <- 
    mu +
    S[i,-i] %*% 
    matrix.power( S[-i,-i], -1 ) %*%
    ( Y[-i] - one[-i] * mu )
}

# create data.table with Yhat and labels
atv3_Yhat <- 
  data.table(
    neighbcode = paste0( 'B', 1:n ),
    neighborhood = names(labels),
    Yhat = Yhat 
  ) 


# merge atv3_Yhat data with shape file
map2_dt <-
  rj_sf %>%
  merge(
    atv3_Yhat[ , .( neighborhood, Yhat ) ],
    by.x = 'NOME',
    by.y = 'neighborhood',
    all.x = T
  ) 
## create plot
map2_res <- 
  ggplot() +
  geom_sf( data = map2_dt,
           aes( fill = Yhat ),
           color = 'gray90',
           size = 0.001 ) +
  scale_fill_gradient(  low = 'steelblue',
                        high = 'tomato3',
                        name = 'Tiroteios\nEstimados (SAR)' ) +
  labs( title = '',
        subtitle = '',
        x = '',
        y = '' ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    strip.background = element_rect( color = 'black',
                                     fill = 'gray90' ),
    strip.text = element_text( size = 15, color = 'black' ),
    legend.text = element_text( size = 12, color = 'black' ),
    legend.title = element_text( size = 13, color = 'black' ),
    axis.ticks = element_blank(),
    axis.text = element_text( size = 10, color = 'black' ),
    plot.title = element_text( size = 17, color = 'black', 
                               face = 'bold' ),
    plot.subtitle = element_text( size = 16, color = 'black' ),
    plot.caption = element_text( size = 13, color = 'black', 
                                 hjust = 0 )
  )

### save descriptive map1
pdf( file = 'atividade3/rmd/map2_res.pdf',
     height = 5, width = 7 )
print( map2_res )
dev.off()
###################################
