### Codes for Final Presentation - MI602-2s2020
### TOPALS relational model for smoothing death 
### probabilities
### Author: josehcms
### Last update: 2021-01-04

### Housekeeping #----------------------
rm( list = ls( ) )
graphics.off( )

req_packs <- 
  c( 'dplyr', 'data.table', 'ggplot2', 'splines' )

lapply( req_packs, require, character.only = TRUE )
#########################################

### Read Data #--------------------------
topals_input <- 
  fread( 'apresentacao/data/topals_input.csv' )

#########################################

### Setup TOPALS function #--------------

# Source: https://github.com/schmert/TOPALS
topals_fit <- 
  function(
    mx, Nx, 
    age = 0 : 99,
    mx_std,
    age_knots = c( 0, 1, 10, 20, 40, 70 )
  ){
    
    D = rpois( 100, Nx * mx )
    
    logmx_std = log( mx_std )
    B = bs( age, knots = age_knots, degree = 1 ) # linear B-spline basis
    
    # Penalized likelihood function
    Q = function(alpha) {
      logmx.hat = as.numeric( logmx_std + B %*% alpha)
      penalty = sum(diff(alpha)^2)
      return( sum(D * logmx.hat - Nx * exp(logmx.hat)) - penalty)
    }
    
    ## expected deaths function
    Dhat = function(alpha) {
      lambda.hat = logmx_std + B %*% alpha
      return( as.numeric( Nx * exp(lambda.hat) ))
    }
    
    ## S matrix for penalty
    S = matrix( 0, 6, 7 )
    diag(S[,1:6]) = -1
    diag(S[,2:7]) = +1
    SS = crossprod(S)
    
    next_alpha = function(alpha) {
      dhat = Dhat(alpha)
      M = solve ( t(B) %*% diag(dhat) %*% B + 2*SS)
      v = t(B) %*% (D - dhat) - 2* (SS %*% alpha)
      return( alpha + M %*% v)
    }
    
    a = matrix(0, 7, 10)
    for (i in 2:ncol(a)) { a[,i] = next_alpha(a[,i-1])}
    # round(a, 4)
    
    Q.iter = apply( a, 2, Q )
    fitted_logmx = logmx_std + B %*% a[,ncol(a)]
    
    return( fitted_logmx )
  }

#########################################


### Run TOPALS #--------------------------

set.seed( 123456 )

topals_fitted_list <- 
  lapply(
    topals_input$municode %>% unique,
    function( x ){
      mx     = topals_input[ municode == x & sex == 'm' ]$mx_obs
      mx_std = topals_input[ municode == x & sex == 'm' ]$mx_std
      Nx     = topals_input[ municode == x & sex == 'm' ]$pop_obs
      name   = topals_input[ municode == x & sex == 'm' ]$muniname %>% unique
      
      fit_topals <- 
        data.table(
          municode = x,
          muniname = name,
          age = 0:99,
          mx_fit = exp( topals_fit( mx = mx, Nx = Nx, mx_std = mx_std ) ),
          mx_obs = mx,
          mx_std = mx_std
        )
      
      return( fit_topals )
      
    }
  )

topals_fitted <- 
  do.call( 
    rbind,
    topals_fitted_list 
    )

### Save plots #----------------
resplot1 <- 
  topals_fitted[ muniname %in% 
                   c( 'Corumbá de Goiás', 'Cocalzinho de Goiás' ) ] %>%
  melt(
    id.vars = c( 'municode', 'muniname', 'age', 'mx_obs' ),
    measure.vars = c( 'mx_fit.V1', 'mx_std' ),
    value.name = 'mx',
    variable.name = 'type'
  ) %>%
  .[ , muniname := factor( muniname,
                           levels = c( 'Corumbá de Goiás', 'Cocalzinho de Goiás',
                                       'Águas Lindas de Goiás',
                                       'Aparecida de Goiânia', 'Goiânia' ) ) ] %>%
  ggplot( ) +
  geom_point( aes( x = age, y = mx_obs ),
              shape = 1 ) +
  geom_line( aes( x = age, y = mx, 
                  color = type, 
                  linetype = type,
                  size = type ) ) +
  scale_y_log10( name = 'log(mx)' ) +
  scale_color_manual( values = c( 'tomato3', 'skyblue' ),
                      labels = c( 'mx_std' = 'Brasil 2010',
                                  'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_linetype_manual( values = c( 'solid', 'longdash' ),
                         labels = c( 'mx_std' = 'Brasil 2010',
                                     'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_size_manual( values = c( 1, 0.75 ),
                     labels = c( 'mx_std' = 'Brasil 2010',
                                 'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_x_continuous( breaks = seq( 0, 100, 10 ),
                      name = '' ) +
  facet_wrap( ~ muniname ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line( linetype = 'longdash',
                                     color = 'gray25',
                                     size = 0.05 ),
    panel.grid.minor = element_line( linetype = 'longdash',
                                     color = 'gray25',
                                     size = 0.05 ),
    axis.title = element_text( color = 'black', size = 13 ),
    axis.text  = element_text( color = 'black', size = 10 ),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text( color = 'black', size = 13 )
  )

ggsave( filename = 'apresentacao/rmd/restopals1.png',
        width = 6, height = 3 )

resplot2 <- 
  topals_fitted[ muniname %in% 
                   c( 'Águas Lindas de Goiás' ) ] %>%
  melt(
    id.vars = c( 'municode', 'muniname', 'age', 'mx_obs' ),
    measure.vars = c( 'mx_fit.V1', 'mx_std' ),
    value.name = 'mx',
    variable.name = 'type'
  ) %>%
  .[ , muniname := factor( muniname,
                           levels = c( 'Corumbá de Goiás', 'Cocalzinho de Goiás',
                                       'Águas Lindas de Goiás',
                                       'Aparecida de Goiânia', 'Goiânia' ) ) ] %>%
  ggplot( ) +
  geom_point( aes( x = age, y = mx_obs ),
              shape = 1 ) +
  geom_line( aes( x = age, y = mx, 
                  color = type, 
                  linetype = type,
                  size = type ) ) +
  scale_y_log10( name = 'log(mx)' ) +
  scale_color_manual( values = c( 'tomato3', 'skyblue' ),
                      labels = c( 'mx_std' = 'Brasil 2010',
                                  'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_linetype_manual( values = c( 'solid', 'longdash' ),
                         labels = c( 'mx_std' = 'Brasil 2010',
                                     'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_size_manual( values = c( 1, 0.75 ),
                     labels = c( 'mx_std' = 'Brasil 2010',
                                 'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_x_continuous( breaks = seq( 0, 100, 10 ),
                      name = '' ) +
  facet_wrap( ~ muniname ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line( linetype = 'longdash',
                                     color = 'gray25',
                                     size = 0.05 ),
    panel.grid.minor = element_line( linetype = 'longdash',
                                     color = 'gray25',
                                     size = 0.05 ),
    axis.title = element_text( color = 'black', size = 13 ),
    axis.text  = element_text( color = 'black', size = 10 ),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text( color = 'black', size = 13 )
  )

ggsave( filename = 'apresentacao/rmd/restopals2.png',
        width = 4, height = 3 )

resplot3 <- 
  topals_fitted[ muniname %in% 
                   c( 'Aparecida de Goiânia', 'Goiânia' ) ] %>%
  melt(
    id.vars = c( 'municode', 'muniname', 'age', 'mx_obs' ),
    measure.vars = c( 'mx_fit.V1', 'mx_std' ),
    value.name = 'mx',
    variable.name = 'type'
  ) %>%
  .[ , muniname := factor( muniname,
                           levels = c( 'Corumbá de Goiás', 'Cocalzinho de Goiás',
                                       'Águas Lindas de Goiás',
                                       'Aparecida de Goiânia', 'Goiânia' ) ) ] %>%
  ggplot( ) +
  geom_point( aes( x = age, y = mx_obs ),
              shape = 1 ) +
  geom_line( aes( x = age, y = mx, 
                  color = type, 
                  linetype = type,
                  size = type ) ) +
  scale_y_log10( name = 'log(mx)' ) +
  scale_color_manual( values = c( 'tomato3', 'skyblue' ),
                      labels = c( 'mx_std' = 'Brasil 2010',
                                  'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_linetype_manual( values = c( 'solid', 'longdash' ),
                         labels = c( 'mx_std' = 'Brasil 2010',
                                     'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_size_manual( values = c( 1, 0.75 ),
                     labels = c( 'mx_std' = 'Brasil 2010',
                                 'mx_fit.V1' = 'Ajuste TOPALS' ) ) +
  scale_x_continuous( breaks = seq( 0, 100, 10 ),
                      name = '' ) +
  facet_wrap( ~ muniname ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line( linetype = 'longdash',
                                     color = 'gray25',
                                     size = 0.05 ),
    panel.grid.minor = element_line( linetype = 'longdash',
                                     color = 'gray25',
                                     size = 0.05 ),
    axis.title = element_text( color = 'black', size = 13 ),
    axis.text  = element_text( color = 'black', size = 10 ),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text( color = 'black', size = 13 )
  )

ggsave( filename = 'apresentacao/rmd/restopals3.png',
        width = 6, height = 3 )
#######################################

### Compute e0 #-----------------------
e0 = function(mx) {
  px = exp(-mx)
  lx = c(1,cumprod(px))
  return( sum(head(lx,-1) + tail(lx,-1)) / 2)
}

topals_fitted[,lapply( .SD, e0 ), by = .( muniname ), .SDcols = 'mx_fit.V1']
topals_fitted[,lapply( .SD, e0 ), by = .( muniname ), .SDcols = 'mx_obs']

#######################################
