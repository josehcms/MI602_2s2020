### Codes for Final Presentation - MI602-2s2020
### TOPALS relational model for smoothing death 
### probabilities
### Author: josehcms
### Last update: 2021-01-02

### Housekeeping #----------------------
rm( list = ls( ) )
graphics.off( )

req_packs <- 
  c( 'dplyr', 'data.table', 'lubridate',
     'microdatasus' )
# devtools::install_github("rfsaldanha/microdatasus")
lapply( req_packs, require, character.only = TRUE )

#########################################

### Read pop data #----------------------
dat_popbr <- 
  fread( 'apresentacao/data/pop_br2010_sidra.csv' ) %>%
  .[ , age := ifelse( is.na( age ), 
                      NA,
                      ifelse( age >= 100, 
                              100, age ) ) ] %>%
  .[ , .( pop = sum( pop, na.rm = T ) ),
     .( sex, age ) ]

dat_popgo <- 
  fread( 'apresentacao/data/pop_census_goias.csv' ) %>%
  .[ , municode := as.numeric( substr( municode, 1, 6 ) ) ] %>%
  .[  municode %in% c( 520551, 520140, 520580, 520870, 520025 ) ]  %>%
  .[ , age := ifelse( is.na( age ), 
                      NA,
                      ifelse( age >= 100, 
                              100, age ) ) ] %>%
  .[ , .( pop = sum( pop, na.rm = T ) ),
     .( municode, sex, age ) ]
#########################################

### Read mortality data #----------------

dat_deathbr <- 
  fread( 'apresentacao/data/deaths_sim_br_processed.csv' )  %>%
  .[ , age := ifelse( is.na( age ), 
                      NA,
                      ifelse( age >= 100, 
                              100, age ) ) ] %>%
  .[ , .( deaths = sum( deaths, na.rm = T ) ),
     .( year, sex, age ) ]

dat_deathgo <- 
  fread( 'apresentacao/data/deaths_sim_go_processed.csv' ) %>%
  .[ , age := ifelse( is.na( age ), 
                      NA,
                      ifelse( age >= 100, 
                              100, age ) ) ] %>%
  .[ year == 2010, 
     .( deaths = sum( deaths, na.rm = T ) ),
     .( municode, muniname, year, sex, age ) ]

#########################################

### Prepare standard qx #----------------

std_mx <- 
  dat_deathbr[ , 
               .( deaths = mean( deaths ) ),
               .( sex, age ) ] %>%
  .[ , age := ifelse( is.na( age ), NA,
                      ifelse( age >= 100, 
                              100, 
                              age ) ) ] %>%
  .[ , .( deaths = round( sum( deaths ) ) ),
     .( sex, age ) ] %>%
  merge(
    dat_popbr,
    by = c( 'sex', 'age' ),
    all = T
  ) %>%
  .[ , mx_std := round( deaths / pop, 6 ) ] %>%
  .[ age < 100, 
     .( sex, age, pop_std = pop, mx_std ) ]

#########################################


### Prepare observed qx #----------------

dat_mx_list <- 
  lapply( c( 520551, 520140, 520580, 
             520870, 520025 ),
          function( x ){
          
          aux_deaths <- 
            dat_deathgo[ municode == x & year == 2010 ] %>% copy
          
          aux_pop <- 
            dat_popgo[ municode == x ] %>% copy
          
          name = dat_deathgo[ municode == x & year == 2010 ]$muniname %>% unique
          
          aux <- 
            merge(
              aux_pop,
              aux_deaths,
              by  = c( 'municode', 'sex', 'age' ),
              all = T
            ) %>%
            .[ , .( sex, 
                    age, 
                    pop_obs = ifelse( is.na( pop ),
                                      0, 
                                      pop ), 
                    mx_obs = ifelse( is.na( deaths / pop ),
                                     0,
                                     deaths / pop ) ) ]
          out <- 
            std_mx %>%
            merge(
              aux,
              by = c( 'sex', 'age' ),
              all = T
            ) %>%
            .[ , mx_obs := ifelse( is.na( mx_obs ),
                                   0,
                                   mx_obs ) ] %>%
            .[ , pop_obs := ifelse( is.na( pop_obs ),
                                    0,
                                    pop_obs ) ] %>%
            .[ age < 100,
               .( municode = x,
                  muniname = name,
                  sex, age,
                  pop_std,
                  mx_std,
                  pop_obs,
                  mx_obs ) ]
          
          return( out )
        } )

dat_mx <- 
  do.call( rbind,
           dat_mx_list )

write.table( dat_mx,
             file = 'apresentacao/data/topals_input.csv',
             row.names = F )
#########################################

### Plots 1 #----------------------------
require(ggplot2)

plot1 <- 
  dat_mx[ municode %in% c( 520551, 520870 ) & sex == 'm' ] %>%
  copy %>%
  .[ , sex := ifelse( sex == 'm', 
                      'Homens', 
                      'Mulheres' ) ]%>%
  .[ , muniname := ifelse( municode == 520551, 
                           'Cocalzinho de Goiás', 'Goiânia' ) ] %>%
  ggplot( ) +
  geom_point( size = 1.25, 
              aes( x = age, y = mx_obs ) ) +
  geom_line( size = 0.75,
             aes( x = age, y = mx_std ),
             linetype = 'dotted' ) +
  scale_y_log10( name = 'log(mx)') +
  scale_x_continuous( breaks = seq( 0, 100, 10 ),
                      name = 'Idade' ) +
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
    axis.text  = element_text( color = 'black', size = 12 ) 
  )

ggsave( 'apresentacao/rmd/motivacao.png',
        height = 3, width = 8,
        plot = plot1 )

#########################################
