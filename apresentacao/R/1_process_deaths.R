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

### Function to prepare age_sex data format #------
process_sexage_deaths <- 
  function( dat_sim ){

    setDT( dat_sim )
    
    # adjust age in years
    dat_sim[ , age := ifelse( !is.na( IDADEhoras ) | 
                                !is.na( IDADEmeses ) | 
                                !is.na( IDADEdias ),
                              0,
                              ifelse( !is.na( IDADEanos ),
                                      IDADEanos,
                                      NA ) ) %>%
               as.numeric ]
    
    # adjust sex labels
    dat_sim[ , sex := ifelse( SEXO == 'Feminino', 
                              'f', 
                              ifelse( SEXO == 'Masculino', 
                                      'm', 
                                      NA ) ) ]
    # add year lab
    dat_sim[ , year := year( DTOBITO ) ]
    
    # prepare output
    out_sim <- 
      dat_sim[ , 
               list(
                 deaths = .N 
                 ),
               .( year, sex, age ) 
               ]
    
    return( out_sim )
  }
#########################################

### Read and process data for 2009-2011 standard #---------
dat_sim_list <- 
  lapply( 
    2009 : 2011, 
    function( x ){
    
      dat_sim <- 
        fetch_datasus( year_start = x, 
                       year_end = x,
                       information_system = "SIM-DO" ) %>%
        process_sim %>%
        process_sexage_deaths
      
      return( dat_sim )
      
  } )

dat_sim <- 
  do.call( rbind,
           dat_sim_list ) %>%
  setorder( year, sex, age ) %>%
  .[ !is.na( age ) & !is.na( sex ) ]

write.table( dat_sim,
             file = 'apresentacao/data/deaths_sim_br_processed.csv',
             row.names = FALSE )

#########################################

### Function to prepare age_sex data format - specific loc #------
process_sexagemuni_deaths <- 
  function( dat_sim ){
    
    setDT( dat_sim )
    
    # adjust age in years
    dat_sim[ , age := ifelse( !is.na( IDADEhoras ) | 
                                !is.na( IDADEmeses ) | 
                                !is.na( IDADEdias ),
                              0,
                              ifelse( !is.na( IDADEanos ),
                                      IDADEanos,
                                      NA ) ) %>%
               as.numeric ]
    
    # adjust sex labels
    dat_sim[ , sex := ifelse( SEXO == 'Feminino', 
                              'f', 
                              ifelse( SEXO == 'Masculino', 
                                      'm', 
                                      NA ) ) ]
    # add year lab
    dat_sim[ , year := year( DTOBITO ) ]
    
    # add municode
    dat_sim[ , municode := CODMUNRES ]
    dat_sim[ , muniname := munResNome ]
    
    # prepare output
    out_sim <- 
      dat_sim[ , 
               list(
                 deaths = .N 
               ),
               .( municode, muniname, year, sex, age ) 
      ]
    
    return( out_sim )
  }
#########################################


### Read and process data for 2009-2011 - specific loc #---------
dat_sim_list <- 
  lapply( 
    2009 : 2011, 
    function( x ){
      
      dat_sim <- 
        fetch_datasus( year_start = x, 
                       year_end = x,
                       uf = 'GO',
                       information_system = "SIM-DO" ) %>%
        process_sim %>%
        process_sexagemuni_deaths
      
      return( dat_sim )
      
    } )

dat_sim <- 
  do.call( rbind,
           dat_sim_list ) %>%
  setorder( municode, year, sex, age ) %>%
  .[ !is.na( age ) & !is.na( sex ) ] %>%
  .[ municode %in% c( 520551, 520140, 520580, 520870,
                      520025 ) ]

write.table( dat_sim,
             file = 'apresentacao/data/deaths_sim_go_processed.csv',
             row.names = FALSE )

#########################################
