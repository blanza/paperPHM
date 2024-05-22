#######################################
### Script: Lee-Carter for RBE article 
### Author: Jose H C Monteiro da Silva
### Last update: 2020-06-22
#######################################


### Prepare R script #-----------------
rm( list = ls() )
gc( reset = T )

# work dir
dir <- 'C:/Users/Dell/Documents/Backup/backup/bayesiano/Arquivos Victor/ajuste_bigdf_corrigido/projeções do Zé'

# load auxiliar functions
source( file.path( dir, 'R/aux_functions.R'))

# list and load packages
pack_list <- c( 'geobr', 'dplyr', 'ggplot2', 'data.table' )
open_install_package(pack_list)

# dictionary for regions
dic_reg <- 
  c( 'N'  = '1',
     'NE' = '2',
     'SE' = '3',
     'S'  = '4',
     'CO' = '5',
     'BR' = '0' )

dic_uf <- 
  c(  '11' = 'RO',
      '12' = 'AC',
      '13' = 'AM',
      '14' = 'RR',
      '15' = 'PA',
      '16' = 'AP',
      '17' = 'TO',
      '21' = 'MA',
      '22' = 'PI',
      '23' = 'CE',
      '24' = 'RN',
      '25' = 'PB',
      '26' = 'PE',
      '27' = 'AL',
      '28' = 'SE',
      '29' = 'BA',
      '31' = 'MG',
      '32' = 'ES',
      '33' = 'RJ',
      '35' = 'SP',
      '41' = 'PR',
      '42' = 'SC',
      '43' = 'RS',
      '50' = 'MS',
      '51' = 'MT',
      '52' = 'GO',
      '53' = 'DF' )
#######################################

### Lee Carter Life Tables - Municipalities #------

for( i in 1:5 ){
  # read Lee Carter Estimation outputs
  lee_carter_logmx <- 
    fread( file.path( dir,
                      'DATA/outputs/lee_carter_municipalities_low_2010_2030_NOVO_so_ahdoc.csv'),
           dec = ',' )
  
  lee_carter_logmx[, uf := as.numeric( substr( municode, 1, 2 ) ) ]
  
  lee_carter_logmx <- lee_carter_logmx[ as.numeric( substr( uf, 1, 1 ) ) == i, ]
  
  uflist <- lee_carter_logmx$uf %>% unique
  
  for( uf_sel in uflist ){
    
    ufabb <- dic_uf[ as.character( uf_sel ) ]
    
    # run life tables for each sex and municipality
    ltdat_forecast <- data.table()
    age <- 0:99
    munilist <- lee_carter_logmx[uf == uf_sel,]$municode %>% unique
    sexlist <- c( 'm', 'f' )
    yearlist <- lee_carter_logmx[uf == uf_sel,]$year %>% unique
    countmuni <- 1
    
    nmuni <- length( munilist )
    countmuni <- 1
    
    for( municode_sel in munilist ){
      gc( reset = T )
      cat( paste0( '\nMunicipality ', countmuni, ' of ', nmuni ) )
      
      for( sex_sel in sexlist ){
        
        for( year_sel in yearlist ){
          mx_sel <- exp( lee_carter_logmx[ sex == sex_sel & year == year_sel & 
                                             uf == uf_sel &
                                             municode == municode_sel,]$logmx )
          aux <- 
            life.table( x = age, nMx = mx_sel ) %>%
            .[ , list( municode = municode_sel,
                       year = year_sel,
                       sex = sex_sel,
                       age, 
                       nMx, nqx, lx, ndx, nLx, Tx, ex ) ]
          
          ltdat_forecast <- 
            rbind( ltdat_forecast, aux )
          
          rm( aux )
        }
      }
      countmuni <- countmuni + 1
    }
    
    write.csv2( 
      ltdat_forecast, 
      file = file.path( dir,
                        paste0( 'DATA/outputs/',
                                'life_tables_forecasts_by_uf_municipalities_NOVO_so_ahdoc/',
                                'life_tables_lee_carter_municipalities_low_2010_2030_NOVO_so_ahdoc',
                                ufabb,
                                '.csv' ) ),
      row.names = FALSE) 
    
  }
  
}


#######################################
