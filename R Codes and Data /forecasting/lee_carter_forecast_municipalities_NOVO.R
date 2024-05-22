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
pack_list <- c( 'dplyr', 'ggplot2', 'data.table' )
open_install_package(pack_list)

# dictionary for regions
dic_reg <-
  c( 'N'  = '1',
     'NE' = '2',
     'SE' = '3',
     'S'  = '4',
     'CO' = '5',
     'BR' = '0' )
#######################################

### Lee Carter - Municipalities #------

# read bayes adjusted life table data for municipalities
logmxdat_munics <-
  rbind(
    # all municipalities but not DF
    fread( file.path( dir,
                      'DATA/inputs/tabuas_muni_ajustadas_2010_so_adhoc.csv' ),
           encoding = 'Latin-1' ) %>%
      .[,.( municode, muniname,
            sex, age,
            q50 = log( nMx ) ) ],
    fread( file.path( dir,
                      'DATA/inputs/tabuas_brasilia_ajustadas_2010_NOVO.csv' ),
           encoding = 'Latin-1' ) %>%
      .[,.( municode = 530010, 
            muniname = 'Distrito Federal',
            sex, age,
            q50 = logmx.med ) ]
  )

# read kt parameters for regions
kt_par_region <-
  fread( file.path( dir,
                    'DATA/aux_files/lee_carter_kt_parameters_regions.csv' ) ) %>%
  .[, regcode := as.numeric( dic_reg[ as.character( region ) ] ) ]

bxax_par_region <-
  fread( file.path( dir,
                    'DATA/aux_files/lee_carter_bxax_parameters_regions.csv' ) ) %>%
  .[, regcode := as.numeric( dic_reg[ as.character( region ) ] ) ]

munilist <- logmxdat_munics$municode %>% unique
sexlist  <- c( 'm', 'f' )

# average variant estimation
# using lapply

nmuni <- length( munilist )
countmuni <- 1

dat_lc_forecast_munics <-
  do.call(
    rbind,
    lapply(
      sort( munilist ),
      function( x ){

        cat( paste0( '\nMunicipality ',
                     which( x == sort( munilist ) ),
                     ' of ',
                     nmuni ) )

        regcode_sel <- as.numeric( substr( x, 1, 1 ) )

        out_lc_forecast_munics <- data.table()

        for( sex_sel in c( 'm', 'f' ) ){

          ax_sel <- logmxdat_munics[ municode == x & sex == sex_sel,]$q50
          kt_sel <- kt_par_region[ regcode == regcode_sel & sex == sex_sel, ]$kt
          bx_sel <- bxax_par_region[ regcode == regcode_sel & sex == sex_sel, ]$model.bx

          aux <-
            LeeCarterForecast( ax = ax_sel,
                               kt = kt_sel,
                               bx = bx_sel,
                               variant_low  = FALSE,
                               variant_high = FALSE,
                               year_start = 2010,
                               year_end = 2030 ) %>%
            as.data.table %>%
            .[ , list( municode = x,
                       year,
                       sex = sex_sel,
                       age,
                       logmx ) ]

          out_lc_forecast_munics <-
            rbind( out_lc_forecast_munics, aux )

          rm( aux )
        }

        return( out_lc_forecast_munics )
      }
    )
  )

write.csv2( dat_lc_forecast_munics,
            file = file.path( dir,
                              'DATA/outputs/lee_carter_municipalities_average_2010_2030_NOVO.csv' ),
            row.names = FALSE )

gc( reset = T )
gc( reset = T )
gc( reset = T )

# low variant estimation
nmuni <- length( munilist )
countmuni <- 1

dat_lc_forecast_munics_low <-
  do.call(
    rbind,
    lapply(
      sort( munilist ),
      function( x ){
        
        cat( paste0( '\nMunicipality ',
                     which( x == sort( munilist ) ),
                     ' of ',
                     nmuni ) )
        
        regcode_sel <- as.numeric( substr( x, 1, 1 ) )
        
        out_lc_forecast_munics <- data.table()
        
        for( sex_sel in c( 'm', 'f' ) ){
          
          ax_sel <- logmxdat_munics[ municode == x & sex == sex_sel,]$q50
          kt_sel <- kt_par_region[ regcode == regcode_sel & sex == sex_sel, ]$kt
          bx_sel <- bxax_par_region[ regcode == regcode_sel & sex == sex_sel, ]$model.bx
          
          aux <-
            LeeCarterForecast( ax = ax_sel,
                               kt = kt_sel,
                               bx = bx_sel,
                               variant_low  = TRUE,
                               variant_high = FALSE,
                               year_start = 2010,
                               year_end = 2030 ) %>%
            as.data.table %>%
            .[ , list( municode = x,
                       year,
                       sex = sex_sel,
                       age,
                       logmx ) ]
          
          out_lc_forecast_munics <-
            rbind( out_lc_forecast_munics, aux )
          
          rm( aux )
        }
        
        return( out_lc_forecast_munics )
      }
    )
  )

write.csv2( dat_lc_forecast_munics_low,
            file = file.path( dir,
                              'DATA/outputs/lee_carter_municipalities_low_2010_2030_NOVO.csv' ),
            row.names = FALSE )

gc( reset = T )
gc( reset = T )
gc( reset = T )

# high variant estimation
nmuni <- length( munilist )
countmuni <- 1

dat_lc_forecast_munics_high <-
  do.call(
    rbind,
    lapply(
      sort( munilist ),
      function( x ){
        
        cat( paste0( '\nMunicipality ',
                     which( x == sort( munilist ) ),
                     ' of ',
                     nmuni ) )
        
        regcode_sel <- as.numeric( substr( x, 1, 1 ) )
        
        out_lc_forecast_munics <- data.table()
        
        for( sex_sel in c( 'm', 'f' ) ){
          
          ax_sel <- logmxdat_munics[ municode == x & sex == sex_sel,]$q50
          kt_sel <- kt_par_region[ regcode == regcode_sel & sex == sex_sel, ]$kt
          bx_sel <- bxax_par_region[ regcode == regcode_sel & sex == sex_sel, ]$model.bx
          
          aux <-
            LeeCarterForecast( ax = ax_sel,
                               kt = kt_sel,
                               bx = bx_sel,
                               variant_low  = FALSE,
                               variant_high = TRUE,
                               year_start = 2010,
                               year_end = 2030 ) %>%
            as.data.table %>%
            .[ , list( municode = x,
                       year,
                       sex = sex_sel,
                       age,
                       logmx ) ]
          
          out_lc_forecast_munics <-
            rbind( out_lc_forecast_munics, aux )
          
          rm( aux )
        }
        
        return( out_lc_forecast_munics )
      }
    )
  )

write.csv2( dat_lc_forecast_munics_high,
            file = file.path( dir,
                              'DATA/outputs/lee_carter_municipalities_high_2010_2030_NOVO.csv' ),
            row.names = FALSE )
