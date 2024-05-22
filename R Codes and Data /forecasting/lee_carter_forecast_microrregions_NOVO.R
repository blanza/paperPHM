#######################################
### Script: Lee-Carter for micros
### Author: Jose H C Monteiro da Silva
### Last update: 2021-06-04
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

### Lee Carter - Microrregions #------

# read bayes adjusted life table data for microrregions
logmxdat_micros <-
  rbind(
    # all micros but not DF
    fread( file.path( dir,
                      'DATA/inputs/tabuas_micro_ajustadas_2010_so_ADHOC.csv' ) ) %>%
      .[,.( microcode, microname,
            sex, age,
            q50 = log( mx ) ) ],
    fread( file.path( dir,
                      'DATA/inputs/tabuas_brasilia_ajustadas_2010_NOVO.csv' ),
           encoding = 'Latin-1' ) %>%
      .[,.( microcode = 53001, 
            microname = 'Distrito Federal',
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

microlist <- logmxdat_micros$microcode %>% unique
sexlist  <- c( 'm', 'f' )

# average variant estimation
dat_lc_forecast_micros <- data.table()

nmicro <- length( microlist )
countmicro <- 1
for( microcode_sel in microlist ){
  cat( paste0( '\nMicrorregion ', countmicro, ' of ', nmicro ) )

  regcode_sel <- as.numeric( substr( microcode_sel, 1, 1 ) )

  for( sex_sel in sexlist ){

    ax_sel <- logmxdat_micros[ microcode == microcode_sel & sex == sex_sel,]$q50
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
      .[ , list( microcode = microcode_sel,
                 year,
                 sex = sex_sel,
                 age,
                 logmx ) ]

    dat_lc_forecast_micros <-
      rbind( dat_lc_forecast_micros, aux )

    rm( aux )
  }
  countmicro <- countmicro + 1
}

write.csv2( dat_lc_forecast_micros,
            file = file.path( dir,
                              'DATA/outputs/lee_carter_microrregions_average_2010_2030_NOVO.csv' ),
            row.names = FALSE )

