### Prepare R script #-----------------
rm( list = ls() )
gc( reset = T )

# work dir
dir <- 'C:/Users/Dell/Documents/Backup/backup/bayesiano/Arquivos Victor/ajuste_bigdf_corrigido/projeções do Zé'

# load auxiliar functions
source( file.path( dir, 'R/aux_functions.R'))

# list and load packages
pack_list <- c( 'geobr', 'dplyr', 'ggplot2', 'data.table', 'rcartocolor', 'splines' )
open_install_package(pack_list)


# dictionary for regions
dic_reg <- 
  c('1' = 'N',
    '2' = 'NE',
    '3' = 'SE',
    '4' = 'S',
    '5' = 'CO')

###Banco com as projeções da das logmx por lee-carter

### read bayes lc forecasts (average)
logmx_muni_average <- 
  fread( file.path( dir,
                    'DATA/outputs/lee_carter_municipalities_average_2010_2030_NOVO.csv' ),
         dec = ',' ) %>%
  .[ , mx_average := exp( logmx ) ] %>% 
  .[ , sex2 := ifelse( sex == 'm', 'male', 'female' ) ]

# calculando a e(0)average
require( MortalityLaws )
ltlcdat_muni_average <- 
  do.call(
    rbind,
    lapply( unique( logmx_muni_average$municode ),
            function( muni ){
              
              aux <- 
                logmx_muni_average[ municode == muni ]
              
              out <- data.table()
              
              for( s in c( 'male', 'female' ) ){
                aux2 <- aux[ sex2 == s ] 
                
                for( y in c( 2010:2030 ) ){
                  
                  aux3 <- aux2[ year == y ] 
                  
                  nMx = aux3$mx
                  x   = aux3$age
                  
                  out <- 
                    rbind(
                      out,
                      data.table(
                        municode = muni,
                        year = y,
                        sex  = substr( s, 1, 1 ),
                        MortalityLaws::LifeTable( x = x, 
                                                  mx = nMx,
                                                  sex = s )$lt %>%
                          setDT %>%
                          .[ x == 0,
                             .( e0_average = ex ) ] )
                    )
                  
                }
              }
              
              return( out )
            } )
  )


write.csv2(ltlcdat_muni_average,
           file = file.path( dir,'DATA/outputs/',
                             'e0_2010_2030_municipalities_NOVO.csv'),
           row.names = FALSE)

#Verificação das e0 média do BR por sexo de 2030
e0_2030_1 <-
   fread( file.path( dir,
                    'DATA/outputs/e0_2010_2030_municipalities_NOVO.csv' ),
         dec = ',' ) %>%
  filter(year == 2030) %>% 
  group_by(sex) %>%
  summarise(mean(e0_average))
  




