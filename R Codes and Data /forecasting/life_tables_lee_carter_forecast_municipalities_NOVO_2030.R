#################
# Construção das tábuas de vida projetadas para 2030 a partir das projeções das Mx (average, high e low)
#################


### Prepare R script #-----------------
rm( list = ls() )
gc( reset = T )

# work dir
dir <- 'C:/Users/Dell/Documents/Backup/backup/bayesiano/Arquivos Victor/ajuste_bigdf_corrigido/projeções do Zé'
outputs_dir <- 'C:/Users/Dell/Documents/Backup/backup/bayesiano/Arquivos Victor/ajuste_bigdf_corrigido/projeções do Zé/DATA/outputs'

# load auxiliar functions
source( file.path( dir, 'R/aux_functions.R'))

# list and load packages
pack_list <- c( 'geobr', 'dplyr', 'ggplot2', 'data.table', 'readr' )
open_install_package(pack_list)


# Ler as estimativas das taxas feitas pelo Lee-Carter (AVERAGE)
lee_carter_logmx_average <- 
  fread( file.path( dir,
                    'DATA/outputs/lee_carter_municipalities_average_2010_2030_NOVO.csv'),
         dec = ',' ) %>% 
  filter(year == 2030) %>% 
  mutate(mx = exp(logmx))


# Construção da tábua de vida projetada 2023 (average)
tabuas_muni <- list()
contador = 1

# life tables are constructed for each sex from each year by each muni:
for ( muni in unique(lee_carter_logmx_average$municode)){
  for ( sexo in c('m','f')){
    
    lt_data <- 
      lee_carter_logmx_average[municode==muni & sex==sexo] 
    
    
    tabua_vida <-  
      life.table( x = lt_data$age,
                  nMx = lt_data$mx) %>%
      .[,`:=`(
        municode   = muni,
        sex         = sexo
      )]
    
    tabuas_muni[[contador]] <- 
      tabua_vida
    
    contador <-
      contador+1
  }
}


tabuas_muni <- 
  rbindlist(tabuas_muni)

write.csv2(tabuas_muni,
           file = file.path(outputs_dir,
                            "lt_2030_average_bayes_ahdoc_NOVO.csv"),
           row.names = F)

#Verificando a média da e0 do BR por sexo
e0_2030_2 <- 
  fread( file.path( dir,
                    'DATA/outputs/lt_2030_average_bayes_ahdoc_NOVO.csv' ),
         dec = ',' ) %>% 
  filter(age == 0) %>% 
  group_by(sex) %>% 
  summarise(mean(ex))



############

# Ler as estimativas das taxas feitas pelo Lee-Carter (HIGH)
lee_carter_logmx_high <- 
  fread( file.path( dir,
                    'DATA/outputs/lee_carter_municipalities_high_2010_2030_NOVO.csv'),
         dec = ',' ) %>% 
  filter(year == 2030) %>% 
  mutate(mx = exp(logmx))


# Construção da tábua de vida projetada 2023 (HIGH)
tabuas_muni <- list()
contador = 1

# life tables are constructed for each sex from each year by each muni:
for ( muni in unique(lee_carter_logmx_high$municode)){
  for ( sexo in c('m','f')){
    
    lt_data <- 
      lee_carter_logmx_high[municode==muni & sex==sexo] 
    
    
    tabua_vida <-  
      life.table( x = lt_data$age,
                  nMx = lt_data$mx) %>%
      .[,`:=`(
        municode   = muni,
        sex         = sexo
      )]
    
    tabuas_muni[[contador]] <- 
      tabua_vida
    
    contador <-
      contador+1
  }
}


tabuas_muni <- 
  rbindlist(tabuas_muni)

write.csv2(tabuas_muni,
           file = file.path(outputs_dir,
                            "lt_2030_high_bayes_ahdoc.csv"),
           row.names = F)


############

# Ler as estimativas das taxas feitas pelo Lee-Carter (LOW)
lee_carter_logmx_low <- 
  fread( file.path( dir,
                    'DATA/outputs/lee_carter_municipalities_low_2010_2030_NOVO.csv'),
         dec = ',' ) %>% 
  filter(year == 2030) %>% 
  mutate(mx = exp(logmx))


# Construção da tábua de vida projetada 2023 (LOW)
tabuas_muni <- list()
contador = 1

# life tables are constructed for each sex from each year by each muni:
for ( muni in unique(lee_carter_logmx_low$municode)){
  for ( sexo in c('m','f')){
    
    lt_data <- 
      lee_carter_logmx_low[municode==muni & sex==sexo] 
    
    
    tabua_vida <-  
      life.table( x = lt_data$age,
                  nMx = lt_data$mx) %>%
      .[,`:=`(
        municode   = muni,
        sex         = sexo
      )]
    
    tabuas_muni[[contador]] <- 
      tabua_vida
    
    contador <-
      contador+1
  }
}


tabuas_muni <- 
  rbindlist(tabuas_muni)

write.csv2(tabuas_muni,
           file = file.path(outputs_dir,
                            "lt_2030_low_bayes_ahdoc.csv"),
           row.names = F)
