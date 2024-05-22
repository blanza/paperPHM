rm(list=ls())
gc(reset=T)

require(data.table)
require(dplyr)
require(readr)

# working dir #-------
work_dir   <- 'C:/Users/Dell/Documents/Backup/backup/bayesiano/Arquivos Victor/ajuste_bigdf_corrigido'
inputs_dir <- 'C:/Users/Dell/Documents/Backup/backup/bayesiano/Arquivos Victor/ajuste_bigdf_corrigido'
outputs_dir <- 'C:/Users/Dell/Documents/Backup/backup/bayesiano/Arquivos Victor/ajuste_bigdf_corrigido'
#####################

### 1) Carrega bases de dados #--------------------------

# funcao tabua de vida:
source(file.path(work_dir,
                 "funcao_tabua_vida.R"))

# base gerada pelo topals:
data_logmx_muni <- 
  rbind(
    fread(file.path(inputs_dir,
                    'lt.all.muni.2010.females.csv'),
          encoding = 'Latin-1') %>%
      .[,sex:='f'],
    fread(file.path(inputs_dir,
                    'lt.all.muni.2010.males.csv'),
          encoding = 'Latin-1') %>%
      .[,sex:='m']
  ) %>%
  .[,mx:=mx.corrigida.med] %>%
  .[,pop:=expos] 

# dados micro - bayesiano (DE MARCOS E CARL)

data_logmx_micro <-
  rbind(
    fread(file.path(inputs_dir,
                    'female.logmx.micro.csv'),  
          sep=',', 
          encoding = 'Latin-1') %>%
      .[,sex:='f'],
    fread(file.path(inputs_dir,
                    'male.logmx.micro.csv'),
          sep=',', 
          encoding = 'Latin-1') %>%
      .[,sex:='m']) %>%
  .[,.(microcode,
       age,
       sex,
       logmx=q50,
       logmx.sup=q90)]

# dados UF - bayesiano (DE MARCOS E CARL)
data_logmx_uf <- 
  rbind(
    fread(file.path(inputs_dir,
                    'Female_Bayes_result_by_UF.csv'), 
          sep=',',
          encoding = 'Latin-1') %>%
      .[,sex:='f'],
    fread(file.path(inputs_dir,
                    'Male_Bayes_result_by_UF.csv'), 
          sep=',',
          encoding = 'Latin-1') %>%
      .[,sex:='m']
  ) %>%
  .[,.(ufcode,
       age,
       sex,
       logmx=q50,
       logmx.sup=q90)]

# dicionario uf, municipio, micro, meso:
muni.url <-  
  'http://schmert.net/topals-mortality/data/muni.df.csv'

muni.df  <- 
  fread(muni.url,
        encoding='Latin-1') %>%
  .[,.(ufcode,mesocode,microcode,municode)]

# nomes micro e nomes municipios
names_micro <- 
  data_logmx_muni[,.(microcode,microname)] %>%
  unique

names_muni <- 
  data_logmx_muni[,.(ufname,microcode,microname,municode,muniname)] %>%
  unique
########################################################

###########################################################

### 5. ADHOC - MICROS ########################

### 5.1 Roda tabua de vida para micros #------------------
# life tables para micros:
tabuas_micro <- list()
contador = 1

for ( micro in unique(data_logmx_micro$microcode)){
  for ( sexo in c('m','f')){
    
    lt_data <- 
      data_logmx_micro[microcode==micro & sex==sexo] %>%
      .[,mx := exp(logmx)]%>%
      .[,mx_sup:=exp(logmx.sup)]
    
    
    tabua_vida <-  
      cbind(
        life.table( x   = lt_data$age,
                    nMx = lt_data$mx),
        life.table( x   = lt_data$age,
                    nMx = lt_data$mx_sup) %>%
          .[,.(nMx.sup = nMx,
               ex.inf  = ex)]) %>%
      .[,`:=`(
        microcode   = micro,
        sex         = sexo
      )]
    
    tabuas_micro[[contador]] <- 
      tabua_vida
    
    contador <-
      contador+1
  }
}

tabuas_micro <- 
  rbindlist(tabuas_micro)
##########################################################

### 5.2 Identifica outliers entre as micros #--------

e0_data <- 
  tabuas_micro[age==0]

e0_list <- 
  list()

contador = 1
for(sexo in unique(e0_data$sex)){
  e0_data_new <- 
    e0_data[sex==sexo] 
  
  # compute Q1 and Q3
  Q1_Q3 <- 
    quantile(e0_data_new$ex, 
             probs = c(0.25,0.75))
  
  # compute IQR
  IQR <- 
    diff(Q1_Q3)
  
  # define outlier upperbound
  out_up <- 
    IQR*1.5+Q1_Q3[2]
  
  e0_data_new[,outlier := ifelse(ex > out_up,
                                 1,
                                 0)]
  
  e0_data_new[,out_ub := IQR*1.5+Q1_Q3[2]]
  e0_data_new[,Q3     := Q1_Q3[2]]
  
  e0_list[[contador]] <- e0_data_new
  
  contador <- 
    contador + 1
}

fator_list <-
  rbindlist(e0_list) %>%
  .[,.(microcode,sex,outlier,out_ub,Q3,ex)]

# micros sao outliers independentemente do sexo, para ajustar diferentes patamares!
fator_list_unsex <-
  fator_list[,.(outlier=sum(outlier)),
             .(microcode)] %>%
  unique

# adiciona à base
# tabuas_micro <- 
#   merge(tabuas_micro,
#         fator_list[,.(microcode,sex,outlier)],
#         by=c('microcode','sex'))

tabuas_micro <-
  merge(tabuas_micro,
        fator_list_unsex,
        by='microcode',
        all.x=T) 
##########################################################

### 5.3 Para as micros outliers, faz mx.micro = mx.uf #------

data_logmx_micro_outliers <- 
  merge(tabuas_micro[outlier>0,
                     .(ufcode=as.numeric( substr( microcode , 1 , 2 ) ),
                       microcode,
                       sex,
                       age,
                       outlier)],
        data_logmx_uf[,.(ufcode,sex,age,logmx)],
        by=c('ufcode','sex','age'),
        all.x=T)


# roda novamente a tabua de vida somente para os outliers
tabuas_micro_outliers <- list()
contador = 1

for ( micro in unique(data_logmx_micro_outliers$microcode)){
  for ( sexo in c('m','f')){
    
    lt_data <- 
      data_logmx_micro_outliers[microcode==micro & sex==sexo] %>%
      .[,mx := exp(logmx)]
    
    
    tabua_vida <-  
      life.table( x   = lt_data$age,
                  nMx = lt_data$mx) %>%
      .[,`:=`(
        microcode   = micro,
        sex         = sexo
      )]
    
    tabuas_micro_outliers[[contador]] <- 
      tabua_vida
    
    contador <-
      contador+1
  }
}

tabuas_micro_outliers <- 
  rbindlist(tabuas_micro_outliers)

# rbind dos outliers com os nao outliers
tabuas_micro_ajustadas <- 
  rbind(tabuas_micro[outlier==0,
                     .(microcode,sex,age,mx=nMx,ex,outlier)],
        tabuas_micro_outliers[,.(microcode,sex,age,mx=nMx,ex,outlier=1)])

# adiciona nomes:
tabuas_micro_ajustadas <- 
  merge(tabuas_micro_ajustadas,
        names_micro,
        by='microcode',
        all.x=T)

write.table(tabuas_micro_ajustadas,
            file = file.path(outputs_dir,
                             'tabuas_micro_ajustadas_2010_so_AHDOC.csv'),
            row.names = F)

###############################################

### 6. ADHOC - MUNICIPIOS ########################

### 6.1 Identifica outliers SEM bayes empirico #------
#tabuas_muni <- 
  #bayes_uf_outliers[,.(municode,sex,age,mx,ex)]

tabuas_muni <- 
  data_logmx_muni[,.(municode,sex,age,mx,ex.med)] %>% 
  rename(ex = ex.med)

e0_data <- 
  tabuas_muni[age==0]

e0_list <- 
  list()

contador = 1
for(sexo in unique(e0_data$sex)){
  e0_data_new <- 
    e0_data[sex==sexo] 
  
  # compute Q1 and Q3
  Q1_Q3 <- 
    quantile(e0_data_new$ex, 
             probs = c(0.25,0.75))
  
  # compute IQR
  IQR <- 
    diff(Q1_Q3)
  
  # define outlier upperbound
  out_up <- 
    IQR*1.5+Q1_Q3[2]
  
  e0_data_new[,outlier := ifelse(ex > out_up,
                                 1,
                                 0)]
  
  e0_data_new[,out_ub := IQR*1.5+Q1_Q3[2]]
  e0_data_new[,Q3     := Q1_Q3[2]]
  
  e0_list[[contador]] <- e0_data_new
  
  contador <- 
    contador + 1
}

fator_list <-
  rbindlist(e0_list) %>%
  .[,.(municode,sex,outlier,out_ub,Q3,ex)]

# teste: torna a lista de outliers por micro, independente do sexo -> sexos com o mesmo patamar de mortalidade!
fator_list_unsex <-
  fator_list[,.(outlier=sum(outlier)),
             .(municode)] %>%
  unique

#Colocando a informação da grande região
fator_list_unsex_reg <- fator_list_unsex %>% 
  mutate(reg = case_when( municode < 200000 ~ "NO",
                          municode >= 200000 & municode < 300000 ~ "NE",
                          municode >= 300000 & municode < 400000 ~ "SE",
                          municode >= 400000 & municode < 500000 ~ "SU",
                          municode >= 500000 & municode < 530000 ~ "CO",
                          municode >= 530000 ~"DF"))


write.csv2(fator_list_unsex_reg,
            file = file.path(outputs_dir,
                             'municipios_outliers_list.csv'))

#Quantidade de municípios outliers em cada grande região (o resultado 2 é quando o município teve a estimativa outlier para os dois sexos)
outliers_by_region <- table(fator_list_unsex_reg$reg, fator_list_unsex_reg$outlier)
write.csv2(outliers_by_region,
           file = file.path(outputs_dir,
                            'outliers_by_region.csv'))


# adiciona à base
# tabuas_muni <- 
#   merge(tabuas_muni,
#         fator_list[,.(municode,sex,outlier)],
#         by=c('municode','sex'))


tabuas_muni <-
  merge(tabuas_muni,
        fator_list_unsex,
        by='municode',
        all.x=T) %>%
  merge(muni.df,
        by='municode')

# Salvando a lista de municípios por sexo e idade com identificação outliers
outliers_muni_sex_age <- tabuas_muni %>% 
  select(municode, sex, age, outlier)

write.csv2(outliers_muni_sex_age,
           file = file.path(outputs_dir,
                            'outliers_muni_sex_age.csv'))

###############################################

### 6.2 Para os municipios outliers, faz mx.muni = mx.micro #-------
data_logmx_muni_outliers <- 
  merge(tabuas_muni[outlier>0,
                    .(microcode,
                      municode,
                      sex,
                      age,
                      outlier)],
        tabuas_micro_ajustadas[,.(microcode,sex,age,mx)],
        by=c('microcode','sex','age'),
        all.x=T)


# reestima tabuas de vida municipais:
data_muni <- 
  rbind(
    tabuas_muni[outlier==0,
                .(municode,sex,age,mx,outlier)],
    data_logmx_muni_outliers[,.(municode,sex,age,mx,outlier)])

tabuas_muni <- list()
contador = 1

# life tables are constructed for each sex from each year by each muni:
for ( muni in unique(data_muni$municode)){
  for ( sexo in c('m','f')){
    
    lt_data <- 
      data_muni[municode==muni & sex==sexo] 
    
    
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

# adiciona nomes dos municipios:
tabuas_muni <- 
  merge(tabuas_muni,
        names_muni,
        by='municode',
        all.x = T)

write.table(tabuas_muni,
            file = file.path(outputs_dir,
                             'tabuas_muni_ajustadas_2010_so_adhoc.csv'),
            row.names = F)

###############################################

### Juntando as estimativas das TV sem ajuste de Valter (topalz), uma variável identificando se o município foi outlier, as estimativas da TV municipais após o ajuste AH-DOC e as estimativas intervalares das micro
lt_completa <- data_logmx_muni %>% 
  left_join(outliers_muni_sex_age, by = c("municode", "sex", "age")) %>% 
  left_join(tabuas_muni, by = c("municode", "sex", "age")) %>% 
  select(-V1, -ufname.y, -microcode.y, -microname.y, -muniname.y) %>% 
  rename(ufname = ufname.x, microcode = microcode.x, muniname = muniname.x, microname = microname.x,
         nqx_ahdoc = nqx, lx_ahdoc = lx, ndx_ahdoc = ndx, nLx_ahdoc = nLx, 
         Tx_ahdoc = Tx, ex_ahdoc = ex, nMx_ahdoc = nMx)

#Lendo as estimastivas pontual e intervalares das micros
logmx_micro <-
  rbind(
    fread(file.path(inputs_dir,
                    'female.logmx.micro.csv'),  
          sep=',', #mudei aqui. antes estava ";" em vez de ","
          encoding = 'Latin-1') %>%
      .[,sex:='f'],
    fread(file.path(inputs_dir,
                    'male.logmx.micro.csv'),
          sep=',', #mudei aqui. antes estava ";" em vez de ","
          encoding = 'Latin-1') %>%
      .[,sex:='m']) %>%
  .[,.(microcode,
       age,
       sex,
       Mx.micro.sup=exp(q90),
       Mx.micro.inf=exp(q10))]


lt_completa <-  lt_completa %>% 
  left_join(logmx_micro, by = c("microcode", "sex", "age"))

write.csv2(lt_completa,
           file = file.path(outputs_dir,
                            'lt_completa.csv'))

###############################################

tabuas_muni <- fread("tabuas_muni_ajustadas_2010_so_adhoc.csv")


data_plot_muni <- 
  tabuas_muni[age==0] %>%
  copy %>%
  .[,sex:=ifelse(sex=='f','Mulheres','Homens')]

require(ggplot2)
ggplot(data=data_plot_muni) +
  geom_boxplot(aes(x=0,y=ex))+
  facet_wrap(~sex)+
  theme_bw()+
  labs(title = 'Boxplot - e0 Municípios',
       y = 'e0',
       x = '')


ggsave(filename = file.path(work_dir,
                            'e0_boxplot_2010_muni_ajustado_so_adhoc_NOVO.png'),
       dpi    = 300,
       width  = 10,
       height = 6)
