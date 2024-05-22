library(dplyr)
library(rstan)
library(withr)

rm(list=ls())

graphics.off()
if (.Platform$OS.type == 'windows') windows(record=TRUE)

## set up muni data to get MICROREGION NAMES in the same order as in .RData files
## (this was an omission in the main model.R program: should have kept the 
##  region names)


## Observed 2010 deaths and exposures for municipios
big.df <- read.csv("new.big.df.csv", head = TRUE, stringsAsFactors = FALSE) %>% 
  select(-X)


## municipio-level geographic data (names and codes for state, meso, micro, etc.)
muni.url = 'http://schmert.net/topals-mortality/data/muni.df.csv'
muni.df  = read.csv(muni.url , header=TRUE, stringsAsFactors = FALSE)


BR      = left_join( big.df, muni.df, by='munid') %>% 
  group_by(microcode) %>%
  summarize(ufcode=ufcode[1],
            ufabb=ufabb[1],
            ufname=ufname[1],
            mesocode=mesocode[1],
            mesoname=mesoname[1])


########################


ufcodes = c(ro=11 , ac=12 , am=13 , rr=14 , pa=15 , ap=16 , to=17 , ma=21 , 
            pi=22 , ce=23 , rn=24 , pb=25 , pe=26 , al=27 , se=28 , ba=29 , 
            mg=31 , es=32 , rj=33 , sp=35 , pr=41 , sc=42 , rs=43 , ms=50 , 
            mt=51 , go=52 , df=53)

ufnames = as.character(c('Rond?nia', 'Acre','Amazonas','Roraima','Par?', 'Amap?', 
                         'Tocantins', 'Maranh?o','Piau?', 'Cear?', 'Rio Grande do Norte', 
                         'Para?ba','Pernambuco', 'Alagoas', 'Sergipe', 'Bahia',
                         'Minas Gerais','Esp?rito Santo','Rio de Janeiro','S?o Paulo',
                         'Paran?','Santa Catarina', 'Rio Grande do Sul','Mato Grosso do Sul',
                         'Mato Grosso', 'Goi?s','Distrito Federal'))

ufabbs = toupper(names(ufcodes))

uf.df = data.frame( ufabb=ufabbs, 
                    ufcode=ufcodes, 
                    ufname=ufnames, 
                    stringsAsFactors = FALSE) %>%
  arrange(ufcodes)                        # sort for later merging

## for each prior there is a set of .RData output files
##  named Stan-(TIMESTAMP)-micro-within-(UFABB)-(SEX).RData
##  and two output files
##  named Stan-(TIMESTAMP)-state-within-Brazil-(SEX).RData

sex_vals = c('Female', 'Male')

uf_vals = head(ufabbs, -1)   # don't include DF (only 1 microregion)

## microregion file information
info = expand.grid(sex   = sex_vals,
                   uf    = uf_vals,
                   stringsAsFactors = FALSE) %>%
  mutate(    filename  = paste0('Stan-23abr22-1955-micro-within-',
                                uf,'-',sex,'.RData'))

## append state file information
more_info = expand.grid(sex   = sex_vals,
                        uf    = 'Brazil',
                        stringsAsFactors = FALSE) %>%
  mutate(         filename  = paste0('Stan-23abr22-1955-state-within-Brazil-',
                                     sex,'.RData'))

file_info = rbind(info,
                  more_info)



########################
# trapez approx of life expectancy from a logmx schedule over ages 0..99
e0 = function(logmx) {
  mx = exp(logmx)
  px = exp(-mx)
  lx = c(1,cumprod(px))
  return( sum(head(lx,-1) + tail(lx,-1)) / 2)
}


uf.list = c('DF', head( unique(file_info$uf) , -1))    # DF, plus all UFs except 'Brazil'

result = data.frame()




for (this.sex in c('Male','Female')) {
  
  for (this.uf in uf.list) {
    
    print(this.uf)
    
    ## get the male e0 information 
    if (this.uf != 'DF') {
      fitfile = (file_info %>% filter(uf==this.uf, sex==this.sex))$filename
    } else {
      fitfile = (file_info %>% filter(uf=='Brazil', sex==this.sex))$filename
    }
    
    load(file=fitfile)
    
    ## microregion codes in the same order as .RData
    these.microregions = unique( filter(BR, ufabb==this.uf)$microcode)
    
    ## distributions under assumptions of full and partial death registration
    R = stanDataList$R
    K = stanDataList$K
    B = stanDataList$B
    lambda_star = stanDataList$std_schedule
    
    a = as.matrix(fit, 'alpha')  # sims x (K*R)
    a = array(a,  c(nrow(a), K, R))
    
    a_full = array(ML$theta_tilde,  c(nrow(ML$theta_tilde), K, R))
    
    ## which small regions have the microregional estimates (DF is a special case: state within Brazil)
    if (this.uf != 'DF') {
      ix = 1:R
    } else {
      ix = which(colnames(stanDataList$N) == '53')
    }
    
    
    ## loop over this state's microregions
    for (r in ix) {
      
      if (this.uf != 'DF') {
        this.microcode = these.microregions[r]
      } else {
        this.microcode = 53001
      }
      
      print(this.microcode)
      
      L = lambda_star + B %*% t(a[,,r])
      e = apply(L,2,e0)
      
      Q = round( quantile( e, probs=c(.05,.10,.25,.50,.75,.90,.95)), 2)
      
      L_full = lambda_star + B %*% t(a_full[,,r])
      e_full = apply(L_full,2,e0)
      
      Q_full = round( quantile( e_full, probs=c(.05,.10,.25,.50,.75,.90,.95)), 2)
      
      tmp = expand.grid(ufabb           = this.uf, 
                        microcode       = this.microcode,
                        sex             = tolower(substr(this.sex,1,1)),    # 'm' or 'f'
                        pctile          = c(5,10,25,50,75,90,95),
                        coverage_model  = c('partial','full'))
      
      tmp$e0  = c(Q, Q_full)  # fill in the column (partial,full) x pctile
      
      result = rbind(result,
                     tmp)
      
    } # for r
    
    
    
  } # for this.uf
  
} # for this.sex

write.csv(result, file='e0_microregion_summary.csv', row.names=FALSE)
