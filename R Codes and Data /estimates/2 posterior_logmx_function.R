library(tidyverse)
library(rstan)
library(ggplot2)
library(magrittr)
library(dplyr)

rm(list=ls())

graphics.off()
if (.Platform$OS.type == 'windows') windows(record=TRUE)


geo = read.csv('C:/Users/Dell/OneDrive/Material GitHuh - Paper PHM/Input Data/IBGE geography.csv', stringsAsFactors = FALSE)

## municipio-level geographic data (names and codes for state, meso, micro, etc.)
muni.url = 'http://schmert.net/topals-mortality/data/muni.df.csv'
muni.df  = read.csv(muni.url , header=TRUE, stringsAsFactors = FALSE)

ibge = read.csv(text='
                UF,male_e0,female_e0
                RO,67.0,73.6
                AC,68.5,75.2
                AM,67.3,73.8
                RR,66.9,72.5
                PA,67.5,74.6
                AP,69.2,75.4
                TO,68.7,74.9
                MA,65.1,72.8
                PI,66.1,73.8
                CE,68.5,76.4
                RN,70.2,78.0
                PB,67.4,75.0
                PE,66.8,75.5
                AL,64.6,74.0
                SE,66.9,75.2
                BA,67.7,76.4
                MG,72.5,78.3
                ES,72.9,79.5
                RJ,70.4,77.6
                SP,72.6,79.2
                PR,72.0,78.6
                SC,73.7,79.9
                RS,72.6,79.2
                MS,70.4,77.6
                MT,69.5,75.9
                GO,70.1,76.4
                DF,72.5,79.7
                ', as.is=TRUE)



# trapez approx of life expectancy from a logmx schedule over ages 0..99
e0 = function(logmx) {
  mx = exp(logmx)
  px = exp(-mx)
  lx = c(1,cumprod(px))
  return( sum(head(lx,-1) + tail(lx,-1)) / 2)
}

file_info = expand.grid(sex=c('Female','Male'), uf= head(geo$abbrev,-1), stringsAsFactors = FALSE) %>%
  mutate( filename=paste0('Stan-23abr22-1955-micro-within-', 
                          uf,'-',sex,'.RData'))

more_info = expand.grid(sex=c('Female','Male'), uf= 'Brazil', stringsAsFactors = FALSE) %>%
  mutate( filename=paste0('Stan-23abr22-1955-state-within-Brazil-', sex,'.RData'))

file_info = rbind(file_info,
                  more_info)

### Início da Função ###

plot_TOPALS_fit = function(this.uf, this.sex, output='onscreen') {
  
  if (this.uf != 'Brazil') {  
    tmp = muni.df %>%
      filter(ufabb == this.uf) %>%
      group_by(microcode) %>%
      summarize(microname = microname[1])
    
    region.name = tmp$microname
  } else {
    tmp = muni.df %>%
      group_by(ufcode) %>%
      summarize(ufname = ufname[1])
    region.name = tmp$ufname
  }
  
  fitfile = (file_info %>% filter(uf==this.uf, sex==this.sex))$filename
  
  load(file=fitfile)
  
  
  ## distributions under assumptions of full and partial death registration
  
  R = stanDataList$R
  K = stanDataList$K
  B = stanDataList$B
  N = stanDataList$N
  D = stanDataList$D
  
  a = as.matrix(fit, 'alpha')  # sims x (K*R)
  a = array(a,  c(nrow(a), K, R))
  
  lambda_star = stanDataList$std_schedule
  
  if (output == 'pdf') pdf(file=paste0(paste('TOPALS_fit',this.uf,this.sex,sep="-"),'.pdf'))
  
  Topals_fit = as.data.frame(matrix(NA, 0, 8))
  
  for (r in 1:R) {
    L = lambda_star + B %*% t(a[,,r])
    
    e = apply(L,2,e0)
    e10 = sprintf("%.1f",   round( quantile(e, .10),1))
    e90 = sprintf("%.1f",   round( quantile(e, .90),1))
    
    q10    = apply(L, 1, quantile,.10)
    q50    = apply(L, 1, quantile,.50)    
    q90    = apply(L, 1, quantile,.90)
    
  if (this.uf!="Brazil") {
    df      = data.frame(microcode = rep(tmp$microcode[r], 100),
                         microname = rep(tmp$microname[r], 100),
                         age   = 0:99,
                         logDN = log(D[,r] / N[,r]),
                         q10   = q10,
                         q50   = q50,
                         q90   = q90,
                         zero  = (D[,r] == 0))
  } else {
    df      = data.frame(ufcode = rep(tmp$ufcode[r], 100),
                         ufname = rep(tmp$ufname[r], 100),
                         age   = 0:99,
                         logDN = log(D[,r] / N[,r]),
                         q10   = q10,
                         q50   = q50,
                         q90   = q90,
                         zero  = (D[,r] == 0))
  }
    
    Topals_fit = rbind(Topals_fit, df)
    
    hue = c(Male='royalblue', Female='deeppink')
    
    G1 = ggplot(data=df, aes(x=age, y=logDN)) +
      xlab('Age') + ylab('Log mortality rate') +
      scale_x_continuous(breaks=seq(0,100,20)) +
      scale_y_continuous(breaks=seq(-10,1,1),minor_breaks=NULL,limits=c(-10,1)) +
      geom_point(data=filter(df, !zero), shape='+', size=4, color=grey(.25)) +
      geom_point(aes(x=age, y=q50), color=hue[this.sex])+
      geom_segment(aes(x=age, y=q10, xend=age, yend=q90), color=hue[this.sex]) +
      geom_rug(data=filter(df,zero)) +
      theme_bw() +
      theme(title=element_text(face='bold')) +
      labs(title= paste(region.name[r], this.uf, this.sex, 
                        '\n 80% interval for e0 = [', e10, '-', e90, ']'))
    
    print(G1)
    
  } # for r
  
  if (output == 'pdf') dev.off()
  
  return(Topals_fit)
  
} # plot_TOPALS_fit

### Fim da Função ###


# Resultados para as UFs:
BR = plot_TOPALS_fit("Brazil", "Female", 'pdf')
write.csv(BR, file='Female_Bayes_result_by_UF.csv', row.names=FALSE)

BR = plot_TOPALS_fit("Brazil", "Male", 'pdf')
write.csv(BR, file='Male_Bayes_result_by_UF.csv', row.names=FALSE)

# Resultados para as micros:
MG = plot_TOPALS_fit("MG", "Female", output='onscreen')
MG = plot_TOPALS_fit("MG", "Male", output='onscreen')

# Resultados para as micros:
RN = plot_TOPALS_fit("RN", "Female", 'pdf')
RN = plot_TOPALS_fit("RN", "Male", 'pdf')



#--------------------
#logmx por UF
#--------------------

UF.fits.f = plot_TOPALS_fit('Brazil', "Female")
UF.fits.f$sex = 'f'
UF.fits.m = plot_TOPALS_fit('Brazil', "Male")
UF.fits.m$sex = 'm'
UF.fits = rbind(UF.fits.f, UF.fits.m)
UF.fits = UF.fits %>%
  select(ufcode, ufname, sex, age, logDN, q10, q50, q90)
write.csv(UF.fits, file='logmx.uf.2010.csv', row.names=FALSE)

#--------------------
#logmx por micro
#--------------------
uf.abb = unique(muni.df$ufabb)
flogmx.micro.f = as.data.frame(matrix(NA, 0, 9))
for (this.uf in uf.abb) {
  
  tmp = plot_TOPALS_fit(this.uf, "Female")
  tmp$sex = 'f'
  flogmx.micro.f = rbind(flogmx.micro.f, tmp)
  
}

flogmx.micro.m = as.data.frame(matrix(NA, 0, 9))
for (this.uf in uf.abb) {
  
  tmp = plot_TOPALS_fit(this.uf, "Male")
  tmp$sex = 'm'
  flogmx.micro.m = rbind(flogmx.micro.m, tmp)
  
}

Micro.fits = rbind(flogmx.micro.f, flogmx.micro.m)
Micro.fits = Micro.fits %>%
  select(microcode, microname, sex,  age, logDN, q10, q50, q90)
write.csv(Micro.fits, file='logmx.micro.2010.csv', row.names=FALSE)

