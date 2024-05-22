#----------------------------------------
#Applying TOPALS_fit to municipios
#---------------------------------------

library(dplyr)
library(splines)
library(data.table)

rm(list = ls())

graphics.off()
if (.Platform$OS.type == 'windows') windows(record=TRUE)


## Observed 2010 deaths and exposures for municipios
big.df <- read.csv("new.big.df.csv", head = TRUE, stringsAsFactors = FALSE) %>% 
  select(-X) 

## municipio-level geographic data (names and codes for state, meso, micro, etc.)
muni.url = 'http://schmert.net/topals-mortality/data/muni.df.csv'
muni.df  = read.csv(muni.url , header=TRUE, stringsAsFactors = FALSE)


#### Topals fit fuction ######

TOPALS_fit = function( N, D, std,
                       max_age        = 99,
                       knot_positions = c(0,1,10,20,40,70), 
                       smoothing_k    = 1,
                       max_iter       = 20,
                       alpha_tol      = .00005,
                       details        = FALSE) {
  
  require(splines)
  
  ## single years of age from 0 to max_age
  age = 0:max_age
  
  ## B is an Ax7 matrix. Each column is a linear B-spline basis function
  B      = splines::bs( age, knots=knot_positions, degree=1 )
  nalpha = ncol(B) 
  
  ## penalized log lik function
  Q = function(alpha) {
    lambda.hat = as.numeric( std + B %*% alpha)
    penalty    = smoothing_k * sum( diff(alpha)^2 )
    return( sum(D * lambda.hat - N * exp(lambda.hat)) - penalty)
  }
  
  ## expected deaths function
  Dhat = function(alpha) {
    lambda.hat = std + B %*% alpha
    return(  as.numeric( N * exp(lambda.hat) ))
  }      
  
  ## S matrix for penalty
  S = matrix(0,nalpha-1,nalpha) 
  diag(S[, 1:(nalpha-1)]) = -1
  diag(S[, 2:(nalpha)  ]) = +1
  SS = crossprod(S)
  
  #------------------------------------------------
  # iteration function: 
  # next alpha vector as a function of current alpha
  #------------------------------------------------
  next_alpha = function(alpha) {
    dhat = Dhat(alpha)
    M = solve ( t(B) %*% diag(dhat) %*% B + 2*smoothing_k *SS)
    v = t(B) %*% (D - dhat) - 2* (smoothing_k * (SS %*% alpha))
    return( alpha + M %*% v)
  }
  
  ## main iteration:     
  a = rep(0, nalpha)
  
  niter = 0
  repeat {
    niter      = niter + 1
    last_param = a
    a          = next_alpha( a )  # update
    change     = a - last_param
    
    converge = all( abs(change) < alpha_tol)
    overrun  = (niter == max_iter)
    
    if (converge | overrun) { break }
    
  } # repeat
  
  if (details | !converge | overrun) {
    if (!converge) print('did not converge')
    if (overrun) print('exceeded maximum number of iterations')
    
    dhat = Dhat(a)
    covar = solve( t(B) %*% diag(dhat) %*% B + 2*smoothing_k *SS)
    
    return( list( alpha    = a, 
                  covar    = covar,
                  Qvalue   = Q(a),
                  converge = converge, 
                  maxiter  = overrun))
  } else return( a) 
  
} # TOPALS_fit



## calculate a "Brazil standard" by sex: result is a 100 x 2 matrix with
##   female and male stds in 'f' and 'm' columns, respectively
tmp = big.df %>%
  group_by(sex,age) %>%
  summarize(lambda = log( sum(death)/sum(expos)))

BR = matrix(tmp$lambda, ncol=2, dimnames=list(0:99, c('f','m')))

## smoothed version: project onto a cubic spline with closely-spaced knots
basis = bs(0:99, knots=seq(0,99,2))
Proj  = basis %*% solve(crossprod(basis)) %*% t(basis)
BRstd = Proj %*% BR

## check standard:
plot(BR[,2]) # 1 = female # 2 = male
lines(BRstd[,2])




# IN?CIO DO LOOPING: rodar para cada sexo

# Par?metors necess?rios:
age = 0:99
B   = bs( 0:99, knots=c(0,1,10,20,40,70), degree=1 )
muni.list = unique(muni.df$municode)
dados <- left_join(big.df, muni.df, by='munid')
#this.std = BRstd[,1] # 1 = female e 2 = male
sex.id = c('f', 'm')

BRstd = as.data.frame(BRstd)


## Topals fit por sexo: ---------------------------------------------------------------------

#Homens

topals.fit = as.data.frame(matrix(NA, 0, 12)) #data frame para guardar os resultados do looping

system.time({
  for(this.muni in muni.list){
   
      tmp = dados %>%
        filter(sex=='m', municode==this.muni) %>% #mudar o sexo
        group_by(ufname, microcode, microname, municode, muniname, sex, age) %>%
        summarise(death = sum(death), expos = sum(expos))
      
      tmp = as.data.frame(tmp)
      
      N = tmp$expos
      D = tmp$death
      this.std = BRstd %>% select(m) #mudar o sexo
      this.std = as.matrix(this.std)
      
      fit = TOPALS_fit(N, D, this.std, details=TRUE)
      
      fitted_logmx = this.std + B %*% fit$a
      
      se_logmx = sqrt( diag (B %*% fit$covar %*% t(B)) )
      logmx.sup = fitted_logmx -1.96 * se_logmx
      logmx.inf = fitted_logmx +1.96 * se_logmx
      
      result = cbind(tmp, logmx.inf, fitted_logmx, logmx.sup)
      
      topals.fit = rbind(topals.fit, result)
      
      print(this.muni)

  }
})

colnames(topals.fit) <- c("ufname","microcode","microname","municode",
                       "muniname","sex","age","death","expos",    
                       "logmx.inf","logmx.med", "logmx.sup")
topals.fit.male <- topals.fit
write.csv(topals.fit.male,"topals.fit.male.csv")
# write_rds(topals.fit.male,"topals.fit.male.rds")




#Mulheres

topals.fit = as.data.frame(matrix(NA, 0, 12)) #data frame para guardar os resultados do looping

system.time({
  for(this.muni in muni.list){
    
    tmp = dados %>%
      filter(sex=='f', municode==this.muni) %>% #mudar o sexo
      group_by(ufname, microcode, microname, municode, muniname, sex, age) %>%
      summarise(death = sum(death), expos = sum(expos))
    
    tmp = as.data.frame(tmp)
    
    N = tmp$expos
    D = tmp$death
    this.std = BRstd %>% select(f) #mudar o sexo
    this.std = as.matrix(this.std)
    
    fit = TOPALS_fit(N, D, this.std, details=TRUE)
    
    fitted_logmx = this.std + B %*% fit$a
    
    se_logmx = sqrt( diag (B %*% fit$covar %*% t(B)) )
    logmx.sup = fitted_logmx -1.96 * se_logmx
    logmx.inf = fitted_logmx +1.96 * se_logmx
    
    result = cbind(tmp, logmx.inf, fitted_logmx, logmx.sup)
    
    topals.fit = rbind(topals.fit, result)
    
    print(this.muni)
    
  }
})

colnames(topals.fit) <- c("ufname","microcode","microname","municode",
                          "muniname","sex","age","death","expos",    
                          "logmx.inf","logmx.med", "logmx.sup")
topals.fit.female <- topals.fit
write.csv(topals.fit.female,"topals.fit.female.csv")
# write_rds(topals.fit.female,"topals.fit.female.rds")

#Topals ambos os sexos

topals.fit.both <- rbind(topals.fit.male, topals.fit.female)

write.csv(topals.fit.both,"topals.fit.2sexos.csv")
# write_rds(topals.fit.both,"topals.fit.2sexos.rds")


# Estimativas municipais - Topals corrigido pelos fatores de correção das micro --------
# Obs.: Foi necessário primeiro calcular os graus de cobertura das micro
# (arquivo: Graus de cobertura - microrregião) para depois seguir nessas estimativas.

gc_micro <- read.csv("graus.cobertura.micro.csv", sep = ",")
fit.muni <- read.csv("topals.fit.2sexos.csv", sep = ',')

fit.muni2 <- left_join(fit.muni,gc_micro, by = c("microcode", "sex", "age")) %>% 
  select(ufname.x,microcode,microname.x,municode,
         muniname,sex,age,death.x,expos.x,logmx.inf.x,
         logmx.med.x,logmx.sup.x, grau_cobertura)
colnames(fit.muni2) <- c("ufname","microcode","microname","municode","muniname","sex",
                         "age","death","expos","logmx.inf","logmx.med","logmx.sup",   
                         "grau_cobertura")
head(fit.muni2)

# Calculando as taxas de mortalidade municipais com corre??o de sub-registro

fit.muni2 <- fit.muni2 %>% 
  filter(ufname != "Distrito Federal") %>%   
  mutate(mx.corrigida.sup = exp(logmx.inf)*(1/grau_cobertura)) %>% 
  mutate(mx.corrigida.med = exp(logmx.med)*(1/grau_cobertura)) %>% 
  mutate(mx.corrigida.inf = exp(logmx.sup)*(1/grau_cobertura)) %>% 
  mutate(logmx.corrigida.sup = log(mx.corrigida.inf)) %>% 
  mutate(logmx.corrigida.med = log(mx.corrigida.med)) %>% 
  mutate(logmx.corrigida.inf = log(mx.corrigida.sup))


write.csv(fit.muni2, "taxa.mortalidade.muni.corrigidas.csv")

# Funçao para gerar as funçõe da tábua de vida ---------------------------

life.table <- function( x, nMx){
  # simple lifetable using Keyfitz and Flieger separation factors and 
  # exponential tail of death distribution (to close out life table)
  b0 <- 0.07;   b1<- 1.7;      
  nmax <- length(x)
  #nMx = nDx/nKx   
  n <- c(diff(x), 999)          		  # width of the intervals
  nax <- n/2;		            	        # default to .5 of interval
  nax[1] <- b0 + b1 * nMx[1]    		  # from Keyfitz & Flieger(1968)
  nax[nmax] <- 1/nMx[nmax] 	  	      # e_x at open age interval
  nqx <- (n * nMx) / (1 + (n - nax) * nMx)
  nqx<-ifelse(nqx > 1, 1, nqx);		    # necessary for high nMx
  nqx[nmax] <- 1.0
  lx <- c(1, cumprod(1 - nqx));   	  # survivorship lx
  lx <- lx[1:length(nMx)]
  ndx <- lx * nqx;
  nLx <- n * lx - nax * ndx;      	 # equivalent to n*l(x+n) + (n-nax)*ndx
  nLx[nmax] <- lx[nmax] * nax[nmax]
  Tx <- rev(cumsum(rev(nLx)))
  ex <- ifelse( lx[1:nmax] > 0, Tx/lx[1:nmax], NA);
  lt <- data.frame(Ages = x, nqx = nqx, lx = lx, ndx = ndx, nLx = nLx, Tx = Tx, ex = ex, nMx = nMx)
  return(lt)
}

mx.muni.f = fit.muni2 %>% filter(sex == 'f')
mx.muni.m = fit.muni2 %>% filter(sex == 'm')

muni.list = unique(fit.muni2$municode)
age = c(0:99)

lt.muni.f = as.data.frame(matrix(NA, 0, 56))
lt.muni.m = as.data.frame(matrix(NA, 0, 56))

## life tables all municipalities:

## Masculino

system.time({
  
  for(this.muni in muni.list) {
    
    tmp = mx.muni.m %>% # mudar sex
      filter(municode == this.muni)
    
    #  colnames(tmp) = c('microcode',  'microname', 'sex', 'age', 'logmx.obs', 'logmx10', 'logmx50', 
    #                    'logmx90', 'mx10', 'mx50', 'mx90')
    
    mx.inf = tmp$mx.corrigida.inf
    mx.med = tmp$mx.corrigida.med
    mx.sup = tmp$mx.corrigida.sup
    
    # mx10 = tmp$stdmx10.muni
    # mx50 = tmp$stdmx50.muni
    # mx90 = tmp$stdmx90.muni
    
    lt.sup = life.table(age, mx.inf)
    lt.sup = lt.sup %>% select(nqx, lx, ndx, nLx, Tx, ex)
    colnames(lt.sup) = c('nqx.inf', 'lx.sup', 'ndx.inf', 'nLx.sup', 'Tx.sup', 'ex.sup')
    lt.med = life.table(age, mx.med)
    lt.med = lt.med %>% select(nqx, lx, ndx, nLx, Tx, ex)
    colnames(lt.med) = c('nqx.med', 'lx.med', 'ndx.med', 'nLx.med', 'Tx.med', 'ex.med')
    lt.inf = life.table(age, mx.sup)
    lt.inf = lt.inf %>% select(nqx, lx, ndx, nLx, Tx, ex)
    colnames(lt.inf) = c('nqx.sup', 'lx.inf', 'ndx.sup', 'nLx.inf', 'Tx.inf', 'ex.inf')
    
    # lt90 = life.table(age, mx10)
    # lt90 = lt90 %>% select(nqx, lx, ndx, nLx, Tx, ex)
    # colnames(lt90) = c('nqx10', 'lx90', 'ndx10', 'nLx90', 'Tx90', 'ex90')
    # lt50 = life.table(age, mx50)
    # lt50 = lt50 %>% select(nqx, lx, ndx, nLx, Tx, ex)
    # colnames(lt50) = c('nqx50', 'lx50', 'ndx50', 'nLx50', 'Tx50', 'ex50')
    # lt10 = life.table(age, mx90)
    # lt10 = lt10 %>% select(nqx, lx, ndx, nLx, Tx, ex)
    # colnames(lt10) = c('nqx90', 'lx10', 'ndx90', 'nLx10', 'Tx10', 'ex10')
    
    tmp = as.data.frame(tmp)
    result = cbind(tmp, lt.inf, lt.med, lt.sup)
    lt.muni.m = rbind(lt.muni.m, result) # mudar "lt.muni.m" para lt.muni.f ou vice-versa
    lt.muni.m.final = lt.muni.m
  }
  
})

write.csv(lt.muni.m.final, 'lt.all.muni.2010.males.csv')
#saveRDS(lt.muni.m.final, "lt.all.muni.2010.males.rds")


## Feminino

system.time({
  
  for(this.muni in muni.list) {
    
    tmp = mx.muni.f %>% # mudar sex
      filter(municode == this.muni)
    
    #  colnames(tmp) = c('microcode',  'microname', 'sex', 'age', 'logmx.obs', 'logmx10', 'logmx50', 
    #                    'logmx90', 'mx10', 'mx50', 'mx90')
    
    mx.inf = tmp$mx.corrigida.inf
    mx.med = tmp$mx.corrigida.med
    mx.sup = tmp$mx.corrigida.sup
    
    # mx10 = tmp$stdmx10.muni
    # mx50 = tmp$stdmx50.muni
    # mx90 = tmp$stdmx90.muni
    
    lt.sup = life.table(age, mx.inf)
    lt.sup = lt.sup %>% select(nqx, lx, ndx, nLx, Tx, ex)
    colnames(lt.sup) = c('nqx.inf', 'lx.sup', 'ndx.inf', 'nLx.sup', 'Tx.sup', 'ex.sup')
    lt.med = life.table(age, mx.med)
    lt.med = lt.med %>% select(nqx, lx, ndx, nLx, Tx, ex)
    colnames(lt.med) = c('nqx.med', 'lx.med', 'ndx.med', 'nLx.med', 'Tx.med', 'ex.med')
    lt.inf = life.table(age, mx.sup)
    lt.inf = lt.inf %>% select(nqx, lx, ndx, nLx, Tx, ex)
    colnames(lt.inf) = c('nqx.sup', 'lx.inf', 'ndx.sup', 'nLx.inf', 'Tx.inf', 'ex.inf')
    
    # lt90 = life.table(age, mx10)
    # lt90 = lt90 %>% select(nqx, lx, ndx, nLx, Tx, ex)
    # colnames(lt90) = c('nqx10', 'lx90', 'ndx10', 'nLx90', 'Tx90', 'ex90')
    # lt50 = life.table(age, mx50)
    # lt50 = lt50 %>% select(nqx, lx, ndx, nLx, Tx, ex)
    # colnames(lt50) = c('nqx50', 'lx50', 'ndx50', 'nLx50', 'Tx50', 'ex50')
    # lt10 = life.table(age, mx90)
    # lt10 = lt10 %>% select(nqx, lx, ndx, nLx, Tx, ex)
    # colnames(lt10) = c('nqx90', 'lx10', 'ndx90', 'nLx10', 'Tx10', 'ex10')
    
    tmp = as.data.frame(tmp)
    result = cbind(tmp, lt.inf, lt.med, lt.sup)
    lt.muni.f = rbind(lt.muni.f, result) # mudar "lt.muni.m" para lt.muni.f ou vice-versa
    lt.muni.f.final = lt.muni.f
  }
  
})

write.csv(lt.muni.f.final, 'lt.all.muni.2010.females.csv')
#saveRDS(lt.muni.f.final, "lt.all.muni.2010.females.rds")

