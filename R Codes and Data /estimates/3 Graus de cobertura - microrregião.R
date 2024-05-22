
############ Graus de cobertura - microrregião V2 ##################
############ A proposta é obter os graus de cobertura pelo quociente entre as mx.micro.topals puro e mx.micro.Bayesiano

# Passo 1: Calcular mx.micro.topals puro (taxas de mortalidade pelo modelo Topals puro) -------------

library(tidyverse)
library(magrittr)
library(splines)
library(dplyr)
library(tidyr)


rm(list = (ls()))


## Observed 2010 deaths and exposures for municipios
big.df <- read.csv("C:/Users/Dell/OneDrive/Material GitHuh - Paper PHM/Input Data/new.big.df.csv", head = TRUE, stringsAsFactors = FALSE) %>% 
  select(-X)

## municipio-level geographic data (names and codes for state, meso, micro, etc.)
muni.url = 'http://schmert.net/topals-mortality/data/muni.df.csv'
muni.df  = read.csv(muni.url , header=TRUE, stringsAsFactors = FALSE)

df <- left_join(big.df, muni.df)

## Agregando por microrregião

big.df2 <- df %>%
  group_by(ufname,microcode,microname,sex,age) %>%
  summarise(death = sum(death), expos = sum(expos)) %>% 
  ungroup()


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
tmp = big.df2 %>%
  group_by(sex,age) %>%
  summarize(lambda = log( sum(death)/sum(expos)))

BR = matrix(tmp$lambda, ncol=2, dimnames=list(0:99, c('f','m')))

## smoothed version: project onto a cubic spline with closely-spaced knots
basis = bs(0:99, knots=seq(0,99,2))
Proj  = basis %*% solve(crossprod(basis)) %*% t(basis)
BRstd = Proj %*% BR

## check standard:
plot(BR[,1]) # 1 = female # 2 = male
lines(BRstd[,1])


# INÍCIO DO LOOPING: rodar para cada sexo ---------------------------------

# Parâmetors necessários:
age = 0:99
B   = bs( 0:99, knots=c(0,1,10,20,40,70), degree=1 )
micro.list = unique(muni.df$microcode)
sex.id = c('f', 'm')

BRstd = as.data.frame(BRstd)

## Por sexo: ---------------------------------------------------------------------

# Homens

topals.fit = as.data.frame(matrix(NA, 0, 10)) #data frame para guardar os resultados do looping

system.time({
  for(this.micro in micro.list){
    
    tmp = big.df2 %>%
      filter(sex=='m', microcode==this.micro) %>% #mudar o sexo
      group_by(ufname, microcode, microname, sex, age) %>%
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
    
    print(this.micro)
    
  }
})

colnames(topals.fit) <- c("ufname","microcode","microname",
                          "sex","age","death","expos",    
                          "logmx.inf","logmx.med", "logmx.sup")
topals.fit.micro.male <- topals.fit
write.csv(topals.fit.micro.male,"topals.fit.micro.male.csv")
# write_rds(topals.fit.male,"topals.fit.male.rds")



# Mulheres

topals.fit = as.data.frame(matrix(NA, 0, 10)) #data frame para guardar os resultados do looping

system.time({
  for(this.micro in micro.list){
    
    tmp = big.df2 %>%
      filter(sex=='f', microcode==this.micro) %>% #mudar o sexo
      group_by(ufname, microcode, microname, sex, age) %>%
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
    
    print(this.micro)
    
  }
})

colnames(topals.fit) <- c("ufname","microcode","microname",
                          "sex","age","death","expos",    
                          "logmx.inf","logmx.med", "logmx.sup")
topals.fit.micro.female <- topals.fit
write.csv(topals.fit.micro.female,"topals.fit.micro.female.csv")
# write_rds(topals.fit.male,"topals.fit.male.rds")

#Topals ambos os sexos - microrregião

topals.fit.micro.both <- rbind(topals.fit.micro.male, topals.fit.micro.female)

write.csv(topals.fit.micro.both,"topals.fit.micro.2sexos.csv")
# write_rds(topals.fit.both,"topals.fit.2sexos.rds")

# Passo 2: Juntando as bases do Topals para micro e modelo bayesiano para micro ------------

topals.fit.micro.both = read.csv("topals.fit.micro.2sexos.csv", sep = ',')
bayesiano = read.csv("logmx.micro.2010.csv", sep = ",")
bayesiano = bayesiano %>% 
  select(microcode,sex,age,q50)
gc_micro = left_join(topals.fit.micro.both,bayesiano)



# Passo 3: Calculando o grau de cobertura em cada idade por sexo e microrregião ------------

gc_micro <- gc_micro %>%
  mutate(grau_cobertura= exp(logmx.med)/exp(q50))
head(gc_micro)

write.csv(gc_micro, "graus.cobertura.micro.csv")

# Passo 4: Cálculo do grau de cobertura médio ponderado pela proporção de óbitos em cada idade por sexo e micro -----------------------------------
# OBS.: Apenas para construção do mapa

proporcao <- gc_micro %>%
  group_by(microcode,sex) %>%
  mutate(prop = death/sum(death))

gc_ponderado <- proporcao %>%
  group_by(microname, microcode, sex) %>%
  summarise(gc_pond = sum(grau_cobertura * prop))


# Mapas -------------------------------------------------------------------

install.packages("brazilmaps")
library(brazilmaps)
library(sf)
library(ggplot2)

# base de dados dos mapas
mapa_micro <- get_brmap("MicroRegion")

#Juntando as bases do estudo com as bases dos mapas

micro_graudecobertura <- mapa_micro %>%
  inner_join(gc_ponderado, c("MicroRegion" ="microcode"))
glimpse(micro_graudecobertura)

# Construção do mapa

map1 = micro_graudecobertura %>%
  mutate(gc_pond = cut(gc_pond,breaks =  c(quantile(gc_pond, probs = seq(0, 1, by = 0.20), na.rm = TRUE)),
                       include.lowest=TRUE, right = FALSE)) %>%
  ggplot() +
  geom_sf(aes(fill = gc_pond),
          # ajusta tamanho das linhas
          colour = "black", size = 0.05) +
  geom_sf(data = get_brmap("State"),
          fill = "transparent",
          colour = "white", size = 0.05) +
  # muda escala de cores
  scale_fill_viridis_d(option = 2, begin = 0.2, end = 0.8, name="Grau de cobertura\nPercentis 20,40,60,80,100", direction = -1) +
  #Título
  ggtitle("Grau de cobertura nas microrregi?es brasileiras por sexo, 2010") +
  # tira sistema cartesiano
  theme(panel.grid = element_line(colour = "transparent"),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  facet_wrap(vars(sex))

ggsave("graus_cobertura_micro_novo.png", map1, width=8, height=6, dpi=300)
