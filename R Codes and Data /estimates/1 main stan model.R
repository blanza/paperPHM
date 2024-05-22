################################################################
# Created for Carl 
#
#-----------------------
# Housekeeping
#----------------------
rm(list=ls())   
graphics.off()
if (.Platform$OS.type == 'windows') windows(width=8, height=8, record=TRUE)

set.seed(6447100)

library(withr)
library(dplyr)
library(splines)
library(rstan)


#-----------------------
# Job Parameters
#----------------------
save.fit  = TRUE
nchains   = 4
niter     = 1200
warmup    = 200

#------------------------------------
# Input Data and other basic info
#------------------------------------

## Observed 2010 deaths and exposures for municipios
big.df <- read.csv("C:/Users/Dell/OneDrive/Material GitHuh - Paper PHM/Input Data/new.big.df.csv", head = TRUE, stringsAsFactors = FALSE) %>% 
  select(-X)

str(big.df)

## municipio-level geographic data (names and codes for state, meso, micro, etc.)
muni.url = 'http://schmert.net/topals-mortality/data/muni.df.csv'
muni.df  = read.csv(muni.url , header=TRUE, stringsAsFactors = FALSE)

##  Death Distribution estimates of STATE-level pi[30+]
coverage.df = read.csv('C:/Users/Dell/OneDrive/Material GitHuh - Paper PHM/Input Data/estimated-coverage-by-UF2010.csv', 
                       skip=7, header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(a=NA, b=NA)

## calculate Method of Moments estimates of beta distribution parameters (a,b)
## for each set of coverage estimates in coverage.df

for (i in 1:nrow(coverage.df)) {
  est = coverage.df[i,] %>%
    select(CNPq.Project:Queiroz.seg.adj.2017) %>%
    as.numeric()
  
  xbar = mean(est)
  s2   = var(est)
  
  K = xbar*(1-xbar)/s2 - 1
  coverage.df$a[i] = xbar * K
  coverage.df$b[i] = (1-xbar)*K
}


## read the busca ativa estimates (ba = busca ativa estimated deaths, rep = reported deaths,
##   .inf is infants,  .tot is total at all ages)
busca.ativa.df = read.csv('C:/Users/Dell/OneDrive/Material GitHuh - Paper PHM/Input Data/busca-ativa-estimates-and-reported-deaths.csv',
                          header=TRUE, stringsAsFactors = FALSE) %>%
  left_join(select(muni.df, munid,ufabb,ufcode,contains('code')), by='munid')


ufcodes = c(ro=11 , ac=12 , am=13 , rr=14 , pa=15 , ap=16 , to=17 , ma=21 , 
            pi=22 , ce=23 , rn=24 , pb=25 , pe=26 , al=27 , se=28 , ba=29 , 
            mg=31 , es=32 , rj=33 , sp=35 , pr=41 , sc=42 , rs=43 , ms=50 , 
            mt=51 , go=52 , df=53)

ufnames = as.character(c('Rondonia', 'Acre','Amazonas','Roraima','Para', 'Amapa', 
                         'Tocantins', 'Maranhao','Piaui', 'Ceara', 'Rio Grande do Norte', 
                         'Paraiba','Pernambuco', 'Alagoas', 'Sergipe', 'Bahia',
                         'Minas Gerais','Espirito Santo','Rio de Janeiro','Sao Paulo',
                         'Parana','Santa Catarina', 'Rio Grande do Sul','Mato Grosso do Sul',
                         'Mato Grosso', 'Goias','Distrito Federal'))

ufabbs = toupper(names(ufcodes)) 

uf.df = data.frame( ufabb=ufabbs, 
                    ufcode=ufcodes, 
                    ufname=ufnames, 
                    stringsAsFactors = FALSE) %>%
  arrange(ufcodes)                        # sort for later merging


## calculate a "Brazil standard" by sex: result is a 100 x 2 matrix with
##   female and male stds in 'f' and 'm' columns, respectively
tmp = big.df %>%
  group_by(sex,age) %>%
  summarize(lambda = log( sum(death)/sum(expos)))

BRstd = matrix(tmp$lambda, ncol=2, dimnames=list(0:99, c('f','m')))

## smoothed version: project onto a cubic spline with closely-spaced knots
basis = bs(0:99, knots=seq(0,99,2))
Proj  = basis %*% solve(crossprod(basis)) %*% t(basis)

BRstd = Proj %*% BRstd


## proportions of reported deaths, by sex. in age groups 0, 1-29, 30+
frac.deaths           = cbind(m=c(.035, .109, .856), 
                              f=c(.037, .047, .916))
rownames(frac.deaths) = c('0', '1-29', '30+')

### build linear B-spline matrix (100 x 7)
age = 0:99
B   = bs( age, knots=c(0,1,10,20,40,70), degree=1)

nalpha = ncol(B)

#######################################################
# STAN MODEL WITH INCOMPLETE REGISTRATION
#######################################################

model.name='incomplete registration'

stanModelText = "
data {
  int<lower=0> R;       // number of regions
  int<lower=0> S;       // number of states (either 1 or 27)
  int<lower=0> A;       // number of ages (usually 100)
  int<lower=0> K;       // number of basis functions (usually 7)
  
  matrix[A,K] B;             // spline basis functions (usually 100 x 7)
  vector[A] std_schedule;    // std logmu schedule for TOPALS

  vector[R] agg_wts;       // wts for aggregating pi[3]s up to state level (won't be used if state is the 'small' geog)
  vector[3] death_frac;    // fraction of deaths at ages 0, 1-29, 30+
  
  vector<lower=0,upper=1>[R] pihat1;    //  pi[r,1] ~ beta( K1*pihat1[r], K1*(1-pihat1[r]))
  vector<lower=0,upper=1>[R] pihat_all;  //  pi[r,ALL] ~ beta( Kall*pihat_all[r], Kall*(1-pihat_all[r]))
  
  vector[S] a3 ;                      //  pibar[3] ~ beta(a3,b3)  state agg
  vector[S] b3 ;          
  
  matrix<lower=0>[A,R] N      ;  // age- and region-specific exposure
  int<lower=0>         D[A,R] ;  // observed deaths by age and region

}

transformed data {
  matrix[A,R] lambda_star;   // std logmu schedule for TOPALS in each column
  
  lambda_star = rep_matrix(std_schedule, R);
}

parameters {
  matrix[K,R] alpha;   // TOPALS offsets
  
  simplex[4] eta[R]; // used to construct (pi1, pi2, pi3)
  
  real<lower=0> K1;    // precision (hyper)parameter for regional pi1 
  real<lower=0> Kall;  // precision (hyper)parameter for regional piall
}

transformed parameters {
  matrix<lower=0,upper=1>[R,3] pi; 
  vector<lower=0,upper=1>[R]   piall;    // coverage probabilities at all ages combined
  real<lower=0, upper=1>       agg_pi3;  // state-level coverage for deaths at 30+
  
  for (region in 1:R) {
  
    pi[region,1] = eta[region,1];
    pi[region,2] = eta[region,1] + eta[region,2] + eta[region,3];  
    pi[region,3] = eta[region,1] + eta[region,2];    
  
  } # for region
  
  piall = pi * death_frac;   // result is an R-vector of avg covg (all ages) by region
  
  agg_pi3 = dot_product( pi[,3], agg_wts); 
  
}

model {
  matrix[A,R]    lambda;               // log mx rates
  matrix[K-1,R]  alpha_diff;           // differences in offsets over adjacent groups
  
  lambda      = lambda_star + B * alpha ;   // regional logmu schedules ages 0,1,..99

  alpha_diff  = alpha[2:K,] - alpha[1:(K-1),];   // first diffs in alphas
  
  ## LIKELIHOOD: Poisson distribution of deaths    
  ## log(N) + lambda + log(pi) = ln(expect death reports)       
  
  for (region in 1:R) {
    D[1,region]    ~ poisson_log( log(N[1,region])    + lambda[1,region]     + log(pi[region,1])) ;  
    D[2:30,region] ~ poisson_log( log(N[2:30,region]) + lambda[2:30,region]  + log(pi[region,2])) ;  
    D[31:A,region] ~ poisson_log( log(N[31:A,region]) + lambda[31:A, region] + log(pi[region,3])) ;  
  }
  
  ## PRIORS
  to_vector(alpha)       ~ normal(0,4);
  to_vector(alpha_diff)  ~ normal(0, 1/sqrt(2));

  K1      ~ exponential(0.05) ;   // truncated exponential for precision of K1
  Kall    ~ exponential(0.05) ;   // truncated exponential for precision of Kall
  
  pi[,1] ~ beta(     (5+K1)*pihat1, (5+K1)*(1-pihat1) );
  
  piall  ~ beta( (5+Kall)*pihat_all, (5+Kall)*(1-pihat_all)  ) ;    // vectorized: prior for total covg regions 1..R
  
  if (S  > 1) pi[,3]  ~ beta(a3, b3 );    // prior for state-level pi3 values (if state is the SMALL area)
  if (S == 1) agg_pi3 ~ beta( a3, b3) ;   // prior for state-level pi3 value  (if state is the LARGE area)
  
  

} //model
"
#######################################################
# STAN MODEL WITH 100% REGISTRATION (all pi=1)
#######################################################

stanModelText_complete_reg = "
data {
  int<lower=0> R;       // number of regions
  int<lower=0> S;       // number of states (either 1 or 27)
  int<lower=0> A;       // number of ages (usually 100)
  int<lower=0> K;       // number of basis functions (usually 7)
  
  matrix[A,K] B;             // spline basis functions (usually 100 x 7)
  vector[A] std_schedule;    // std logmu schedule for TOPALS
  
  matrix<lower=0>[A,R] N      ;  // age- and region-specific exposure
  int<lower=0>         D[A,R] ;  // observed deaths by age and region
}

transformed data {
  matrix[A,R] lambda_star;   // std logmu schedule for TOPALS in each column

  lambda_star = rep_matrix(std_schedule, R);
}

parameters {
  matrix[K,R] alpha;   // TOPALS offsets
}

model {
  matrix[A,R]   lambda;               // log mx rates
  matrix[K-1,R] alpha_diff;           // differences in offsets over adjacent groups
  
  lambda      = lambda_star + B * alpha ;        // regional logmu schedules ages 0,1,..99

  alpha_diff  = alpha[2:K,] - alpha[1:(K-1),];   // first diffs in alphas
  
  ## LIKELIHOOD: Poisson distribution of deaths    
  ## log(N) + lambda + log(pi) = ln(expect death reports)       
  
  for (region in 1:R) {
    D[,region]    ~ poisson_log( log(N[,region]) + lambda[,region]  ) ;  
  }
  
  ## PRIORS
  ## very weak prior on alpha differences (stabilizes estimates when exposure is very small)

  to_vector(alpha)       ~ normal(0,4);
  to_vector(alpha_diff)  ~ normal(0, 1/sqrt(2));

} //model
"

## compile the models

part_reg_model  = stan_model( model_code = stanModelText)
full_reg_model  = stan_model( model_code = stanModelText_complete_reg)


#==============================================================================
# MAIN LOOP OVER MANY CASES
#==============================================================================


## generate the cases for each state/sex combination at microregion level
case = expand.grid(sex=c('f','m'), 
                   big_area = ufabbs,  
                   small_area='micro'  ) %>%
  filter(big_area != 'DF')

## append the cases for all of Brazil, with states as small areas
case = rbind(case,
             expand.grid(sex=c('f','m'), 
                         big_area='Brazil',
                         small_area='state'))


####################################

time.stamp = format(Sys.time(), '%d%b%y-%H%M')

for (casenum in 1:nrow(case)) {
  
  print(case[casenum,])
  
  this.sex = as.character( case$sex[casenum])
  
  geog   = list(big   = case$big_area[casenum], 
                small = case$small_area[casenum])
  
  geog.label = paste0('-',geog$small,'-within-',geog$big)
  sex.label  = c(m='-Male',f='-Female')[this.sex]
  
  fitfile    = paste0('Stan-',time.stamp, geog.label, sex.label,'.RData')
  
  ## aggregate busca ativa estimates to the chosen 'small' geog
  if (geog$small == 'state') busca.ativa.df$selected.code = busca.ativa.df$ufcode
  if (geog$small == 'meso')  busca.ativa.df$selected.code = busca.ativa.df$mesocode
  if (geog$small == 'micro') busca.ativa.df$selected.code = busca.ativa.df$microcode
  
  
  field_audit.df =  busca.ativa.df %>%
    group_by(ufcode,ufabb,selected.code) %>%
    summarize( pi1.hat   = max( .005, min(.995, sum(rep.inf)/ sum(ba.inf))),
               piall.hat = max( .005, min(.995, sum(rep.tot)/ sum(ba.tot)))) %>%
    ungroup
  
  if (geog$big %in% ufabbs) { field_audit.df = filter(field_audit.df, ufabb == geog$big) }
  
  ################################################################
  # construct the data and priors for small areas
  ################################################################
  
  ## DATA 
  ## step 1. read death and exposure data, keep selected sex and all municipios in selected big area
  
  tmp = left_join( big.df, muni.df, by='munid') %>%
    filter(sex==this.sex)
  
  if (geog$big %in% ufabbs) tmp = filter(tmp, ufabb == geog$big)
  
  ## step 2. aggregate the remaining data into the selected small areas
  if (geog$small == 'state') tmp$selected.code = tmp$ufcode
  if (geog$small == 'meso')  tmp$selected.code = tmp$mesocode
  if (geog$small == 'micro') tmp$selected.code = tmp$microcode
  
  tmp = tmp %>%
    group_by(selected.code, ufabb, ufcode, age) %>%
    summarize(D = sum(death), N=sum(expos)) 
  
  D = matrix(tmp$D, nrow=100)   # a 100 ages x R regions matrix of death counts
  N = matrix(tmp$N, nrow=100)   # a 100 ages xR regions matrix of exposure 
  
  colnames(D) = colnames(N) = unique(tmp$selected.code)
  
  ## small region names for later labeling
  if (geog$small == 'state') region.names = ufabbs
  if (geog$small != 'state') region.names = colnames(N)
  
  
  ## PRIOR INFO
  
  ## beta parameters for coverage at 30+, state level
  ab = filter(coverage.df, sex==this.sex) %>%
    select(a:b) %>%
    as.matrix()
  
  rownames(ab) = ufcodes
  
  ## if using all states, keep all coverage estimates. Otherwise only keep the estimates for state that is the big area
  if (geog$big == 'Brazil') {
    prior.a.pi3bar = ab[,1]
    prior.b.pi3bar = ab[,2]
  } else {
    prior.a.pi3bar = ab[which(ufabbs == geog$big),1]
    prior.b.pi3bar = ab[which(ufabbs == geog$big),2]
  }
  
  ## weighted point estimates for coverage at ages 0 
  ## and at all ages, from busca ativa model
  ## MAKE SURE THAT THESE ARE IN THE SAME ORDER AS THE ROWS of D and N
  
  o = match( colnames(D), field_audit.df$selected.code)
  
  prior.mean.pi1   = field_audit.df$pi1.hat[o]      # R x 1
  prior.mean.piall = field_audit.df$piall.hat[o]    # R x 1
  
  
  #######################################################
  # STAN Data
  #######################################################
  
  R=ncol(N)
  S=ifelse( geog$big=='Brazil', 27, 1)  # number of states
  A=nrow(N)
  K=ncol(B)
  
  stanDataList <- list(
    R=R,
    S=S,
    A=A,
    K=K,
    
    B            = B,
    std_schedule = BRstd[,this.sex],
    
    agg_wts = colSums( N[age>29,] ) / sum(N[age>29,]),
    death_frac = frac.deaths[,this.sex],
    
    pihat1    = prior.mean.pi1,    
    pihat_all = prior.mean.piall,
    a3    = array(prior.a.pi3bar, dim=S), #array() is necessary because Stan doesn't handle a 'vector' of length 1 easily
    b3    = array(prior.b.pi3bar, dim=S),
    
    N = pmax(N, .01),
    D = D
  )
  
  #######################################################
  # STAN initialization function
  #######################################################
  
  stanInits = function() {
    p1 = rbeta(R, 100*prior.mean.pi1, 100*(1-prior.mean.pi1))
    p3 = runif(R, p1, 1)
    p2 = runif(R, p3, 1)
    
    eta1 = p1
    eta2 = p3-p1
    eta3 = p2-p3
    eta4 = 1-p2
    
    a  = matrix(runif(K*R, -0.10, +0.10), K, R)
    
    list(
      alpha  = a,
      eta  = cbind(eta1,eta2,eta3,eta4),
      
      K1    = rexp(1,.05),
      Kall  = rexp(1,.05)
    )
  } # stanInits
  
  
  #######################################################
  # SAMPLE FROM POSTERIOR FOR FULL MODEL, INCL
  # REGISTRATION PROB <1
  #######################################################
  
  fit <- sampling(part_reg_model,
                  data    = stanDataList,
                  init    = stanInits,
                  seed    = 6447100,
                  pars    = c('alpha','pi','piall'),
                  control = list(max_treedepth=10),
                  iter    = niter, 
                  warmup  = warmup,
                  thin    = 1, 
                  chains  = nchains)
  
  #######################################################
  # APPROXIMATE POSTERIOR FOR ABRIDGED MODEL WITH PI=1
  # (several tests showed that these results are essentially
  # identical to MCMC sampling from posterior for the pi=1
  # model)
  #######################################################
  ML = optimizing(full_reg_model, 
                  data  = stanDataList, 
                  init  = stanInits, 
                  seed  = 6447100,
                  draws = 2000)
  
  if (save.fit) save( fit, stanDataList, ML, region.names, file=fitfile)
  
} # for casenum
