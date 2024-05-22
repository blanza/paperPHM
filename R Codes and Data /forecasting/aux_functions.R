####################################
######## FUNCOES AUXILIARES ########
####################################

### Carrega/instala pacotes faltantes #------
open_install_package <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}
##############################################

### Tabua de vida #---------------------------
life.table <- function( x, nMx ){
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
  lt <- data.table(age = x, nMx = nMx, nqx = nqx, lx = lx, ndx = ndx, nLx = nLx, Tx = Tx, ex = ex)
  return(lt)
}
##############################################

### Bayes Empirico #--------------------------

tx.bemp<-function(obs,ni){


# (obs=vetor dos eventos por peque area i;ni= vetor das pop nas peq areas i)

  n<-sum(ni)
  media<-sum(obs)/n
  xi<-obs/ni
  areas<-length(xi)
  s2<- (sum(ni*((xi-media)^2)))/n
  nmed<-n/areas
  A<-s2-(media/nmed)
  ci<-NULL
  for(i in 1:areas){

    aux1<- A/(A+(media/ni[i]))
    if (A>=0) {ci[i]=aux1}
    else {ci[i]=0}
  }
  ci[is.na(ci)]=0

  xi.bayes<- media+(ci*(xi-media))

  result<-list(taxa.bayesiana=xi.bayes,taxa.original=xi,ci=ci)
  return(result)
}


tx.bemp.modificada<-function(obs,ni){


  # (obs=vetor dos eventos por peque area i;ni= vetor das pop nas peq areas i)

  n<-sum(ni)
  media<-sum(obs)/n
  xi<-obs/ni
  areas<-length(xi)
  s2<- (sum(ni*((xi-media)^2)))/n
  nmed<-n/areas
  A<-s2-(media/nmed)
  ci<-NULL
  for(i in 1:areas){

    aux1<- A/(A+(media/ni[i]))
    if (A>=0) {ci[i]=aux1}
    else {ci[i]=0}
  }
  ci[is.na(ci)]=0

  xi.bayes<- media+(ci*(xi-media))

  result <-
    data.table(taxa.bayesiana=xi.bayes,taxa.original=xi,ci=ci)

  result[,cor.bayes:=taxa.bayesiana/taxa.original]
  return(result)
}
##############################################

### Lee-Carter #------------------------------


LeeCarterForecast <- function( ax,
                               kt,
                               bx,
                               variant_low  = FALSE,
                               variant_high = FALSE,
                               year_start,
                               year_end ){

  time_span <- year_end - year_start

  # loop:
  age <- c(0:99)

  # Preparando kt para projecao:
  kt_diff <- diff(kt)
  summary_kt <- summary( lm( kt_diff ~ 1 ) )
  kt_drift <- summary_kt$coefficients[1,1]
  sec <- summary_kt$coefficients[1,2]
  see <- summary_kt$sigma

  # Projetando kt usando um modelo ARIMA, ruido aleatorio com tendencia:
  h <- seq(0, time_span) # anos projetados (0 = 2010, 1 = 2011, 20 = 2030)
  kt_stderr <- ( ( h * see ^ 2 ) + ( h * sec ) ^ 2 ) ^ .5
  kt_forecast <- tail( kt, 1 ) + ( h * kt_drift )
  kt_forecast_lo <- kt_forecast - ( 1.96 * kt_stderr )
  kt_forecast_hi <- kt_forecast + ( 1.96 * kt_stderr )

  # normalizando um novo kt para que em 2010 seja zero
  kt_forecast_std <- as.matrix( ( kt_forecast - head( kt_forecast, 1 ) ) )

  if( variant_low == TRUE & variant_high == FALSE){
    kt_forecast_std <- kt_forecast_std - ( 1.96 * kt_stderr )
  }

  if( variant_high == TRUE & variant_low == FALSE ){
    kt_forecast_std <- kt_forecast_std + ( 1.96 * kt_stderr )
  }

  forecast_logmx <- matrix( nrow = length( kt_forecast_std ), ncol = length( age ) )

  for (j in 1 : length( kt_forecast_std ) ){
    forecast_logmx[j, ] <- ( bx * kt_forecast_std[ j ] ) + ax
  }

  colnames( forecast_logmx ) = paste0('logmx.', age)

  forecast_logmx_wide <-
    data.frame( year = year_start:year_end,
                forecast_logmx )

  forecast_logmx_long <-
    reshape( forecast_logmx_wide,
             varying = list(2:101),
             timevar = 'age',
             times = age,
             v.names = 'logmx',
             idvar = 'year',
             direction = 'long' )

  forecast_logmx_long <-
    forecast_logmx_long[ order( forecast_logmx_long[,1],
                                forecast_logmx_long[,2]), ]

  return(forecast_logmx_long)
}

##############################################
