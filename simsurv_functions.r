ibs_pec <- function(data, tte, cens, x){
  if(all(data[, cens] == 1)){
    print('pec IBS calculation seems wrong when the input data has no censor')
  }
  
  form <- as.formula(paste("Surv(",tte, ",", cens,")~", x))
  coxfit <- coxph(form,data=data, y=TRUE, x=TRUE)
  
  form_censor <- as.formula(paste("Surv(",tte,",",cens,")~1"))
  
  e1 <-  pec(object=list('MyModel' = coxfit), formula=form_censor, data=data, exact=TRUE, 
             cens.model="marginal",splitMethod="none",B=0, verbose=F)
  i <- unclass(crps(e1))
  c(original = i['MyModel',1], scaled  = 1 - i['MyModel',1] / i['Reference',1])
}

ibs_strata <- function(data, tte, cens, x, strat.var){
  form <- as.formula(paste("Surv(",tte, ",", cens,")~", x, '+ strata(',strat.var,')'))
  coxfit <- coxph(form,data=data, y=TRUE, x=TRUE)
  
  form_censor <- as.formula(paste("Surv(",tte,",",cens,")~ 1 + strata(",strat.var,')'))
  stratRef <- coxph(form_censor, data=data, y=TRUE, x=TRUE)
  
  e1 <-  pec(object=list('MyModel' = coxfit, 'stratRef' = stratRef), formula=form_censor, data=data, exact=TRUE, cens.model="nonpar",splitMethod="none",B=0, verbose=F)
  
  i <- unclass(crps(e1))
  c(original =   i['MyModel',1], scaled  = 1 - i['MyModel',1] / i['Reference',1], scaled2 =1 - i['MyModel',1] / i['stratRef',1] )
}


p1 <- c(gamma = 1, lambda =0.05597537,   trt = - 0.3)


weibull_surv <- function(t, lambda=1, gamma=1, linearPredictor){
  exp(-1 * lambda * (t^gamma) * exp(linearPredictor))
}

if(F){
  summary(weibull_surv(t=1,linearPredictor= rnorm(n=100, sd=3 )))
  summary(weibull_surv(t=1, linearPredictor= rnorm(n=100, sd=0.3 )))
  
  # survival probability increase when sd becomes smaller
}


#====decide follow up time==========

taylor_exp <- function(sd, a){
  # to roughly calculate the expected value of exp(a*x), where x is normal with mean 0 and sd
  # E(exp(a*x)) = 1 + ax + (ax)^2/2 + (ax)^3/3! + (ax)^4/4! + (ax)^5/5!
  set.seed(27519)
  1 + a^2 /2 * sd^2 + a^4/2/4 * sd^4  
  #https://math.stackexchange.com/questions/92648/calculation-of-the-n-th-central-moment-of-the-normal-distribution-mathcaln
}


# this function seems to only work with alpha = 0. it assumes the hazard at the landmark time point is the same across time point, which dramatically over-estimates the hazard at the beginning.
fupByEPS_weibull <- function(targetEPS=0.75, maxFollowUp=100,  para=c(), tol1=0.01){

  for (v in setdiff(names(para_template), names(para))) para[v] <- para_template[v]
  
  # only requires the active arm to reach the designed EPS 
  # use taylor exp
#   f <- function(t) log(1-targetEPS) * -1 - (t ^ para["gamma"] ) * taylor_exp(sd = para['sd_long'], a = para['alpha']) *   exp(para["betaEvent_intercept"] +  para["betaEvent_trt"] + para["alpha"] * para['bm_trt']) 

  # use simulation
  set.seed(27519)
  
   f <- function(t) mean( log(1-targetEPS) * -1 - (t ^ para["gamma"] ) * exp(para["betaEvent_intercept"] +  para["betaEvent_trt"] + para["alpha"] * ( para['bm_trt'] + rnorm(n=1e5, sd=para['sd_long']) )) )
   fup <- uniroot(f, interval=c(0,maxFollowUp), tol=tol1)$root    
    fup
}

#ftime <- fupByEPS_weibull (targetEPS=0.75, maxFollowUp=200,  para=c(alpha=5, sd_long=0.3))  


censorTimeByEPS <- function(ds, tte = 'eventtime', cens = 'status', targetEPS = 0.75){
  mean_eps <- mean(ds[, cens])
  if(   mean_eps >= targetEPS ) {
    maxt <- quantile(ds[, tte], probs=targetEPS)
    sel <- ds[, tte] >= maxt
    ds[sel, tte] <- maxt
    ds[sel, cens] <- 0
  }else{
    stop(paste('observed EPS is', round(mean_eps,2),'; no censoring is applied'))
  }
  
  ds
}


#=========apply random censor ============
survival_lambda_exponential <- function(targetRate = 0.075, t ){
  log(1- targetRate)/ -1 /t
}

#survival_lambda_exponential(targetRate = 0.25, t=30) 



#=========main functions




rfsrc_cindex <- function(...){
#  randomForestSRC:::cindex(...) # 2.8.0 or early
  randomForestSRC::get.cindex(...) # 2.9.0 or late
}


haz <- function(t, x, betas, k=2, returnTrajectoryUpto=0 , addRandomEffectTo = 'a', ...) {
  
  a <- ( betas[['betaLong_trt']] * x[['trt']] +   betas[['betaLong_soc']] * (1 - x[['trt']]) ) * -1
  
  stopifnot(!is.null(betas[['longRandomIntercept']]))
  
  if (addRandomEffectTo == 'a') a <- a + betas[['longRandomIntercept']]
  # this way, betaLong_trt = -1 will make a = 1 for the treatment arm
  
  if(returnTrajectoryUpto > 0){
    s1 <- sapply(seq(0,returnTrajectoryUpto, length.out=10), function(t1) {
      (2/(1 + exp(k * t1)) - 1) *  a 
    })
    
    colnames(s1) <- seq(0,returnTrajectoryUpto, length.out=10)
    
    s1 + betas[['longRandomIntercept']] * ifelse(addRandomEffectTo == 'y', 1, 0)
    
  }else{
    
    for (v in c("gamma", "betaEvent_intercept", "betaEvent_binary", "betaEvent_assoc")) {
      stopifnot(!is.null(betas[[v]]))
     # print( paste(v, round(unique(betas[[v]]),3), sep=' : '))
    }
    
    betas[["gamma"]] * (t ^ (betas[["gamma"]] - 1)) * exp(
      betas[["betaEvent_intercept"]] +
        betas[["betaEvent_binary"]] * x[["trt"]] +
        betas[["betaEvent_assoc"]] * (
          (2/(1 + exp(k * t)) - 1) * a + betas[['longRandomIntercept']] * ifelse(addRandomEffectTo == 'y', 1, 0)
        )
    )
    
  }
}

# bi-exponential trajector
haz_biexp <- function(t, x, betas, returnTrajectoryUpto=0 , ...) {
  
  ks <- ( betas[['ks_trt']] * x[['trt']] +   betas[['ks_soc']] * (1 - x[['trt']]) ) 
  kg <- ( betas[['kg_trt']] * x[['trt']] +   betas[['kg_soc']] * (1 - x[['trt']]) ) 
  

  if(returnTrajectoryUpto > 0){
    s1 <- sapply(seq(0,returnTrajectoryUpto, length.out=10), function(t1) {
      exp(-1* ks * t1) + exp(kg * t1) - 2
    })
    
    colnames(s1) <- seq(0,returnTrajectoryUpto, length.out=10)
    s1
  }else{
    
    for (v in c("gamma", "betaEvent_intercept", "betaEvent_binary", "betaEvent_assoc")) {
      stopifnot(!is.null(betas[[v]]))
      # print( paste(v, round(unique(betas[[v]]),3), sep=' : '))
    }
    
    betas[["gamma"]] * (t ^ (betas[["gamma"]] - 1)) * exp(
      betas[["betaEvent_intercept"]] +
        betas[["betaEvent_binary"]] * x[["trt"]] +
        betas[["betaEvent_assoc"]] * (
          exp(-1* ks * t) + exp(kg * t) - 2
        )
    )
    
  }
}


# simplied hazard
haz1 <- function(t, para1, k=2,  ...) {
  
  para1["gamma"] * (t ^ (para1["gamma"] - 1)) * exp(
    para1[["betaEvent_intercept"]] +   para1["alpha"] * (
      (2/(1 + exp(k * t)) - 1) * para1['bm_trt'] * -1  
    )
  )
  
}




if(F){
  lpRange <- function(lm=2, maxt =180, para=c()){
    #print(seed)
    for (v in setdiff(names(para_template), names(para))) para[v] <- para_template[v]
    
    cumHaz <- list( 
      integrate(haz1, lower=0, upper=lm,   para1=para),
      integrate(haz1, lower=0, upper=maxt, para1=para) )
    1 - exp(-1 * sapply(cumHaz, function(x)x$value))
  }
  
  lpRange(para=c())
  
  pList <- expand.grid(bm_trt=seq(-0.8,0, 0.2),  betaEvent_intercept= -4:0, alpha = 0:6)
  eventRate <- apply(pList, 1, function(p1) {
    
    lpRange(para=c(alpha=unname(p1['alpha']), betaEvent_intercept = unname(p1['betaEvent_intercept']), bm_trt=unname(p1['bm_trt']))) })
  
  pList2 <- cbind(pList, t(eventRate))
  
  write.csv(pList2, file = '~/aPDL1/oak/results/linearPredictor_range.csv', row.names=F)
}

# maxt is necessary because simsurv may complain that it could not find a reasonable time for some patients when the hazard is really low
sim1 <- function(lm=2, seed = 27519, N = 100, maxt =180, targetEPS=0.75, uniformCensor=F, operationSummaryOnly=F, para=c()){
  #print(seed)
  for (v in setdiff(names(para_template), names(para))) para[v] <- para_template[v]
  
  set.seed(seed)
  
  betas1  = data.frame(
    # gamma = rep(p1['gamma'], N),
    gamma = rep(para['gamma'], N),
    betaEvent_intercept = rep(para['betaEvent_intercept'], N),
    betaEvent_binary = rep(para['betaEvent_trt'], N),
    #betaEvent_assoc = rep(para['alpha'], N),
    betaLong_trt = rep(para['bm_trt'], N),
    betaLong_soc = rep(para['bm_ctrl'], N) 
  )
  
  covdat <- data.frame(
    trt = c(rep(1, N/2), rep(0, N/2))
  )

  if('alpha' %in% names(para)){
    betas1$betaEvent_assoc = rep(para['alpha'], N)
  }else{
    stopifnot(all(c('alpha1', 'alpha0') %in% names(para)))
    betas1$betaEvent_assoc = ifelse(covdat$trt == 1, 
                          rep(para['alpha1'], N), rep(para['alpha0'], N))
  }
    
  betas1[['longRandomIntercept']] <- rnorm(n= N, mean=0, sd=para['sd_long'])
  
  s2 <- simsurv(hazard = haz, x = covdat, betas = betas1, maxt = maxt)

  # apply practicall follow-up time
  if(!is.null(targetEPS)) s2 <- censorTimeByEPS(s2,  targetEPS=targetEPS)
  
  censor_maxfup <-  tapply(1 - s2[,'status'], covdat[,'trt'], sum)
  
  # TODO add censoring, which may affect atrisk status, event time, and event label
  dropout_censor <- 0
  readout_time <- max(s2[,'eventtime'])
  
  if(uniformCensor ){
    dropout_time <- runif(n =N, min = 0, max=readout_time)
    dropout_censor <- dropout_time < s2[,'eventtime'] 
    s2[dropout_censor, 'eventtime'] <- dropout_time[dropout_censor]
    
    dropout_censor <- dropout_censor & s2[, 'status']  == 1
    s2[dropout_censor, 'status'] <- 0
    
  }
  
  if( F && dropRateAtReadout > 0){
 # need update   
    dropout_lambda <- survival_lambda_exponential(targetRate = dropRateAtReadout, t = readout_time)
    dropout_time <- rexp(n=N, rate =dropout_lambda)  
     dropout_censor <- dropout_time < s2[,'eventtime']
    s2[dropout_censor, 'eventtime'] <- dropout_time[dropout_censor]
    s2[dropout_censor, 'status'] <- 0
  }
  
  
  
  # biomarker value, only use SE at the landmark visit
  tr1 <- haz(x = covdat, betas = betas1, returnTrajectoryUpto = lm)
  bm <- tr1[, ncol(tr1)]   # the last column is the biomarker value at the landmark time
  if('measure_sd' %in% names(para) && !is.na(para['measure_sd']) && para['measure_sd'] > 0) bm <- bm + rnorm(n=length(bm), mean=0, sd = para['measure_sd'])
  
  
  ds1 <- cbind(cbind(s2, covdat) , bm= bm )
  ds1$t_lm <- ds1$eventtime - lm
  ds1$atrisk = with(ds1, eventtime > lm)
  
  # summarize censor / early events by arm
 
  cen1 <- with(ds1[ds1$atrisk, ], tapply(1 - status, trt, sum))
#  cen1 <- with(ds1[ds1$atrisk, ], sum(1 - status))
  early_event <- with(ds1, tapply(1 - atrisk, trt, sum))
  
  # alternative ways to code biomarker 
  ds1$bm2 <- as.integer(ds1$bm >= median(ds1$bm))
  
  # impute biomarker values for patients who are not at risk at the landmark time (even though we have the simulated values
  #  ds1$bm_impute <- with(ds1, ifelse(atrisk, bm, max(bm)))
  
  # ingore (not impute) biomarker values that may be missing for patients who are not at risk
  bm_diff <- with(ds1, tapply(bm, trt, mean)) 
  
  bm_t <- t.test(bm ~ trt, data=ds1)$statistic
  
  k1 <- summary(survfit(Surv(eventtime, status) ~ trt, ds1))$table[,'median']
  
  if(operationSummaryOnly){
    return(# about trial operation
      c(seed = seed,
      N_overall_censor_atrisk = cen1, n_early_event = early_event, max_fup = readout_time, 
      N_censor_dropout = sum(dropout_censor), N_censor_maxfup = censor_maxfup, medianTime=k1
    ))
  }
  
  
  c1 <- coxph(Surv(eventtime, status) ~ trt , ds1)
  
  #  c2 <- coxph(Surv(t_lm, status) ~ trt , ds1[ ds1$atrisk, ])
  
  c3 <- coxph(Surv(t_lm, status) ~ bm + trt , ds1[ ds1$atrisk, ])
  c3a <- coxph(Surv(t_lm, status) ~ bm + strata(trt) , ds1[ ds1$atrisk, ])
  c3b <- coxph(Surv(t_lm, status) ~ bm , ds1[ ds1$atrisk, ])
  c3c <- coxph(Surv(t_lm, status) ~ bm , ds1[ ds1$atrisk & ds1$trt==1, ])
  c3d <- coxph(Surv(t_lm, status) ~ bm, ds1[ ds1$atrisk & ds1$trt==0, ])
  
  c4a <- coxph(Surv(t_lm, status) ~ bm2 + strata(trt) , ds1[ ds1$atrisk, ])
  # extract c index, harrell's
  
  #c_h <- concordance(object=c3a)$concordance
  c_rfs <- with(ds1[ ds1$atrisk, ], rfsrc_cindex(t_lm, status, bm * -1))
  c_rfs1 <- with(ds1[ ds1$atrisk & ds1$trt == 1, ], rfsrc_cindex(t_lm, status, bm* -1))
  c_rfs0 <- with(ds1[ ds1$atrisk & ds1$trt == 0, ], rfsrc_cindex(t_lm, status, bm* -1))
  
  
  # extract uno's c index
  
  unoc = with(ds1[ ds1$atrisk, ], UnoC(Surv.rsp = Surv(t_lm, status) , Surv.rsp.new = Surv(t_lm, status) , lpnew =bm ))

  unoc1 = with(ds1[ ds1$atrisk & ds1$trt == 1, ], UnoC(Surv.rsp = Surv(t_lm, status) , Surv.rsp.new = Surv(t_lm, status) , lpnew =bm ))
  unoc0 = with(ds1[ ds1$atrisk & ds1$trt == 0, ], UnoC(Surv.rsp = Surv(t_lm, status) , Surv.rsp.new = Surv(t_lm, status) , lpnew =bm ))
  
  # get IBS from pec
  
  ibs = ibs_pec(data =ds1[ ds1$atrisk, ] , tte = 't_lm', cens = 'status', x = 'bm')
  ibs1 = ibs_pec(data =ds1[ ds1$atrisk  & ds1$trt == 1, ] , tte = 't_lm', cens = 'status', x = 'bm')
  ibs0 = ibs_pec(data =ds1[ ds1$atrisk  & ds1$trt == 0, ] , tte = 't_lm', cens = 'status', x = 'bm')
#  ibs_strata = ibs_strata(data =ds1[ ds1$atrisk, ] , tte = 't_lm', cens = 'status', x = 'bm', strat.var = 'trt')
  
  c(# about trial operation
    seed = seed,
    N_overall_censor_atrisk = cen1, n_early_event = early_event, max_fup = readout_time, 
    N_censor_dropout = sum(dropout_censor), N_censor_maxfup = censor_maxfup, medianTime=k1,
    
    # about performance measurement    
     # c3$coefficients['bm2'], c3a$coefficients['bm2'],c3b$coefficients['bm2'],c3c$coefficients['bm2'],c3d$coefficients['bm2'],  # binary biomarker
    

     #bmHR_adjusted_by_trt = c3$coefficients['bm'], 
     bmHR_stratifi_by_trt = c3a$coefficients['bm'],
#   bmHR_ingore_trt= c3b$coefficients['bm'],
     bmHR_trt1 = c3c$coefficients['bm'],
     bmHR_trt0 = c3d$coefficients['bm'], 

    bm2HR_stratified = c4a$coefficients['bm2'],

    #c_h_strata = c_h, 
    c_rfs = c_rfs, c_rfs_trt1 = c_rfs1, c_rfs_trt0 = c_rfs0,
    unoC=unoc, unoC_trt1 = unoc1, unoC_trt0 = unoc0,
    ibs = ibs,   ibs_trt1 = ibs1, ibs_trt0= ibs0,
    #ibs_strata =ibs_strata 

    # to support meta-regression
     'trt_itt' = c1$coefficients['trt'],    
     #   'trt_atrisk' = c2$coefficients['trt'],
     'bm_diff' = bm_diff['1'] - bm_diff['0'], 
     'bm_t' = bm_t
)
}

norm1 <- function(n, mean, sd){
  if(is.null(sd) || is.na(sd) || sd==0 ) {
    0
  }else{  
    rnorm(n=n, mean=mean, sd=sd)
  }
}

#' simulate tte studies under a biexponential model
#'
#' @param lm landmark time point 
#' @param seed 
#' @param N number patients from both arms (1:1 randomized to 2 arms)
#' @param maxt max follow up time for simulation using \code{simsurv}
#' @param targetEPS target event patient ratio in 2 arms combined. If the simulated study by \code{simsurv} has already a lower EPR, the function will throw an error. Otherwise, some long event times are censored so that the EPR is close to \code{targetEPS}. 
#' @param uniformCensor apply a uniform dropout between start and the time when the study reaches \code{targetEPS}
#' @param operationSummaryOnly only summarise operation characters in the returned vector
#' @param para parameters for simulation. Parameters specified here will overwrite those in para_template_biexp from the global environment
#' @param time_multiplier TGI model parameters use a unit of week, which are multiplied by this number so the time scale become month
#' @param reportPDchangeOnly default is F:  analyze TTE endpoints in both patient level performance metric and trial level association. If T:  do not analyze survival endpoints
#' @param imputeCensoredSE default is T: 1) if a patient is no longer at risk at \code{lm}, the biomarker evalue is imputed as the largest observed values (from both arms); 2) use wilcox test to compare such biomarker. If F: 1) assume we can measure the biomarker at \code{lm} time point even if the patient has an event before; 2) use t test
#'
#' @return
#' @export
#'
#' @examples
sim2_be <- function(lm=2, seed = 27519, N = 100, maxt , targetEPS=0.75, uniformCensor=F, operationSummaryOnly=F, para=c(), time_multiplier = 4, reportPDchangeOnly=F, imputeCensoredSE=T, returnData='' ){
  #print(seed)
  for (v in setdiff(names(para_template_biexp), names(para))) para[v] <- para_template_biexp[v]
  
  set.seed(seed)
  
  betas1  = data.frame(
    # gamma = rep(p1['gamma'], N),
    gamma = rep(para['gamma'], N),
    betaEvent_intercept = rep(para['betaEvent_intercept'], N),
    betaEvent_binary = rep(para['betaEvent_trt'], N),
    betaEvent_assoc = rep(para['alpha'], N),
    ks_trt = exp( rep(log(para['ks_trt']), N) + norm1(n=N, mean=0, sd=para['ks_omega']) )* time_multiplier,
    ks_soc = exp( rep(log(para['ks_soc']), N) + norm1(n=N, mean=0, sd=para['ks_omega']) )* time_multiplier,
    kg_trt = exp( rep(log(para['kg_trt']), N) + norm1(n=N, mean=0, sd=para['kg_omega']) )* time_multiplier,
    kg_soc = exp( rep(log(para['kg_soc']), N) + norm1(n=N, mean=0, sd=para['kg_omega']) )* time_multiplier
  )
  
  covdat <- data.frame(
    trt = c(rep(1, N/2), rep(0, N/2))
  )
  

  s2 <- simsurv(hazard = haz_biexp, x = covdat, betas = betas1, maxt = maxt, interval = c(1e-08, 500) / time_multiplier)
  
  # apply practicall follow-up time
  if(!is.null(targetEPS)) s2 <- censorTimeByEPS(s2,  targetEPS=targetEPS)
  
  
  censor_maxfup <-  tapply(1 - s2[,'status'], covdat[,'trt'], sum)
  
  # TODO add censoring, which may affect atrisk status, event time, and event label
  dropout_censor <- 0
  readout_time <- max(s2[,'eventtime'])
  
  if(uniformCensor ){
    dropout_time <- runif(n =N, min = 0, max=readout_time)
    dropout_censor <- dropout_time < s2[,'eventtime'] 
    s2[dropout_censor, 'eventtime'] <- dropout_time[dropout_censor]
    
    dropout_censor <- dropout_censor & s2[, 'status']  == 1
    s2[dropout_censor, 'status'] <- 0
    
  }
  
  if( F && dropRateAtReadout > 0){
    # need update   
    dropout_lambda <- survival_lambda_exponential(targetRate = dropRateAtReadout, t = readout_time)
    dropout_time <- rexp(n=N, rate =dropout_lambda)  
    dropout_censor <- dropout_time < s2[,'eventtime']
    s2[dropout_censor, 'eventtime'] <- dropout_time[dropout_censor]
    s2[dropout_censor, 'status'] <- 0
  }
  
  
  
  # biomarker value, only use SE at the landmark visit
  tr1 <- haz_biexp(x = covdat, betas = betas1, returnTrajectoryUpto = lm)
  if(returnData == 'trajectory') return(tr1)
  
  bm <- tr1[, ncol(tr1)]   # the last column is the biomarker value at the landmark time
  
  if('measure_sd' %in% names(para) && !is.na(para['measure_sd']) && para['measure_sd'] > 0) bm <- bm + rnorm(n=length(bm), mean=0, sd = para['measure_sd'])
  
  ds1 <- cbind(cbind(s2, covdat) , bm= bm )
  
  ds1$t_lm <- ds1$eventtime - lm
  ds1$atrisk = with(ds1, eventtime > lm)
  
  if(returnData == 'landmark') return(ds1)
  
  if(!any(ds1$atrisk)) {
    stop('no patients at risk')
    
  }
  
  # summarize censor / early events by arm
  
  cen1 <- with(ds1[ds1$atrisk, ], tapply(1 - status, trt, sum))
  #  cen1 <- with(ds1[ds1$atrisk, ], sum(1 - status))
  early_event <- with(ds1, tapply(1 - atrisk, trt, sum))
  
  # alternative ways to code biomarker 
  ds1$bm2 <- as.integer(ds1$bm >= median(ds1$bm))
  
  # impute biomarker values for patients who are not at risk at the landmark time (even though we have the simulated values
  #  ds1$bm_impute <- with(ds1, ifelse(atrisk, bm, max(bm)))
  
  # ingore (not impute) biomarker values that may be missing for patients who are not at risk
  if(imputeCensoredSE){
    ds1$bm_impute <- with(ds1, ifelse(atrisk, bm, max(bm)))
    bm_diff <- with(ds1, tapply(bm_impute, trt, median))

    bm_t <- wilcox.test(bm_impute ~ trt, data=ds1, alternative = 'greater', exact=F)
  }else{  
    bm_diff <- with(ds1, tapply(bm, trt, mean)) 
  
    bm_t <- t.test(bm ~ trt, data=ds1) 
  }
  
  k1 <- summary(survfit(Surv(eventtime, status) ~ trt, ds1))$table[,'median']
  
  r_ops <- c(seed = seed,
             N_overall_censor_atrisk = cen1, n_early_event = early_event, max_fup = readout_time, 
             N_censor_dropout = sum(dropout_censor), N_censor_maxfup = censor_maxfup, medianTime=k1
  )
  
  if(operationSummaryOnly){
    return(# about trial operation
      r_ops
      )
  }
  
  if(reportPDchangeOnly){
    stopifnot(imputeCensoredSE)
   return(c(r_ops,    
            'bm_diff' = bm_diff['1'] - bm_diff['0'] ,  # median different
            'bm_t' =   bm_t$p.value   # pvale of wilcoxon test   
              )
   )
  }
  
  
  c1 <- coxph(Surv(eventtime, status) ~ trt , ds1)
  
  #  c2 <- coxph(Surv(t_lm, status) ~ trt , ds1[ ds1$atrisk, ])
  
  c3 <- coxph(Surv(t_lm, status) ~ bm + trt , ds1[ ds1$atrisk, ])
  c3a <- coxph(Surv(t_lm, status) ~ bm + strata(trt) , ds1[ ds1$atrisk, ])
  c3b <- coxph(Surv(t_lm, status) ~ bm , ds1[ ds1$atrisk, ])
  c3c <- coxph(Surv(t_lm, status) ~ bm , ds1[ ds1$atrisk & ds1$trt==1, ])
  c3d <- coxph(Surv(t_lm, status) ~ bm, ds1[ ds1$atrisk & ds1$trt==0, ])
  
  c4a <- coxph(Surv(t_lm, status) ~ bm2 + strata(trt) , ds1[ ds1$atrisk, ])
  # extract c index, harrell's
  
  #c_h <- concordance(object=c3a)$concordance
  c_rfs <- with(ds1[ ds1$atrisk, ], rfsrc_cindex(t_lm, status, bm * -1))
  c_rfs1 <- with(ds1[ ds1$atrisk & ds1$trt == 1, ], rfsrc_cindex(t_lm, status, bm* -1))
  c_rfs0 <- with(ds1[ ds1$atrisk & ds1$trt == 0, ], rfsrc_cindex(t_lm, status, bm* -1))
  
  
  # extract uno's c index
  
  unoc = with(ds1[ ds1$atrisk, ], UnoC(Surv.rsp = Surv(t_lm, status) , Surv.rsp.new = Surv(t_lm, status) , lpnew =bm ))
  
  unoc1 = with(ds1[ ds1$atrisk & ds1$trt == 1, ], UnoC(Surv.rsp = Surv(t_lm, status) , Surv.rsp.new = Surv(t_lm, status) , lpnew =bm ))
  unoc0 = with(ds1[ ds1$atrisk & ds1$trt == 0, ], UnoC(Surv.rsp = Surv(t_lm, status) , Surv.rsp.new = Surv(t_lm, status) , lpnew =bm ))
  
  # get IBS from pec
  
  ibs = ibs_pec(data =ds1[ ds1$atrisk, ] , tte = 't_lm', cens = 'status', x = 'bm')
  ibs1 = ibs_pec(data =ds1[ ds1$atrisk  & ds1$trt == 1, ] , tte = 't_lm', cens = 'status', x = 'bm')
  ibs0 = ibs_pec(data =ds1[ ds1$atrisk  & ds1$trt == 0, ] , tte = 't_lm', cens = 'status', x = 'bm')
  
  ibs_strata = ibs_strata(data =ds1[ ds1$atrisk, ] , tte = 't_lm', cens = 'status', x = 'bm', strat.var = 'trt')
  
  c(# about trial operation
    r_ops,
    
    # about performance measurement    
    # c3$coefficients['bm2'], c3a$coefficients['bm2'],c3b$coefficients['bm2'],c3c$coefficients['bm2'],c3d$coefficients['bm2'],  # binary biomarker
    
    
    #bmHR_adjusted_by_trt = c3$coefficients['bm'], 
    bmHR_stratifi_by_trt = c3a$coefficients['bm'],
    #   bmHR_ingore_trt= c3b$coefficients['bm'],
    bmHR_trt1 = c3c$coefficients['bm'],
    bmHR_trt0 = c3d$coefficients['bm'], 
    
    bm2HR_stratified = c4a$coefficients['bm2'],
    
    #c_h_strata = c_h, 
    c_rfs = c_rfs, c_rfs_trt1 = c_rfs1, c_rfs_trt0 = c_rfs0,
    unoC=unoc, unoC_trt1 = unoc1, unoC_trt0 = unoc0, 
    ibs = ibs,   ibs_trt1 = ibs1, ibs_trt0= ibs0, 
    ibs_strata =ibs_strata, 
    
    # to support meta-regression
    'trt_itt' = c1$coefficients['trt'],    
    #   'trt_atrisk' = c2$coefficients['trt'],
    'bm_diff' = bm_diff['1'] - bm_diff['0'], 
    'bm_t' = bm_t$statistic,
    'bm_p' = bm_t$p.value
  )
}

#=default parameters=============

#para_template=c('gamma' = 1,betaEvent_intercept = log(0.056),  betaEvent_trt = 0, alpha=-2.5, bm_trt=-0.8, bm_ctrl=-0.4, sd_long= 0.3)

para_template=c('gamma' = 1,betaEvent_intercept = -1.5,  betaEvent_trt = 0, alpha=2.5, bm_trt=-0.8, bm_ctrl=-0.4, sd_long= 0.3)

para_template_biexp=c('gamma' = 1,betaEvent_intercept = -1.5,  betaEvent_trt = 0, alpha=2.5, ks_trt = 0.08,  ks_soc= 0.08, kg_trt=0.02, kg_soc=0.02, ks_omega = sqrt(0.8), kg_omega = sqrt(0.6 ), measure_sd= 0.09)

pair_simulation1 <- function(ds =r1,  f1, returnSimData=F ) {
  #each seed contribute to one different value of bm_trt / beta1
  stopifnot( 'seed' %in% colnames(f1))
  stopifnot(all(colnames(f1) %in% colnames(ds)))  
  x <- merge(ds, f1)     
  
  train <- merge(subset(ds, ! seed %in% f1[, 'seed'] ), f1[, setdiff(colnames(f1),'seed'), drop=F])
  train_sample <- train[sample(nrow(train), size=1),,drop=F]
  
  
  if(returnSimData) {
    print( paste('simulation includes', nrow(x), 'studies') )
    x
  }else{
    
    lm1 <- summary(lm(trt_itt.trt ~	bm_diff.1, data=x))
    
    r3 <- with(x, 
               c( #c_rfs=mean(c_rfs), ibs= mean(ibs), beta=mean(bmHR_stratifi_by_trt.bm), 
                 c_rfs=train_sample[, 'c_rfs'], ibs= train_sample[, 'ibs.original'],ibs.scaled= train_sample[, 'ibs.scaled'], beta = train_sample[, 'bmHR_stratifi_by_trt.bm'], 
                 cor1=cor(trt_itt.trt,	bm_diff.1 ), r.squared=lm1[["r.squared"]], intercept=lm1$coefficients['(Intercept)','Estimate'], slope = lm1$coefficients['bm_diff.1','Estimate'], n=nrow(x) ))
    
    for (e in c('alpha','beta1','bm_trt') ) {
      if(e %in% colnames(x) && length( unique(as.numeric(x[,e])) ) ==1 ) r3[e] <- as.numeric(unique(x[,e]))
    }
    
    r3
  }
}



# this function will assemble simulations stratified by stratum variables (which are likely simulation parameters to affect performance  or SE utility)

pair_simulation2 <- function(ds =r1,  pd_mod.var  ='bm_trt', dup = 1, strata.var=c('alpha','beta1'), returnSimData=F, verbose = 0 ) {
  
  #each seed contribute to one different value of bm_trt / beta1
  stopifnot( all(c(pd_mod.var, strata.var, 'seed') %in% colnames(ds) ))
  
  s1 <- lapply( strata.var, function(v){
      u <- unique( ds[, v] )
      if( length(u ) > 1 ){
        if (verbose >= 1) print(paste(v, 'has the following levels to stratify', paste(u, collapse = ',')))
        ds[, v]
      }else{
        return(NULL)
      }
    } )
  
  names(s1) <- strata.var
  
  ds1 <- split(ds, s1[sapply(s1, length) > 0])
  if(verbose >= 1) print(paste('spliting input data into', length(ds1),'strata'))

  exp_b <- data.frame( unique(ds[, pd_mod.var,drop=F]) )
  colnames(exp_b) <- pd_mod.var
  
  if(verbose >= 2) print(exp_b)
  
  s3 <- do.call(rbind, lapply(ds1, function(ds2){
    sl1 <- sample(unique(ds2$seed), size = nrow(exp_b) * dup, replace=F)
    f1_temp <- cbind(exp_b[rep(1:nrow(exp_b), dup), ], sl1)
    colnames(f1_temp)  <- c(colnames(exp_b),'seed')
    
    if(verbose >=1) print( paste('simulation sample template file has', nrow(f1_temp),'rows' ) )
    pair_simulation1 (ds2, f1=f1_temp, returnSimData = returnSimData )
    
  }))
}


# this function assume all elements of vList have the same assessment time points.
multiTrajectory <- function(vList, legPos='bottomleft', increaseYrange=1, ...){
  min_y <- min( as.vector(unlist(vList) ))
  max_y <- max(as.vector(unlist(vList) ))  
  
  s1 <- sapply(2:length(vList), function(i){all( vList[[i]][2, ] == vList[[1]][2, ])})
  
  stopifnot(all(s1))
  
  tr1 <- as.matrix( sapply(vList, function(x)x[1, ]) )
  x <- as.numeric(row.names(tr1))
  
  plot(x=x, y=vList[[1]][2, ], col='black', type='b', xlim=range(x), ylim=c(min_y, max_y) * increaseYrange, ylab= 'SE value', xlab='month', ...)
  col1 <- rainbow(ncol(tr1))
  for (i in 1:ncol(tr1)) {
    lines(x=x, y=tr1[,i], type ='b', col=col1[i])
  }
  
  if(!is.na(legPos)) legend(legPos, pch=1, ncol=2, col=c('black',col1), legend=c('SOC', names(vList)))
}


library(survival)
library(simsurv)
if(F){

library(survAUC)
library(pec)
library(knitr)
  
  source('~/R/lib_2020/output_functions.r')

#ftime <- fupByEPS_weibull (targetEPS=0.75, maxFollowUp=200,  para=c(alpha=2, sd_long=0.4))  

# with sd_long=0.4, there is reasonable fup for alpha = 0 and alpha = 4

#ftime <- fupByEPS_weibull (targetEPS=0.75, maxFollowUp=200,  para=c(alpha=0, sd_long=0.4)) 

test1 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=300, para=c(sd_long=0.001),operationSummaryOnly=T)) )

test2 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=300, para=c(),operationSummaryOnly=T)) )

t( sapply(1:5, function(i) sim1(seed= i+27519, N=300, para=c(betaEvent_intercept=0),operationSummaryOnly=T)) )

change_a_0 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=300, para=c(alpha=0))) )

change_a_2 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=300, para=c( alpha=2))) )

change_a_4 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=300, para=c( alpha=4))) )

change_a_6 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=300, para=c(betaEvent_intercept=0, alpha=6, bm_trt=-0.6, bm_trt=-0))) )
 

knitr::kable(summarizeMatrix(change_a_0), digits=3)
knitr::kable(summarizeMatrix(change_a_4), digits=3)


change_b_1 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=200, para=c( alpha=2, betaEvent_trt=-1))) )
 
change_b_2 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=200, para=c( alpha=2, betaEvent_trt=-2 ))) )


change_chg_0 <- t( sapply(1:50, function(i) sim1(seed= i+27519, N=200, para=c( alpha=4, bm_trt=-0.4))) )


if(F){
# experiment to change follow up time  

ftime <- fupByEPS_weibull (targetEPS=0.75, maxFollowUp=300,  para=c('betaEvent_trt'= -0.25))  
s2 <- t( sapply(1:50, function(i) sim1(seed= i+27519, maxt =ftime, N=200, para=c('betaEvent_trt'= -1))) )

#experiment to change sd_long. 

s2 <- t( sapply(1:100, function(i) sim1(seed= i+27519, N=100,sd_long=1)) )

'with N=600, change long_sd from 0.1 to 1 brings the following changes
* more event (bad prognosis): less censored events by max followup and much more early event
* cox coefficient for bm remains the same
* cox treatment effect is reduced (from -1.1 to -0.4)
* treatment effect on bm remains the same
* c index improves a lot (0.4-> 0.1) and ibs improves numerically less: 0.1->0.075'

}




}