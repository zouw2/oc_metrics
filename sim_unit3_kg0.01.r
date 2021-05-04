

input <- commandArgs(trailingOnly = TRUE)
print(input)

spec <- list(N=400, landmarkTime=2, home = '~/aPDL1/oak') # time unit is week now

#input <- c('alpha', '1', '27519','10','change_alpha2')


library(survival)
library(simsurv)
library(survAUC)
library(pec)

spec$seed_base <- as.integer(input[3])
spec$niter <- as.integer(input[4])

para_name <- strsplit(input[1], split=',')[[1]]

para_value <- as.numeric(strsplit(input[2], split=',')[[1]])

stopifnot(length(para_name) == length(para_value))

source('~/aPDL1/oak/simsurv_functions.r')


para_template_biexp=c('gamma' = 1,betaEvent_intercept = -2,  betaEvent_trt = 0, alpha=2.5, 
                      ks_trt = 0.02,  ks_soc= 0.02, 
                      kg_trt=0.01, kg_soc=0.01, ks_omega = sqrt(0.8), kg_omega = sqrt(0.6 ), measure_sd= 0.09)


stopifnot(all(para_name %in% names(para_template_biexp)))

for (i in 1:length(para_name)) para_template_biexp[para_name[i]] <- para_value[i]

print('parameter values modified here:')

print(data.frame(value=para_template_biexp[para_name]))

print('conditional these parameter values from the template')
print(data.frame(value= para_template_biexp[setdiff(names(para_template_biexp), para_name)]))


s1 <- t( sapply(1:spec$niter, function(i) sim2_be(lm =spec$landmarkTime, seed= i+spec$seed_base, N=spec$N ,maxt = 120)) )


save(s1,para_template_biexp, file = paste(spec$home,'/results/',input[5],'/', paste(gsub(',','_', input[1:3], fixed=T), collapse = '_'), '.Rdata', sep='' ) )

