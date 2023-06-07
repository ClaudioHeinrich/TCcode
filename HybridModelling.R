

# get skill of hybrid models: predict SSTs (IOD or N34), then regress on that:

rm(list = ls())

devtools::load_all()
library(data.table)
library(ggplot2)

standardize_precip = FALSE # do you want to standardize precip and precip forecasts before working with them?

### get datasets ###

plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/'
data_dir = '/nr/project/stat/CONFER/Data/'

prec_dt = fread(paste0(data_dir,'monthly_mean_prec/OND_ON_prec_fc_init_mon8.csv')) # beware: cmcc and dwd have only forecasts up to November

# reduce to the years we look at:

years = 1993:2016
prec_dt = prec_dt[year %in% years]
prec_dt = prec_dt[season == 'ON']

# mask low precip regions:
prec_dt = mask_precip(prec_dt,value_col = 'obs',bycols = c('lon','lat','season'))
prec_dt = prec_dt[!(mask)][,mask:= NULL]


# standardize forecasts:

if(standardize_precip)
{
  prec_dt = loyo_standardize(prec_dt,SDcols = c('prec','obs'),bycols = c('lon','lat','season','system'))
  prec_dt[,prec_pp := prec]
}

IOD_dt = fread(paste0(data_dir,'IOD.csv'))
IOD_dt = IOD_dt[,c('IOD_west','IOD_east') :=lapply(.SD,function(x) (x - mean(x))/sd(x)),.SDcols = c('IOD_west','IOD_east'),by = month]
IOD_dt[,IOD := IOD_west - IOD_east]

N34_dt = fread(paste0(data_dir,'N34.csv'))
N34_dt = N34_dt[,'N34' := lapply(.SD,function(x) (x - mean(x))/sd(x)),.SDcols = c('N34'),by = month]

SST_dt = fread(paste0(data_dir,'/N34_IOD_forecasts_raw.csv'))
SST_dt = SST_dt[forecast_month == 8][,forecast_month := NULL]  
setnames(SST_dt,2:3,c('year','month'))
SST_dt = SST_dt[month %in% 9:11][year %in% years]
SST_dt = loyo_standardize(SST_dt,
                          SDcols = c('N34','IOD_west','IOD_east'),
                          bycols = c('month','system'))


SST_dt[,IOD := IOD_west-IOD_east]
SST_dt[,c('IOD_west','IOD_east'):=NULL]
SST_dt = SST_dt[,(lapply(.SD,mean,na.rm = T)),.SDcols = c('N34','IOD'),by = .(month,system,year)]
SST_dt_mixed = SST_dt[,(lapply(.SD,mean,na.rm = T)),.SDcols = c('N34','IOD'),by = .(month,year)]
SST_dt_mixed[,system := 'mixed']
SST_dt = rbindlist(list(SST_dt,SST_dt_mixed),use.names = T)

### regress ON precip ###
obs_dt = fread(paste0(data_dir,'CHIRPS_prec_upscaled.csv'))
obs_dt = obs_dt[month %in% 10:11,.(prec = mean(prec)),.(year,lon,lat)]

# reduce to correct region:
obs_dt = obs_dt[!is.na(prec)]
setkey(obs_dt,lon,lat)
obs_dt = obs_dt[GHA_locs()]

obs_dt = mask_precip(obs_dt)[(!mask)][,mask := NULL]

###

prec_dt = merge(prec_dt,IOD_dt[month == 8][,.(year,IOD)],by = c('year'))
setnames(prec_dt,'IOD','IOD_obs_Aug')
prec_dt = merge(prec_dt,N34_dt[month == 8][,.(year,N34)],by = c('year'))
setnames(prec_dt,'N34','N34_obs_Aug')

prec_mean_dt = prec_dt[,lapply(.SD,mean),
                       .SDcols = c('prec','sys_clim','sys_clim_sd','obs','clim','clim_sd','prec_pp','IOD_obs_Aug','N34_obs_Aug'),
                       by = c('year','lon','lat','system')]


### run a bunch of regressions: ###


### reduce IOD error ###
dt_temp = copy(prec_mean_dt)
dt_temp[,error := prec_pp - obs]

test = loyo_dt_lm(dt_temp,formula = error ~ IOD_obs_Aug ,bycols = c('system','lon','lat'),mc_cores = 8)
setkey(test,year,system,lon,lat)
setkey(prec_mean_dt,year,system,lon,lat)

prec_mean_dt[,prec_pp_iod := prec_pp - test[,prediction]]

### reduce N34 error ###
dt_temp = copy(prec_mean_dt)
dt_temp[,error := prec_pp - obs]

test = loyo_dt_lm(dt_temp,formula = error ~ N34_obs_Aug ,bycols = c('system','lon','lat'),mc_cores = 8)
setkey(test,year,system,lon,lat)
setkey(prec_mean_dt,year,system,lon,lat)

prec_mean_dt[,prec_pp_n34 := prec_pp - test[,prediction]]



# pp_dt = loyo_dt_lm(prec_mean_dt,formula = obs~prec,bycols = c('system','lon','lat'),mc_cores = 8)
# setnames(pp_dt,'prediction','prec_pp_new')
# prec_mean_dt = merge(prec_mean_dt,pp_dt[,.SD,.SDcols = c('year','system','lon','lat','prec_pp_new')],c('year','system','lon','lat'))
# 
# IOD_pp_dt = loyo_dt_lm(prec_mean_dt,formula = obs~prec + IOD_obs_Aug,bycols = c('system','lon','lat'),mc_cores = 8)
# setnames(IOD_pp_dt,'prediction','prec_pp_iod')
# prec_mean_dt = merge(prec_mean_dt,IOD_pp_dt[,.SD,.SDcols = c('year','system','lon','lat','prec_pp_iod')],c('year','system','lon','lat'))
# 
# N34_pp_dt = loyo_dt_lm(prec_mean_dt,formula = obs~prec + N34_obs_Aug,bycols = c('system','lon','lat'),mc_cores = 8)
# setnames(N34_pp_dt,'prediction','prec_pp_n34')
# prec_mean_dt = merge(prec_mean_dt,N34_pp_dt[,.SD,.SDcols = c('year','system','lon','lat','prec_pp_n34')],c('year','system','lon','lat'))
# 
### assess skill ###

skill_dt1 = MSE_dt(prec_mean_dt,fc_col = 'prec_pp',obs_col = 'obs',by_cols = c('lon','lat','year','system'))[,pp := 'standard']
skill_dt2 = MSE_dt(prec_mean_dt,fc_col = 'prec_pp_iod',obs_col = 'obs',by_cols = c('lon','lat','year','system'))[,pp := 'iod']
skill_dt3 = MSE_dt(prec_mean_dt,fc_col = 'prec_pp_n34',obs_col = 'obs',by_cols = c('lon','lat','year','system'))[,pp := 'n34']

skill_dt = rbindlist(list(skill_dt1,
                          skill_dt2,
                          skill_dt3))


test_dt = skill_dt[,.(MSE = mean(MSE)),by = c('system','pp')]

test_dt[system == 'dwd']
test_dt[system == 'ecmwf']
test_dt[system == 'meteo_france']

bt_skill = bootstrap_scores_dt()

#### try eof regression for this ####
temp_dt = data.table()
for(sys in unique(prec_dt[,system]))
{
  ttemp = loyo_eof_regression(dt_temp[system == sys],formula = error ~ 0+IOD_obs_Aug ,nv = 1,mc_cores = 8)
  ttemp[,system := sys]
  temp_dt = rbindlist(list(temp_dt,ttemp))
}

setkey(temp_dt,year,system,lon,lat)
setkey(prec_mean_dt,year,system,lon,lat)

prec_mean_dt[,prec_pp_iod_eof := prec_pp - temp_dt[,prediction]]

# N3.4 

temp_dt = data.table()
for(sys in unique(prec_dt[,system]))
{
  ttemp = loyo_eof_regression(dt_temp[system == sys],formula = error ~ 0+N34_obs_Aug ,nv = 1,mc_cores = 8)
  ttemp[,system := sys]
  temp_dt = rbindlist(list(temp_dt,ttemp))
}

setkey(temp_dt,year,system,lon,lat)
setkey(prec_mean_dt,year,system,lon,lat)

prec_mean_dt[,prec_pp_n34_eof := prec_pp - temp_dt[,prediction]]


skill_dt4 = MSE_dt(prec_mean_dt,fc_col = 'prec_pp_iod_eof',obs_col = 'obs',by_cols = c('lon','lat','year','system'))[,pp := 'iod_eof']
skill_dt5 = MSE_dt(prec_mean_dt,fc_col = 'prec_pp_n34_eof',obs_col = 'obs',by_cols = c('lon','lat','year','system'))[,pp := 'n34_eof']

skill_dt = rbindlist(list(skill_dt1,
                          skill_dt2,
                          skill_dt3,
                          skill_dt4,
                          skill_dt5))


test_dt = skill_dt[,.(MSE = mean(MSE)),by = c('system','pp')]

test_dt[system == 'dwd']
test_dt[system == 'ecmwf']
test_dt[system == 'meteo_france']
test_dt[system == 'cmcc']

#############################################################################
#### run regression models on observed N34 and IOD and append to prec_dt ####
#############################################################################

# ### N34 observed in August ###
# 
# reg_dt = data.table()
# 
# for(yy in years)
# {
#   N34_reg = N34_dt[year >= 1993 & month ==8 & year != yy & year %in% years,N34]
#   
#   coeff_dt = matrix_linear_regression(N34_reg,yy)
#   coeff_dt[,system := 'N34_observed_August'][,year := yy]
#   
#   reg_dt = rbindlist(list(reg_dt,coeff_dt))
# 
# }
# 
# reg_dt = dcast(reg_dt,lon + lat + season + system + year ~ coef)
# 
# reg_dt = merge(reg_dt,N34_dt[year >= 1993 & month ==8 ,.(year,N34)],by = 'year')
# if('intercept' %in% names(reg_dt)) 
# {reg_dt[,prec := intercept + beta1*N34][,prec_pp := prec]
#   reg_dt[,c('beta1','intercept','N34') := NULL]}
# if(!('intercept' %in% names(reg_dt))) 
# {reg_dt[,prec := beta1*N34][,prec_pp := prec]
#   reg_dt[,c('beta1','N34') := NULL]}
# 
# reg_dt = merge(reg_dt,obs_dt,by = c('year','lon','lat','season'))
# 
# prec_dt = rbindlist(list(prec_dt,reg_dt),use.names= T)
# 
# ### IOD observed in August ###
# 
# reg_dt = data.table()
# 
# for(yy in years)
# {
#   IOD_reg = IOD_dt[year >= 1993 & month ==8 & year != yy & year %in% years,IOD]
#   
#   coeff_dt = matrix_linear_regression(IOD_reg,yy)
#   coeff_dt[,system := 'IOD_observed_August'][,year := yy]
#   
#   reg_dt = rbindlist(list(reg_dt,coeff_dt))
#   
# }
# 
# reg_dt = dcast(reg_dt,lon + lat + season + system + year ~ coef)
# 
# reg_dt = merge(reg_dt,IOD_dt[year >= 1993 & year %in% years & month ==8 ,.(year,IOD)],by = 'year')
# if('intercept' %in% names(reg_dt)) 
# {reg_dt[,prec := intercept + beta1*IOD][,prec_pp := prec]
#   reg_dt[,c('beta1','intercept','IOD') := NULL]}
# if(!('intercept' %in% names(reg_dt))) 
# {reg_dt[,prec := beta1*IOD][,prec_pp := prec]
#   reg_dt[,c('beta1','IOD') := NULL]}
# 
# reg_dt = merge(reg_dt,obs_dt,by = c('year','lon','lat','season'))
# 
# prec_dt = rbindlist(list(prec_dt,reg_dt),use.names= T)
# 
# test = prec_dt[,mean((obs - prec_pp)^2),.(system,season)]
# 
# 



############################################################################################
#### predict SSTs (IOD, N34), then run a regression model for precip on predicted SSTs  ####
############################################################################################

prec_dt = prec_dt[system != 'mixed']
SST_dt

iod_dt_temp = dcast(SST_dt,system + year + member ~ month,value.var = 'IOD')
setnames(iod_dt_temp,4:6,paste0('iod_',9:11))
n34_dt_temp = dcast(SST_dt,system + year +member~ month,value.var = 'N34')
setnames(n34_dt_temp,4:6,paste0('n34_',9:11))

prec_dt_new = merge(prec_dt,iod_dt_temp,by =c('system','year','member'))
prec_dt_new = merge(prec_dt_new,n34_dt_temp,by =c('system','year','member'))

# get error correlations by months:

prec_dt_new[,error:= prec_pp - obs]

MSE_pred_pp = data.table()

for(mon in 9:11)
{
  for(tc in c('iod','n34'))
  {
    
  print(paste0(mon,', ',tc))
  ff = as.formula(paste0('error ~ 0 + ',tc,'_',mon))
  temp = loyo_dt_lm(prec_dt_new,formula = ff ,bycols = c('system','lon','lat'),alongcols = c('year','member'),mc_cores = 8)
  temp[,prec_pp := prec_pp - prediction]
  temp = temp[,lapply(.SD,mean),.SDcols = c('obs','prec_pp'),by = c('system','year','lon','lat')]
  MSE_temp = MSE_dt(temp,fc_col = 'prec_pp',obs_col = 'obs')
  MSE_temp[,formula:=as.character(Reduce(paste, deparse(ff)))][,pred_month := mon]
  
  MSE_pred_pp = rbindlist(list(MSE_pred_pp,MSE_temp))
  }
}

test = MSE_pred_pp[,mean(MSE),by = c('system','formula','pred_month')]
test[system == 'ecmwf']

cor_dt = unique(prec_dt_new[,.(system,lon,lat,cor_iod_9,cor_iod_10,cor_iod_11,cor_n34_9,cor_n34_10,cor_n34_11)])
rr = c(-0.75,0.75)

systems = unique(prec_dt_new[,system])


for(ss in systems)

{
  mm = 11
  pp = ggplot_dt(cor_dt[system == ss],paste0('cor_iod_',mm),centering = 0,rr = rr) + ggtitle(ss)
  print(pp)
  
  pp = ggplot_dt(cor_dt[system == ss],paste0('cor_n34_',mm),centering = 0,rr = rr) + ggtitle(ss)
  print(pp)
}

# what if we take ensemble means?

prec_dt_new = merge(prec_dt,iod_dt_temp,by =c('system','year','member'))
prec_dt_new = merge(prec_dt_new,n34_dt_temp,by =c('system','year','member'))


prec_dt_new = prec_dt_new[,lapply(.SD,mean),.SDcols = colnames(prec_dt_new)[7:ncol(prec_dt_new)],by = .(system,year,lon,lat)]
# get error correlations by months:

prec_dt_new[,error:= prec_pp - obs]

for(mon in 9:11)
{
  prec_dt_new[,paste0('cor_iod_',mon) := cor(error,get(paste0('iod_',mon))),by = .(system,lon,lat)]
  prec_dt_new[,paste0('cor_n34_',mon) := cor(error,get(paste0('n34_',mon))),by = .(system,lon,lat)]
}

cor_dt = unique(prec_dt_new[,.(system,lon,lat,cor_iod_9,cor_iod_10,cor_iod_11,cor_n34_9,cor_n34_10,cor_n34_11)])
rr = c(-0.75,0.75)

systems = unique(prec_dt_new[,system])


for(ss in systems)
  
{mm = 11
pp = ggplot_dt(cor_dt[system == ss],paste0('cor_iod_',mm),centering = 0,rr = rr) + ggtitle(ss)
print(pp)

pp = ggplot_dt(cor_dt[system == ss],paste0('cor_n34_',mm),centering = 0,rr = rr) + ggtitle(ss)
print(pp)
}

