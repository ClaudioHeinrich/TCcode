

# get skill of hybrid models: predict SSTs (IOD or N34), then regress on that:

rm(list = ls())

devtools::load_all()
library(data.table)
library(ggplot)

standardize_precip = FALSE # do you want to standardize precip and precip forecasts before working with them?

### get datasets ###

plot_dir = '/nr/project/stat/CONFER/plots/TCmodels/'
data_dir = '/nr/project/stat/CONFER/Data/'

prec_dt = fread(paste0(data_dir,'monthly_mean_prec/OND_ON_prec_fc_init_mon8.csv')) # beware: cmcc and dwd have only forecasts up to November

# reduce to the years we look at:

years = 1993:2016
prec_dt = prec_dt[year %in% years]
prec_dt = prec_dt[season == 'ON']

# mask low precip regions:
prec_dt = mask_precip(prec_dt,value_col = 'obs',bycols = c('lon','lat','season'))
prec_dt = prec_dt[!(mask)][,mask:= NULL]


# recompute (oos) climatology and clim_sd:

prec_dt = loyo(prec_dt,FUN = mean,SDcols = 'obs',bycols = c('lon','lat'))
prec_dt[,clim := NULL]
setnames(prec_dt,'obs_new','clim')

prec_dt = loyo(prec_dt,FUN = sd,SDcols = 'obs',bycols = c('lon','lat','system','member'))
prec_dt[,clim_sd := NULL]
setnames(prec_dt,'obs_new','clim_sd')

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
  SST_dt = SST_dt[month %in% 10:12]
  SST_dt = loyo_standardize(SST_dt,
                            SDcols = c('N34','IOD_west','IOD_east'),
                            bycols = 'month')

  
years = 1993:2016

obs_dt = unique(prec_dt[system == 'ecmwf' & year %in% years,.(lon,lat,season,year,obs,clim,clim_sd)]) # ecmwf has all years

    
matrix_linear_regression = function(covariate_vector,yy, intercept = ifelse(standardize_precip,yes = FALSE,no = TRUE))
{
  # leave-one-year-out matrix linear regression
  
  
  test = dcast(obs_dt[year != yy],year ~ season + lon + lat,value.var = 'obs')
  
  if(intercept) testmod = lm(as.matrix(test[,2:ncol(test)]) ~ covariate_vector)
  if(!intercept) testmod = lm(as.matrix(test[,2:ncol(test)]) ~ 0 + covariate_vector )

    
  reg_coef_dt = as.data.table(testmod$coefficients)
  
  if(intercept) reg_coef_dt[,coef := c('intercept','beta1')]
  if(!intercept) reg_coef_dt[,coef := c('beta1')]
  
  reg_coef_dt = melt(reg_coef_dt,id.vars = 'coef')
  
  reg_coef_dt[,variable:= as.character(variable)]
  
  dt_temp = data.table(matrix(unlist(strsplit(x = reg_coef_dt[,variable],split = '_')),ncol = 3,byrow = T))
  setnames(dt_temp,c('season','lon','lat'))
  dt_temp[,c('lon','lat') := lapply(.SD,as.numeric),.SDcols = c('lon','lat')]
  
  reg_coef_dt = data.table(reg_coef_dt,dt_temp)
  reg_coef_dt[,variable:=NULL]
  
  return(reg_coef_dt)
}


### run a bunch of regressions: ###

# reduce to model means, ensemble verification harder for regression models.

sdcols = c('prec','prec_pp','clim','clim_sd') # note that 'clim' is leave-one-year-out climatology
prec_dt = prec_dt[,(sdcols) := lapply(.SD,mean,na.rm = T),.SDcols = sdcols,by = .(lon,lat,season,year,system)]
prec_dt = unique(prec_dt[,.(lon,lat,season,year,system,prec,prec_pp,obs,clim,clim_sd)])

# get climatology in prediction model format:

clim_dt = unique(prec_dt[,.(lon,lat,season,year,clim,clim_sd,obs)])

clim_dt[,prec := clim]
clim_dt[,prec_pp:=prec]
clim_dt[,system := 'climatology']


prec_dt = rbindlist(list(prec_dt,clim_dt),use.names = T)

setkey(prec_dt,season,system,year,lon,lat)


#############################################################################
#### run regression models on observed N34 and IOD and append to prec_dt ####
#############################################################################

### N34 observed in August ###

reg_dt = data.table()

for(yy in years)
{
  N34_reg = N34_dt[year >= 1993 & month ==8 & year != yy & year %in% years,N34]
  
  coeff_dt = matrix_linear_regression(N34_reg,yy)
  coeff_dt[,system := 'N34_observed_August'][,year := yy]
  
  reg_dt = rbindlist(list(reg_dt,coeff_dt))

}

reg_dt = dcast(reg_dt,lon + lat + season + system + year ~ coef)

reg_dt = merge(reg_dt,N34_dt[year >= 1993 & month ==8 ,.(year,N34)],by = 'year')
if('intercept' %in% names(reg_dt)) 
{reg_dt[,prec := intercept + beta1*N34][,prec_pp := prec]
  reg_dt[,c('beta1','intercept','N34') := NULL]}
if(!('intercept' %in% names(reg_dt))) 
{reg_dt[,prec := beta1*N34][,prec_pp := prec]
  reg_dt[,c('beta1','N34') := NULL]}

reg_dt = merge(reg_dt,obs_dt,by = c('year','lon','lat','season'))

prec_dt = rbindlist(list(prec_dt,reg_dt),use.names= T)

### IOD observed in August ###

reg_dt = data.table()

for(yy in years)
{
  IOD_reg = IOD_dt[year >= 1993 & month ==8 & year != yy & year %in% years,IOD]
  
  coeff_dt = matrix_linear_regression(IOD_reg,yy)
  coeff_dt[,system := 'IOD_observed_August'][,year := yy]
  
  reg_dt = rbindlist(list(reg_dt,coeff_dt))
  
}

reg_dt = dcast(reg_dt,lon + lat + season + system + year ~ coef)

reg_dt = merge(reg_dt,IOD_dt[year >= 1993 & year %in% years & month ==8 ,.(year,IOD)],by = 'year')
if('intercept' %in% names(reg_dt)) 
{reg_dt[,prec := intercept + beta1*IOD][,prec_pp := prec]
  reg_dt[,c('beta1','intercept','IOD') := NULL]}
if(!('intercept' %in% names(reg_dt))) 
{reg_dt[,prec := beta1*IOD][,prec_pp := prec]
  reg_dt[,c('beta1','IOD') := NULL]}

reg_dt = merge(reg_dt,obs_dt,by = c('year','lon','lat','season'))

prec_dt = rbindlist(list(prec_dt,reg_dt),use.names= T)

test = prec_dt[,mean((obs - prec_pp)^2),.(system,season)]



SST_dt[,IOD := IOD_west-IOD_east]

SST_dt[,c('IOD_west','IOD_east'):=NULL]

SST_dt = SST_dt[,(lapply(.SD,mean,na.rm = T)),.SDcols = c('N34','IOD'),by = .(month,system,year)]


######################
#### predict SSTs ####
######################


for(model in c('cmcc','dwd','ecmwf','meteo_france','ukmo'))
    {
      print(model)

############## use predicted SST-values ###########

model_dt = SST_dt[system == model]
model_dt = model_dt[month %in% 10:11,(lapply(.SD,mean)),.SDcols = c('N34','IOD'),by = year]
model_dt[,season := 'ON']

setkey(model_dt,season,year)

### regression model ###

# N34 ON #
reg_dt = data.table()

for(yy in years)
{
  N34_reg = model_dt[season == 'ON' & year >= 1993 & year != yy & year %in% years,N34]
  
  coeff_dt = matrix_linear_regression(N34_reg,yy)
  coeff_dt[,system := paste0('N34_ON_',model)][,year := yy]
  
  reg_dt = rbindlist(list(reg_dt,coeff_dt))
  
}

reg_dt = dcast(reg_dt,lon + lat + season + system + year ~ coef)

reg_dt = merge(reg_dt,model_dt[season == 'ON' & year >= 1993 ,.(year,N34)],by = 'year')
if('intercept' %in% names(reg_dt)) 
  {reg_dt[,prec := intercept + beta1*N34][,prec_pp := prec]
  reg_dt[,c('beta1','intercept','N34') := NULL]}
if(!('intercept' %in% names(reg_dt))) 
  {reg_dt[,prec := beta1*N34][,prec_pp := prec]
  reg_dt[,c('beta1','N34') := NULL]}

reg_dt = merge(reg_dt,obs_dt,by = c('year','lon','lat','season'))

prec_dt = rbindlist(list(prec_dt,reg_dt),use.names= T)




# IOD ON #
reg_dt = data.table()

for(yy in years)
{
  IOD_reg = model_dt[season == 'ON' & year >= 1993 & year != yy & year %in% years,IOD]
  
  coeff_dt = matrix_linear_regression(IOD_reg,yy)
  coeff_dt[,system := paste0('IOD_ON_',model)][,year := yy]
  
  reg_dt = rbindlist(list(reg_dt,coeff_dt))
  
}

reg_dt = dcast(reg_dt,lon + lat + season + system + year ~ coef)

reg_dt = merge(reg_dt,model_dt[season == 'ON' & year >= 1993 ,.(year,IOD)],by = 'year')
if('intercept' %in% names(reg_dt)) 
{reg_dt[,prec := intercept + beta1*IOD][,prec_pp := prec]
  reg_dt[,c('beta1','intercept','IOD') := NULL]}
if(!('intercept' %in% names(reg_dt))) 
{reg_dt[,prec := beta1*IOD][,prec_pp := prec]
  reg_dt[,c('beta1','IOD') := NULL]}

reg_dt = merge(reg_dt,obs_dt,by = c('year','lon','lat','season'))

prec_dt = rbindlist(list(prec_dt,reg_dt),use.names= T)

}

#####################################
### compare results and bootstrap ###
if(standardize_precip) 
  {
  prec_dt[system == 'climatology',prec_pp := 0]
  prec_dt[,prec_pp := prec_pp * clim_sd + clim]
  prec_dt[,obs := obs * clim_sd + clim]
  }

MSE_dt = prec_dt[,.(MSE = (obs - prec_pp)^2),.(system,season,year,lon,lat)]

### bootstrap everything ###

r = 250

pred_mods = unique(MSE_dt[,system])
seasons = unique(MSE_dt[,season])

statistic = function(x,inds){mean(x[inds],na.rm = T)}

bootstrap_all_region_all_years = as.data.table(expand.grid(system = pred_mods,
                                                           season = seasons,
                                                           R = 1:r))

for(pred_mod in pred_mods)
{
  for(ss in seasons)
  {
    print(paste(pred_mod,ss))
    data = MSE_dt[system == pred_mod & season == ss,MSE]
    bootstrap_all_region_all_years[system == pred_mod & season == ss,boots := boot::boot(data,statistic = statistic,R = r)$t]
  }
}

# names = as.character(unique(bootstrap_all_region_all_years[,system]))
# replace_names = c(names[1:7],'N34','IOD','N34_ON','N34_OND','IOD_ON','IOD_OND')
# bootstrap_all_region_all_years[,system := replace_names[match(system,names)]]

# bootstrap_all_region_all_years[,model := factor(model,levels = c('N3.4','IOD','SFE',systems,'N3.4&IOD'))]
# 

model_names_old = c('climatology',
                'cmcc','N34_ON_cmcc','IOD_ON_cmcc',
                'dwd','N34_ON_dwd','IOD_ON_dwd',
                'ecmwf','N34_ON_ecmwf','IOD_ON_ecmwf',
                'meteo_france','N34_ON_meteo_france','IOD_ON_meteo_france',
                'ukmo','N34_ON_ukmo','IOD_ON_ukmo',
                'mixed','N34_observed_August','IOD_observed_August')

model_names_new = c('clim',
                    'cmcc','cmcc_N34','cmcc_IOD',
                    'dwd','dwd_N34','dwd_IOD',
                    'ecmwf','ecmwf_N34','ecmwf_IOD',
                    'mf','mf_N34','mf_IOD',
                    'ukmo','ukmo_N34','ukmo_IOD',
                    'mixed','N34_obs','IOD_obs')

rename = function(x){model_names_new[match(x,model_names_old)]}

bootstrap_all_region_all_years[,system := rename(system)]
setkey(bootstrap_all_region_all_years,system)

theme_set(theme_bw(base_size = 14))

pp = ggplot(bootstrap_all_region_all_years[season == 'ON']) + 
  geom_boxplot(mapping = aes(x = system,y = boots,color = system))
pp = pp + theme(legend.position = 'none') + ylab('MSE bootstrapped')
pp

pp2 = ggplot(bootstrap_all_region_all_years[season == 'ON' & system %in% c('clim',
                                                                           'cmcc','dwd','ecmwf','ecmwf_N34','ecmwf_IOD', 'ukmo','mixed','N34_obs','IOD_obs')]) + 
  geom_boxplot(mapping = aes(x = system,y = boots,color = system))
pp2 = pp2 + theme(legend.position = 'none') + ylab('MSE bootstrapped')
pp2


pdf(paste0(plot_dir,'MSE_bootstrapped.pdf'))
print(pp)
dev.off()



