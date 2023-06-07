
rm(list = ls())

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()

### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo')


init_mons = 8

data_dir = '/nr/project/stat/CONFER/Data/monthly_mean_prec/'
plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/Precip_skill/'

#################################

### assess raw forecast skill ###

prec_dt = fread(paste0(data_dir,'monthly_prec_fc_init_mon',init_mons,'.csv'))

prec_dt = prec_dt[forecast_year < 2021]
prec_dt_temp = prec_dt[target_month %in% c(10:12)][,.(prec = mean(total_precipitation)),by = .(lon,lat,system,forecast_year,member)]
prec_dt_temp[,season := 'OND']

prec_dt_temp2 = prec_dt[target_month %in% c(10:11)][,.(prec = mean(total_precipitation)),by = .(lon,lat,system,forecast_year,member)]
prec_dt_temp2[,season := 'ON']

prec_dt = rbindlist(list(prec_dt_temp,prec_dt_temp2))

setnames(prec_dt,'forecast_year','year')

# convert to mm/day

prec_dt[,prec := prec*1000*3600*24]

# get a mixed ensemble forecast, for which each member is simply the ensemble mean for one of the forecast systems

temp_dt  = prec_dt[,.(prec = mean(prec)),by = .(lon,lat,system,year,season)]
temp_dt[,member := 1:.N,by = .(lon,lat,year)][,system := 'mixed']

prec_dt = rbindlist(list(prec_dt,temp_dt),use.names = T)

systems = c(systems,'mixed')

# get loyo model climatologies and model climatology sds

if(!('sys_clim' %in% colnames(prec_dt))) 
{
  prec_dt = loyo(prec_dt,
                 FUN = mean,
                 SDcols = 'prec',
                 bycols = c('lon','lat','system','season'),
                 na.rm = T)
  setnames(prec_dt,'prec_new','sys_clim')
  
}


if(!('sys_clim_sd' %in% colnames(prec_dt)))
{
  prec_dt = loyo(prec_dt,
                 FUN = sd,
                 SDcols = 'prec',
                 bycols = c('lon','lat','system','season'),
                 na.rm = T)
  
  setnames(prec_dt,'prec_new','sys_clim_sd')
}


### merge with obsercation ###


obs_dt = load_chirps()
obs_dt[,index_coarse := NULL]

obs_dt1 = obs_dt[month %in% 10:12, .(obs = mean(prec)),by = .(lon,lat,year)][,season := 'OND']
obs_dt2 = obs_dt[month %in% 10:11, .(obs = mean(prec)),by = .(lon,lat,year)][,season := 'ON']

obs_dt = rbindlist(list(obs_dt1,obs_dt2))

prec_dt = merge(prec_dt,obs_dt)


# get climatology and clim_sd for observation:

prec_dt = loyo(prec_dt,
               FUN = mean,
               SDcols = 'obs',
               bycols = c('lon','lat','season'),
               na.rm = T)
setnames(prec_dt,'obs_new','clim')

prec_dt = loyo(prec_dt,
               FUN = sd,
               SDcols = 'obs',
               bycols = c('lon','lat','season'),
               na.rm = T)
setnames(prec_dt,'obs_new','clim_sd')


# recalibrate or post-process by using standardized variables:

prec_dt[,prec_pp := clim_sd/sys_clim_sd * (prec - sys_clim) + clim]


fwrite(prec_dt,file = paste0(data_dir,'OND_ON_prec_fc_init_mon',init_mons,'.csv'))