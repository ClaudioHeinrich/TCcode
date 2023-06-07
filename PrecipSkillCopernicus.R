### Assess the precipitation prediction skill for the different systems 
# In this script we look at seasonal mean precipitation over GHA and assess the skill of the Copernicus models before and after simple postprocessing

rm(list = ls())

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()

### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo')

target_mons = 10:12
init_mons = 8

data_dir = '/nr/project/stat/CONFER/Data/monthly_mean_prec/'
plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/Precip_skill/'

#################################

### assess raw forecast skill ###

prec_dt = fread(paste0(data_dir,'monthly_prec_fc_init_mon',init_mons,'.csv'))

prec_dt = prec_dt[forecast_year < 2021]
prec_dt = prec_dt[target_month %in% c(10:12),.(prec = mean(total_precipitation)),by = .(lon,lat,system,forecast_year,member)]

setnames(prec_dt,'forecast_year','year')

# convert to mm/day

prec_dt[,prec := prec*1000*3600*24]

# get a mixed ensemble forecast, for which each member is simply the ensemble mean for one of the forecast systems

temp_dt  = prec_dt[,.(prec = mean(prec)),by = .(lon,lat,system,year)]
temp_dt[,member := 1:.N,by = .(lon,lat,year)][,system := 'mixed']

prec_dt = rbindlist(list(prec_dt,temp_dt),use.names = T)

systems = c(systems,'mixed')

# get loyo model climatologies and model climatology sds

if(!('sys_clim' %in% colnames(prec_dt))) 
{
  prec_dt = loyo(prec_dt,
                 FUN = mean,
                 SDcols = 'prec',
                 bycols = c('lon','lat','system'),
                 na.rm = T)
  setnames(prec_dt,'prec_new','sys_clim')
  
}
  
    
if(!('sys_clim_sd' %in% colnames(prec_dt)))
{
  prec_dt = loyo(prec_dt,
                 FUN = sd,
                 SDcols = 'prec',
                 bycols = c('lon','lat','system'),
                 na.rm = T)
  
  setnames(prec_dt,'prec_new','sys_clim_sd')
}


### merge with obsercation ###


obs_dt = load_chirps()
obs_dt[,index_coarse := NULL]
obs_dt = obs_dt[month %in% 10:12, .(obs = mean(prec)),by = .(lon,lat,year)]

prec_dt = merge(prec_dt,obs_dt)

# get climatology and clim_sd for observation:

prec_dt = loyo(prec_dt,
               FUN = mean,
               SDcols = 'obs',
               bycols = c('lon','lat'),
               na.rm = T)
setnames(prec_dt,'obs_new','clim')

prec_dt = loyo(prec_dt,
               FUN = sd,
               SDcols = 'obs',
               bycols = c('lon','lat'),
               na.rm = T)
setnames(prec_dt,'obs_new','clim_sd')


# recalibrate or post-process by using standardized variables:

prec_dt[,prec_pp := clim_sd/sys_clim_sd * (prec - sys_clim) + clim]

# get ensemble mean:

prec_dt[,prec_mean := mean(prec,na.rm = T),.(lon,lat,year,system)]
prec_dt[,prec_pp_mean := mean(prec_pp,na.rm = T),.(lon,lat,year,system)]

#### plot correlation of ensemble mean ####

cor_dt = unique(prec_dt[,.(lon,lat,year,system,prec_mean,obs)])

cor_dt = mask_precip(cor_dt,value_col = 'obs')

cor_dt[,correlation := cor(obs,prec_mean), by = .(lon,lat,system)]

for(sys in systems)
{
  pp = ggplot_dt(cor_dt[system == sys & year == year[1] & !(mask)],'correlation',
                 colorscale = scale_fill_gradient2(limits = c(-1,1),low = "blue", mid = "white", high = "red", 
                                                   name = 'cor'))
  pp = pp + ggtitle(paste0(sys,' correlation'))
  
  png(paste0(plot_dir,'prec_correlation_',sys,'.png'))
    print(pp)
  dev.off()
}


### get scores ###

prec_dt_masked = mask_precip(prec_dt,value_col = 'obs')[!(mask)][,mask := NULL]

mse_dt_pp = MSESS_ens_fc(prec_dt_masked,fc_cn = 'prec_pp',
                        obs_cn = 'obs',
                        by_cns = c('lon','lat','system'))

crps_dt_pp = CRPSS_ens_fc(prec_dt_masked,fc_cn = 'prec_pp',
                         obs_cn = 'obs',
                         by_cns = c('lon','lat','system'))





for(sys in systems)
{
  pp = ggplot_dt(mse_dt_pp[system == sys],'MSESS', rr= c(-1,1),
                 colorscale = scale_fill_gradient2(limits = c(-1,1),low = "blue", mid = "white", high = "red", 
                                                   name = 'MSES'))
  pp = pp + ggtitle(sys)
  
  png(paste0(plot_dir,'prec_msess_',sys,'.png'))
  print(pp)
  dev.off()
  
  
  pp = ggplot_dt(crps_dt_pp[system == sys],'CRPSS', rr= c(-1,1),
                 colorscale = scale_fill_gradient2(limits = c(-1,1),low = "blue", mid = "white", high = "red", 
                                                   name = 'CRPSS'))
  pp = pp + ggtitle(sys)
  
  png(paste0(plot_dir,'prec_crpss_',sys,'.png'))
  print(pp)
  dev.off()
}



