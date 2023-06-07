rm(list = ls())


library(data.table)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()


### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo')


fc_dir = '/nr/samba/user/claudio/bigdisk/CONFER/Systems_monthly_compiled/'
data_dir = '/nr/project/stat/CONFER/Data/'

### get the Nino 3,4 index ###
# see https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni


obs_dt = load_era_monthly_data(variable = 'sea_surface_temperature',months = 1:12,years = 1970:2021,root_dir = claudio_sfe_dir())


## get the same for predictions ##

EN_dt = copy(obs_dt)

sys = systems[1]

print(sys)
sys_dt = load_dt_monthly(var_name = 'SST',
                         area = 'Nino3.4',
                         system_name = sys,
                         target_mons = target_mons,
                         init_mons = init_mons)

sys_dt[,clim := mean(sea_surface_temperature),.(lon,lat,month,member)]
sys_dt[,paste0(sys,'_ano') := sea_surface_temperature - clim,.(lon,lat,year,month,member)]
sys_dt = sys_dt[,mean(get(paste0(sys,'_ano'))),.(year,month,member)]
setnames(sys_dt,'V1',paste0(sys,'_ano'))
setkey(sys_dt,year,month,member)  

if(length(systems) > 1)
{
  
  for(sys in systems[2 : length(systems)])
  {
    print(sys)
    temp_dt = load_dt_monthly(var_name = 'SST',
                             area = 'Nino3.4',
                             system_name = sys,
                             target_mons = target_mons,
                             init_mons = init_mons)
    
    temp_dt[,clim := mean(sea_surface_temperature),.(lon,lat,month,member)]
    temp_dt[,paste0(sys,'_ano') := sea_surface_temperature - clim,.(lon,lat,year,month,member)]
    temp_dt = temp_dt[,mean(get(paste0(sys,'_ano'))),.(year,month,member)]
    setnames(temp_dt,'V1',paste0(sys,'_ano'))
    setkey(temp_dt,year,month,member)  
    
    sys_dt = merge(sys_dt,temp_dt,by = c('year','month','member'),all = TRUE)
  
  }
  
  EN_dt = merge(EN_dt,sys_dt,by = c('year','month'))
  
}  




### get errors and ensemble member ranks ###

for(sys in systems)
{
  EN_dt[,paste0(sys,'_err') := abs(get(paste0(sys,'_ano')) - ERA_ano)]
  EN_dt[,paste0(sys,'_rk') := rank(get(paste0(sys,'_err')),na.last = 'keep'),.(year,month)]
}

### save ###

fwrite(EN_dt,file = paste0(data_dir,'N3.4.csv'))

