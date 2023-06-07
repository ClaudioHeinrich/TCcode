# here we get a data table for monthly mean precipitation, containing CHIRPS observations and predictions from all copernicus models and the SFE prediction

rm(list = ls())


library(data.table)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

library(PostProcessing)

devtools::load_all()
devtools::document()


### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ncep','ukmo')

target_mons = 10:12
init_mon = 8

fc_dir = '/nr/samba/user/claudio/bigdisk/CONFER/Systems_monthly_compiled/'
SFE_dir = '/nr/user/claudio/bigdisk/SFE/Forecasts/'
data_dir = '/nr/project/stat/CONFER/Data/'

EN_dt = fread(paste0(data_dir,'N3.4.csv'))
load(paste0(data_dir,'ICPAC_region.RData'))

CHIRPS_dt = fread(paste0(data_dir,'CHIRPS_prec_upscaled.csv'))

########################
# 
# 
# lead_times = 0
# 
# prec_dt =load_cds_monthly_data(variable = 'total_precipitation',
#                               forecast_month = init_mon,
#                               years = NULL,
#                               pressure = NA,
#                               systems = global_systems(),
#                               lead_times = lead_times,
#                               lon_subset = global_confer_lon_subset(),
#                               lat_subset = global_confer_lat_subset(),
#                               root_dir = claudio_sfe_dir())

# function for getting precipictation for locations covered by ICPAC predictions:

merge_nwp_precip = function(target_mons,init_mons, systems = 'ecmwf')
{
  
  for(sys in systems)
  {
    print(sys)
    temp = load_dt_monthly('prec','ICPAC',sys,target_mons,init_mons)
    temp = load_dt_monthly('prec','ICPAC',sys,target_mons,init_mons)
    temp[,(sys) := total_precipitation * 1000 * 24*3600]
    temp[,total_precipitation := NULL]
    assign(paste0(sys,'_dt'),
           value = temp)
  }
  
  
  ret_dt = get(paste0(systems[1],'_dt'))
  
  if(length(systems) > 1)
  {
    for(sys in systems[2:length(systems)])
    {
      
      ret_dt = merge(ret_dt,get(paste0(sys,'_dt')),by = c('lon','lat','year','month','member'),all = TRUE)
    }
  }
  
  return(ret_dt)
}


#################

#generate data table:
prec_dt = merge_nwp_precip(target_mons,init_mons,systems)

# add observations:
prec_dt_new = merge(prec_dt,CHIRPS_dt,by = c('year','month','lon','lat'))
prec_dt = prec_dt_new[,index_coarse := NULL]
setnames(prec_dt,'prec','obs')

### get SFE forecast ###

test = fread(paste0(SFE_dir,'forecast_total_precipitation_2005_9.csv'))

fwrite(prec_dt,file = paste0(data_dir,'GHA_OND_prec_08.csv'))

