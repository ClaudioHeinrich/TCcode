rm(list = ls())


library(data.table)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()


### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo')
lead_times = 1:5


fc_dir = '/nr/samba/user/claudio/bigdisk/CONFER/Systems_monthly_compiled/'
data_dir = '/nr/project/stat/CONFER/Data/'

### get the Nino 3,4 index ###
# see https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni


obs_dt = load_era_monthly_data(variable = 'sea_surface_temperature',months = 1:12,years = 1970:2021,root_dir = claudio_sfe_dir())


## get the same for predictions ##

N34_dt = data.table()
IOD_east_dt = data.table()
IOD_west_dt = data.table()

for(mm in 1:12)
{
  print(mm)
temp_dt = load_cds_monthly_data(variable = 'sea_surface_temperature',
                               forecast_month = mm,
                               systems = systems,
                               lead_times = lead_times,
                               root_dir = claudio_sfe_dir())

# extract teleconnection-relevant areas
N34_dt_temp = temp_dt[lon %between% c(-170,-120) & lat %between% c(-5,5),mean(sea_surface_temperature,na.rm = T),by = .(system, forecast_year,forecast_month,target_month,member)]
setnames(N34_dt_temp,'V1','N34')
N34_dt = rbindlist(list(N34_dt,N34_dt_temp))


#IOD:
IOD_west_dt_temp = temp_dt[lon %between% c(50,70) & lat %between% c(-10,10),mean(sea_surface_temperature,na.rm = T),by = .(system, forecast_year,forecast_month,target_month,member)]
IOD_east_dt_temp = temp_dt[lon %between% c(90,110) & lat %between% c(-10,0),mean(sea_surface_temperature,na.rm = T),by = .(system, forecast_year,forecast_month,target_month,member)]
setnames(IOD_west_dt_temp,'V1','IOD_west')
setnames(IOD_east_dt_temp,'V1','IOD_east')
IOD_west_dt = rbindlist(list(IOD_west_dt,IOD_west_dt_temp))
IOD_east_dt = rbindlist(list(IOD_east_dt,IOD_east_dt_temp))

}


# combine

TC_dt = merge(N34_dt,IOD_west_dt,by = c('system','forecast_year','forecast_month','target_month','member'))
TC_dt = merge(TC_dt,IOD_east_dt,by = c('system','forecast_year','forecast_month','target_month','member'))

### save ###

fwrite(TC_dt,file = paste0(data_dir,'N34_IOD_forecasts_raw.csv'))

