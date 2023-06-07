# derive several teleconnection-relevant indices based on SST-averaging, both on ERA5 and system forecasts

rm(list = ls())


library(data.table)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()


### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo')

# where do you want to store stuff?
data_dir = '/nr/project/stat/CONFER/Data/'

### get the Nino 3,4 index ###

obs_dt = load_era_monthly_data(variable = 'sea_surface_temperature',months = 1:12,years = 1970:2021,root_dir = claudio_sfe_dir())

setnames(obs_dt,5,'sst')

obs_dt[,sst := sst - 273.15]

obs_dt = obs_dt[!is.na(sst)]

### get Nino 3.4
# see https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni

N34_dt = obs_dt[lon %between% c(-170,-120) & lat %between% c(-5,5),mean(sst),by = .(year,month)]
setnames(N34_dt,'V1','N34')
fwrite(N34_dt,file = paste0(data_dir,'N34.csv'))

### get IOD

temp_dt1 = obs_dt[lon %between% c(50,70) & lat %between% c(-10,10),mean(sst),by = .(year,month)]
temp_dt2 = obs_dt[lon %between% c(90,110) & lat %between% c(-10,0),mean(sst),by = .(year,month)]

temp_dt1[,IOD_east := temp_dt2[,V1]]
setnames(temp_dt1,'V1','IOD_west')

fwrite(temp_dt1,file = paste0(data_dir,'IOD.csv'))

