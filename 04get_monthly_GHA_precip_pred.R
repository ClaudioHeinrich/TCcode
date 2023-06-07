rm(list = ls())


library(data.table)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()


### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo')
lead_times = 1:5


fc_dir = '/nr/samba/user/claudio/bigdisk/CONFER/Systems_monthly_compiled/'
data_dir = '/nr/project/stat/CONFER/Data/monthly_mean_prec/'

dir.create(data_dir,showWarnings = F)


gha = GHA_locs('half_deg')

gha_lons = unique(gha[,lon])
gha_lats = unique(gha[,lat])
for(mm in 1:12)
{
  cat('month =',mm)
  temp_dt = load_cds_monthly_data(variable = 'total_precipitation',
                                 forecast_month = mm,
                                 systems = systems,
                                 lead_times = lead_times,
                                 lon_subset = gha_lons,
                                 lat_subset = gha_lats,
                                 root_dir = claudio_sfe_dir())

  temp_dt = merge(temp_dt,gha,by = c('lon','lat'))

  gc()

  fwrite(temp_dt,file = paste0(data_dir,'monthly_prec_fc_init_mon',mm,'.csv'))

}

