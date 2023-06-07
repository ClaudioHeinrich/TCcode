

### plot correlation of ON precip and SSTs ###

rm(list = ls())

devtools::load_all()
library(data.table)
library(ggplot2)

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


# recompute (oos) climatology and clim_sd:

prec_dt = loyo(prec_dt,FUN = mean,SDcols = 'obs',bycols = c('lon','lat'))
prec_dt[,clim := NULL]
setnames(prec_dt,'obs_new','clim')

prec_dt = loyo(prec_dt,FUN = sd,SDcols = 'obs',bycols = c('lon','lat','system','member'))
prec_dt[,clim_sd := NULL]
setnames(prec_dt,'obs_new','clim_sd')

########

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




### regress ON precip ###
obs_dt = fread(paste0(data_dir,'CHIRPS_prec_upscaled.csv'))
obs_dt = obs_dt[month %in% 10:11,.(prec = mean(prec)),.(year,lon,lat)]

IOD_Aug = IOD_dt[month == 8,.(year,IOD)]
IOD_Sep = IOD_dt[month == 9,.(year,IOD)]
IOD_Oct = IOD_dt[month == 10,.(year,IOD)]
IOD_Nov = IOD_dt[month == 11,.(year,IOD)]

setnames(IOD_Aug,'IOD','IOD_Aug')
setnames(IOD_Sep,'IOD','IOD_Sep')
setnames(IOD_Oct,'IOD','IOD_Oct')
setnames(IOD_Nov,'IOD','IOD_Nov')

N34_Aug = N34_dt[month == 8,.(year,N34)]
N34_Sep = N34_dt[month == 9,.(year,N34)]
N34_Oct = N34_dt[month == 10,.(year,N34)]
N34_Nov = N34_dt[month == 11,.(year,N34)]

setnames(N34_Aug,'N34','N34_Aug')
setnames(N34_Sep,'N34','N34_Sep')
setnames(N34_Oct,'N34','N34_Oct')
setnames(N34_Nov,'N34','N34_Nov')

obs_dt = merge(obs_dt,IOD_Aug,by = 'year')
obs_dt = merge(obs_dt,IOD_Sep,by = 'year')
obs_dt = merge(obs_dt,IOD_Oct,by = 'year')
obs_dt = merge(obs_dt,IOD_Nov,by = 'year')

obs_dt = merge(obs_dt,N34_Aug,by = 'year')
obs_dt = merge(obs_dt,N34_Sep,by = 'year')
obs_dt = merge(obs_dt,N34_Oct,by = 'year')
obs_dt = merge(obs_dt,N34_Nov,by = 'year')

# reduce to correct region:
obs_dt = obs_dt[!is.na(prec)]
setkey(obs_dt,lon,lat)
obs_dt = obs_dt[GHA_locs()]

obs_dt = mask_precip(obs_dt)[(!mask)][,mask := NULL]
#### plots ####

plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/correlations/'

theme_set(theme_bw(base_size = 20))

fullmon = c('August','September','October','November')
for(mon in c('Aug','Sep','Oct','Nov'))
{
  iod_cor_dt = obs_dt[,cor(prec,get(paste0('IOD_',mon))),by = .(lon,lat)]
  setnames(iod_cor_dt,'V1','cor')
  
  fm = fullmon[which(mon == c('Aug','Sep','Oct','Nov'))]
  
  pp = ggplot_dt(iod_cor_dt,rr = c(-1,1)) + ggtitle(paste0('IOD ',fm,' correlation'))
  png(paste0(plot_dir,'IOD_',mon,'_ON_prec_cor.png'))
    print(pp)
  dev.off()
  
  n34_cor_dt = obs_dt[,cor(prec,get(paste0('N34_',mon))),by = .(lon,lat)]
  setnames(n34_cor_dt,'V1','cor')
  
  pp = ggplot_dt(n34_cor_dt,rr = c(-1,1)) + ggtitle(paste0('N34 ',fm,' correlation'))
  png(paste0(plot_dir,'N34_',mon,'_ON_prec_cor.png'))
  print(pp)
  dev.off()
  
  
}



