rm(list = ls())

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()

### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo')


init_mons = 8

data_dir = '/nr/project/stat/CONFER/Data/'
plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/Precip_skill/'

#################################

### get predictions:

prec_dt = fread(paste0(data_dir,'monthly_mean_prec/monthly_prec_fc_init_mon',init_mons,'.csv'))

prec_dt = prec_dt[forecast_year < 2021]
prec_dt = prec_dt[target_month %in% c(10:12)][,.(prec = mean(total_precipitation)),by = .(lon,lat,system,forecast_year,member)]
prec_dt[,season := 'OND']
setnames(prec_dt,'forecast_year','year')

prec_dt[,prec:= 3600 * 24*1000*prec]


### get observations:

obs_dt = load_chirps()
obs_dt[,index_coarse := NULL]

obs_dt = obs_dt[month %in% 10:12, .(obs = mean(prec)),by = .(lon,lat,year)][,season := 'OND']


prec_dt = merge(prec_dt,obs_dt,by = c('lon','lat','season','year'))

ggplot_dt(prec_dt[year == 1993 & member == 1 & system == 'cmcc'],'prec')
ggplot_dt(prec_dt[ member == 1 & system == 'cmcc',.(mean(obs)),.(lon,lat)])

mean_precs = prec_dt[,.(prec = mean(prec),
                        obs = mean(obs)),
                     by = .(season,year,system, member)]

# get IOD predictions and observations:
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
                          bycols = c('month','system'))

SST_dt[,IOD := IOD_west - IOD_east]

SST_dt = SST_dt[,.(IOD = mean(IOD), N34 = mean(N34)),by = .(system,year,member)]

#
mean_precs = merge(mean_precs,SST_dt,c('year','system','member'))

IOD_dt = IOD_dt[month %in% 10:12][,.(IOD_obs = mean(IOD)),year]
mean_precs = merge(mean_precs,IOD_dt,'year')

N34_dt = N34_dt[month %in% 10:12][,.(N34_obs = mean(N34)),year]
mean_precs = merge(mean_precs,N34_dt,'year')

### get correllations ###

mmean_precs = mean_precs[,lapply(.SD,mean),.SDcols = 5:10,by = .(year,system)]

mmean_precs[,cor(prec,obs),system]

mmean_precs[,cor(prec,IOD),system]
mmean_precs[,cor(obs,IOD_obs),system]
mmean_precs[,cor(obs,IOD),system]
