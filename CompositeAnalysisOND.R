### Assess the precipitation prediction skill for the different systems 
# In this script we look at seasonal mean precipitation over GHA and assess the skill of the Copernicus models before and after simple postprocessing

rm(list = ls())

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

library(data.table)
library(tidyverse)
devtools::load_all()

### which systems do you want to have a look at? ###

systems = c('cmcc','dwd','ecmwf','meteo_france','ukmo','mixed')

data_dir = '/nr/project/stat/CONFER/Data/monthly_mean_prec/'
plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/CompositeAnalysis/'

#################################

# get data

prec_dt = fread(paste0(data_dir,'OND_precip.csv'))

IOD_dt = fread('/nr/project/stat/CONFER/Data/IOD.csv')
  # reduce to OND mean:
  IOD_dt = IOD_dt[month %in% 10:12]
  IOD_dt = IOD_dt[,lapply(.SD,mean),.SDcols = c('IOD_west','IOD_east'),by = year]
  # standardize and get difference:
  IOD_dt[,IOD_west := (IOD_west - mean(IOD_west))/sd(IOD_west)]
  IOD_dt[,IOD_east := (IOD_east - mean(IOD_east))/sd(IOD_east)]
  IOD_dt[,IOD := IOD_west - IOD_east]

N34_dt = fread('/nr/project/stat/CONFER/Data/N34.csv')
  # reduce to OND mean and standardize  
  N34_dt = N34_dt[month %in% 10:12,.(N34 = mean(N34)),year]
  N34_dt[,N34 := (N34 - mean(N34))]


# get and process TC forecasts:  
tcfc_dt = fread('/nr/project/stat/CONFER/Data/N34_IOD_forecasts_raw.csv')
  
  # extract OND mean
  setnames(tcfc_dt,c('forecast_year','target_month'),c('year','month'))
  tcfc_dt = tcfc_dt[forecast_month == 8 & month %in% 10:12][,lapply(.SD,mean),.SDcols = c('N34','IOD_west','IOD_east'),by = .(system,year,member)]
  
  # standardize (makes conversion to Celsius superflous)
  tcfc_dt = tcfc_dt[,c(.(year = year,member = member),lapply(.SD,FUN = function(x) (x-mean(x))/sd(x))),.SDcols = c('N34','IOD_west','IOD_east'),by = system]
  
  tcfc_dt[,IOD:= IOD_west - IOD_east]
  
  # get a mixed ensemble forecast, for which each member is simply the ensemble mean for one of the forecast systems
  temp_dt  = tcfc_dt[,.(N34 = mean(N34), IOD = mean(IOD)),by = .(system,year)]
  temp_dt[,member := 1:.N,by = .(year)][,system := 'mixed']
  tcfc_dt = rbindlist(list(tcfc_dt[,.SD,.SDcols = !c('IOD_west','IOD_east')],temp_dt),use.names = T)
  setkey(tcfc_dt,system,year,member)

###### composite analysis based on OND averaged teleconnections ######
  
plot_dir_temp = paste0(plot_dir,'OND_TC/')
dir.create(plot_dir_temp,showWarnings = F)

# composite analysis for observations:

obs_dt = unique(prec_dt[,.(lon,lat,year,obs,mask)])

N34_ca_obs = composite_analysis(obs_dt[(!mask)],N34_dt,TC_name = 'N34',by_cols = c('lon','lat'),var_name = 'obs')
IOD_ca_obs = composite_analysis(obs_dt[(!mask)],IOD_dt,TC_name = 'IOD',by_cols = c('lon','lat'),var_name = 'obs')

N34_ca_pred = composite_analysis(prec_dt[(!mask)],tcfc_dt,TC_name = 'N34',by_cols = c('lon','lat','system'),var_name = 'prec',average_along_cols = c('year','member'))
IOD_ca_pred = composite_analysis(prec_dt[(!mask)],tcfc_dt,TC_name = 'IOD',by_cols = c('lon','lat','system'),var_name = 'prec',average_along_cols = c('year','member'))

#### plots: ####

#obs:

#N3.4
pp = ggplot_dt(N34_ca_obs,'x_plus',rr = c(-1.5,1.5)) + ggtitle('x_plus, N3.4, observed')
png(paste0(plot_dir_temp,'N34_x_plus_obs.png'))
  print(pp)
dev.off()

pp = ggplot_dt(N34_ca_obs,'x_minus',rr = c(-1.5,1.5)) + ggtitle('x_minus, N3.4, observed')
png(paste0(plot_dir_temp,'N34_x_minus_obs.png'))
  print(pp)
dev.off()

#IOD
pp = ggplot_dt(IOD_ca_obs,'x_plus',rr = c(-1.5,1.5)) + ggtitle('x_plus, IOD, observed')
png(paste0(plot_dir_temp,'IOD_x_plus_obs.png'))
  print(pp)
dev.off()

pp = ggplot_dt(IOD_ca_obs,'x_minus',rr = c(-1.5,1.5)) + ggtitle('x_minus, IOD, observed')
png(paste0(plot_dir_temp,'IOD_x_minus_obs.png'))
  print(pp)
dev.off()

#systems:

for(sys in systems)
{
  #N3.4
  pp = ggplot_dt(N34_ca_pred[system == sys],'x_plus',rr = c(-1.5,1.5)) + ggtitle(paste0('x_plus, N3.4, ',sys))
  png(paste0(plot_dir_temp,'N34_x_plus_',sys,'.png'))
  print(pp)
  dev.off()
  
  pp = ggplot_dt(N34_ca_pred[system == sys],'x_minus',rr = c(-1.5,1.5)) + ggtitle(paste0('x_minus, N3.4, ',sys))
  png(paste0(plot_dir_temp,'N34_x_minus_',sys,'.png'))
  print(pp)
  dev.off()
  
  #IOD
  pp = ggplot_dt(IOD_ca_pred[system == sys],'x_plus',rr = c(-1.5,1.5)) + ggtitle(paste0('x_plus, IOD, ',sys))
  png(paste0(plot_dir_temp,'IOD_x_plus_',sys,'.png'))
  print(pp)
  dev.off()
  
  pp = ggplot_dt(IOD_ca_pred[system == sys],'x_minus',rr = c(-1.5,1.5)) + ggtitle(paste0('x_minus, IOD, ',sys))
  png(paste0(plot_dir_temp,'IOD_x_minus_',sys,'.png'))
  print(pp)
  dev.off()
}

##########################################################################################
###### Impact of Teleconnection state in August (time of issue) on prediction error ######
##########################################################################################

plot_dir_temp = paste0(plot_dir,'ErrorTCConnection/')
dir.create(plot_dir_temp,showWarnings = F)

# get teleconnection state in August:

IOD_dt = fread('/nr/project/stat/CONFER/Data/IOD.csv')
  # reduce to OND mean:
  IOD_dt = IOD_dt[month %in% 8]
  IOD_dt = IOD_dt[,lapply(.SD,mean),.SDcols = c('IOD_west','IOD_east'),by = year]
  # standardize and get difference:
  IOD_dt[,IOD_west := (IOD_west - mean(IOD_west))/sd(IOD_west)]
  IOD_dt[,IOD_east := (IOD_east - mean(IOD_east))/sd(IOD_east)]
  IOD_dt[,IOD := IOD_west - IOD_east]

N34_dt = fread('/nr/project/stat/CONFER/Data/N34.csv')
  N34_dt = N34_dt[month == 8,.(N34 = mean(N34)),year]
  N34_dt[,N34 := (N34 - mean(N34))/sd(N34)]


##############

mean_pred = prec_dt[,.(prec = mean(prec), 
                       obs = mean(obs)),
                    by = .(lon,lat,year,system)]

mean_pred = merge(mean_pred,N34_dt,by = 'year')
mean_pred = merge(mean_pred,IOD_dt[,.(year,IOD)],by = 'year')

mean_pred[,error := prec - obs]
  
for(sys in systems)
{
  temp_dt = mean_pred[system == sys,cor(error,N34),.(lon,lat)]
  setnames(temp_dt,'V1','correlation')
  pp = ggplot_dt(temp_dt,rr = c(-1,1))+ ggtitle(paste0('correlation of ',sys,' error to N3.4'))
  
  png(paste0(plot_dir_temp,'error_cor_N34_',sys,'.png'))
    print(pp) 
  dev.off()
  
  temp_dt = mean_pred[system == sys,cor(error,IOD),.(lon,lat)]
  setnames(temp_dt,'V1','correlation')
  pp = ggplot_dt(temp_dt,rr = c(-1,1))+ ggtitle(paste0('correlation of ',sys,' error to IOD'))
  
  png(paste0(plot_dir_temp,'error_cor_IOD_',sys,'.png'))
  print(pp) 
  dev.off()
}
  
