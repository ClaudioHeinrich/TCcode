# here we get a data table for monthly mean precipitation, containing CHIRPS observations and predictions from all copernicus models and the SFE prediction

rm(list = ls())


library(data.table)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

library(PostProcessing)

devtools::load_all()
devtools::document()


### which systems do you want to have a look at? ###

data_dir = '/nr/project/stat/CONFER/Data/'


load(paste0(data_dir,'ICPAC_region.RData'))

CHIRPS_dt = fread(paste0(data_dir,'CHIRPS_prec_upscaled.csv'))

for(lag in 0:2)
{
  
  print(paste0('lag = ',lag))
#plot_dir = '/nr/www/virtual/files.nr.no/htdocs/samba/CONFER/figures/'
plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/'
dir.create(plot_dir,showWarnings = F)

  ## Nino 3.4 ##
  
  EN_dt = fread(paste0(data_dir,'N34.csv'))
  EN_dt = time_lag_monthly_dt(EN_dt,lag = lag)
  
  ### composite analysis ###
  
  CA_dt = composite_analysis(var_dt = CHIRPS_dt,
                             TC_dt = EN_dt,
                             TC_name = 'N34',
                             by_cols = c('month','lon','lat'),
                             var_name = 'prec')
  
  

  #plot_dir_temp = paste0(plot_dir,'/CompositeAnalysis/obs/lag',lag,'/N34/')
  #dir.create(plot_dir_temp,recursive = T,showWarnings = F)
  
  rr = c(-4,4)
  for(mm in 10:12)
  {
    pp = ggplot_dt(CA_dt[month == mm],data_col = 'x_plus',mn = paste0('Nino 3.4 x_plus, month ',mm,', lag ',lag),rr=rr,
                   colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_plus',limits = rr))
    
    png(file = paste0(plot_dir,'N34_x_plus_m',mm,'_lag',lag,'.png'))
    print(pp)
    dev.off()
    
    pp = ggplot_dt(CA_dt[month == mm],data_col = 'x_minus',mn = paste0('Nino 3.4 x_minus, month ',mm,', lag ',lag),rr=rr,
                   colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_minus',limits = rr))
    
    png(file = paste0(plot_dir,'N34_x_minus_m',mm,'_lag',lag,'.png'))
    print(pp)
    dev.off()
  }
  
  ### correlation plots ###
  
  # plot_dir_temp = paste0(plot_dir,'/correlations/obs/lag',lag,'/N34/')
  # dir.create(plot_dir_temp,recursive = T,showWarnings = F)
  # 
  rr = c(-0.7,0.7)
  dt_temp = merge(CHIRPS_dt,EN_dt, by = c('year','month'))
  dt_temp = mask_precip(dt_temp)
  cor_dt = dt_temp[!(mask),cor(prec,N34),by = .(month,lon,lat)]
  setnames(cor_dt,'V1','correlation')
  
  for(mm in 10:12)
  {
    pp = ggplot_dt(cor_dt[month == mm],data_col = 'correlation',mn = paste0('Nino 3.4 correlation, month ',mm,', lag ',lag),rr=rr,
                  colorscale = scale_fill_gradient2(low = "blue", mid = "white", high = "red",name = 'correlation',limits = rr))
    
    png(file = paste0(plot_dir,'N34_correlation_m',mm,'_lag',lag,'.png'))
    print(pp)
    dev.off()
  }
   
  #########################################
  
  ## IOD ##
  
  IOD_dt = fread(paste0(data_dir,'IOD.csv'))
  IOD_dt = time_lag_monthly_dt(IOD_dt,lag = lag)
  
  IOD_dt[,IOD_west := (IOD_west - mean(IOD_west))/sd(IOD_west),by = month]
  IOD_dt[,IOD_east := (IOD_east - mean(IOD_east))/sd(IOD_east),by = month]
  IOD_dt[,IOD:= IOD_west - IOD_east]
  
  CA_dt = composite_analysis(var_dt = CHIRPS_dt,
                             TC_dt = IOD_dt,
                             TC_name = 'IOD',
                             by_cols = c('month','lon','lat'),
                             var_name = 'prec')
  
  
  plot_dir_temp = paste0(plot_dir,'/CompositeAnalysis/obs/lag',lag,'/IOD/')
  dir.create(plot_dir_temp,recursive = T,showWarnings = F)
  
  
  rr = c(-4,4)
  for(mm in 10:12)
  {
    pp = ggplot_dt(CA_dt[month == mm],data_col = 'x_plus',mn = paste0('IOD x_plus, month ',mm,', lag ',lag),rr=rr,
                   colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_plus',limits = rr))
    
    png(file = paste0(plot_dir,'N34_x_plus_m',mm,'_lag',lag,'.png'))
    print(pp)
    dev.off()
    
    pp = ggplot_dt(CA_dt[month == mm],data_col = 'x_minus',mn = paste0('IOD x_minus, month ',mm,', lag ',lag),rr=rr,
                   colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_minus',limits = rr))
    
    png(file = paste0(plot_dir,'N34_x_minus_m',mm,'_lag',lag,'.png'))
    print(pp)
    dev.off()
  }
  
  ### correlations ###
  
  plot_dir_temp = paste0(plot_dir,'/correlations/obs/lag',lag,'/IOD/')
  dir.create(plot_dir_temp,recursive = T,showWarnings = F)
  
  rr = c(-0.7,0.7)
  dt_temp = merge(CHIRPS_dt,IOD_dt, by = c('year','month'))
  dt_temp = mask_precip(dt_temp)
  cor_dt = dt_temp[!(mask),cor(prec,IOD),by = .(month,lon,lat)]
  setnames(cor_dt,'V1','correlation')
  
  for(mm in 10:12)
  {
    pp = ggplot_dt(cor_dt[month == mm],data_col = 'correlation',mn = paste0('IOD correlation, month ',mm,', lag ',lag),rr=rr,
                   colorscale = scale_fill_gradient2(low = "blue", mid = "white", high = "red",name = 'correlation',limits = rr))
    
    png(file = paste0(plot_dir,'IOD_correlation_m',mm,'_lag',lag,'.png'))
    print(pp)
    dev.off()
  }
  

#   ############################################
#   
#   ### Somali jet ###
#   
#   SJ_dt = fread(paste0(data_dir,'SJ_obs.csv'))
#   SJ_dt = time_lag_monthly_dt(SJ_dt,lag = lag)
#   CA_dt = composite_analysis(var_dt = CHIRPS_dt,
#                              TC_dt = SJ_dt,
#                              TC_name = 'SJI',
#                              by_cols = c('month','lon','lat'),
#                              var_name = 'prec')
#   
#   
#   plot_dir_temp = paste0(plot_dir,'/CompositeAnalysis/obs/lag',lag,'/SJ/')
#   dir.create(plot_dir_temp,recursive = T,showWarnings = F)
#   
#   
#   rr = c(-4,4)
#   for(mm in 1:12)
#   {
#     pp = ggplot_dt(CA_dt[month == mm],data_col = 'x_plus',mn = paste0('SJ x_plus, month ',mm),rr=rr,
#                    colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_plus',limits = rr))
#     
#     png(file = paste0(plot_dir_temp,'x_plus_m',mm,'.png'))
#     print(pp)
#     dev.off()
#     
#     pp = ggplot_dt(CA_dt[month == mm],data_col = 'x_minus',mn = paste0('SJ x_minus, month ',mm),rr=rr,
#                    colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_minus',limits = rr))
#     
#     png(file = paste0(plot_dir_temp,'x_minus_m',mm,'.png'))
#     print(pp)
#     dev.off()
#   }
# 
#   ### correlations ###
#   
#   plot_dir_temp = paste0(plot_dir,'/correlations/obs/lag',lag,'/SJ/')
#   dir.create(plot_dir_temp,recursive = T,showWarnings = F)
#   
#   rr = c(-0.7,0.7)
#   dt_temp = merge(CHIRPS_dt,SJ_dt, by = c('year','month'))
#   dt_temp = mask_precip(dt_temp)
#   cor_dt = dt_temp[!(mask),cor(prec,SJI),by = .(month,lon,lat)]
#   setnames(cor_dt,'V1','correlation')
#   
#   for(mm in 1:12)
#   {
#     pp = ggplot_dt(cor_dt[month == mm],data_col = 'correlation',mn = paste0('Somali Jet Index correlation, month ',mm),rr=rr,
#                    colorscale = scale_fill_gradient2(low = "blue", mid = "white", high = "red",name = 'correlation',limits = rr))
#     
#     png(file = paste0(plot_dir_temp,'correlation_m',mm,'.png'))
#     print(pp)
#     dev.off()
#   }
#   
#   
#   
}
