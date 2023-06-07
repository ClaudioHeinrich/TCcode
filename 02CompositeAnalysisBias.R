# here we get a data table for monthly mean precipitation, containing CHIRPS observations and predictions from all copernicus models and the SFE prediction

rm(list = ls())


library(data.table)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

#library(PostProcessing)

devtools::load_all()
devtools::document()


### which systems do you want to have a look at? ###

data_dir = '/nr/project/stat/CONFER/Data/'

plot_dir = '/nr/project/stat/CONFER/plots/monthlyTCs/'
load(paste0(data_dir,'ICPAC_region.RData'))

CHIRPS_dt = fread(paste0(data_dir,'CHIRPS_prec_upscaled.csv'))
setnames(CHIRPS_dt,'prec','obs')

lead_times = 1:4

for(system in setdiff(global_systems(),c('ukmo','ecmwf')))
{
  print(paste0('###### system = ',system,' ######'))

  for(lag in 0:2)
  {
    
    print(paste0('### lag = ',lag,' ###'))
    #plot_dir = '/nr/www/virtual/files.nr.no/htdocs/samba/CONFER/figures/'
    
    
    for(fm in 1:12)
    {
      
      print(paste0('month = ',fm))
      # get forecast
      
      timer = Sys.time()
      
      fc_dt = load_cds_monthly_data(variable = 'total_precipitation',
                                    forecast_month = fm,
                                    years = NULL,
                                    systems = system,
                                    root_dir = claudio_sfe_dir(),
                                    lon_subset = global_gha_lon_subset(),
                                    lat_subset = global_gha_lat_subset(),
                                    lead_times = lead_times
                                    )
      
      timer = Sys.time() - timer
      
      # reduce to ensemble mean and minor adjustments:
      
      fc_dt = take_ens_mean(fc_dt,pred_vars = 'total_precipitation')
      
      setnames(fc_dt,'total_precipitation','prec')
      #convert m/s to mm/day:
      fc_dt[,prec := prec * 1000 * 3600 *24]
      setnames(fc_dt,c('target_month'),'month')
      # get year, note that it's not necessarily the same as forecast_year, which is when the forecast is initialized
      fc_dt[,year:= forecast_year]
      fc_dt[forecast_month > month,year := year + 1]
      fc_dt[,lead_time := month - forecast_month]
      fc_dt[lead_time <= 0, lead_time := lead_time + 12]
      fc_dt[,c('forecast_year','forecast_month') := NULL]
      
      
      # merge with observations:
      
      
      fc_dt = merge(fc_dt,CHIRPS_dt,by= c('year','month','lon','lat'))
      fc_dt[,bias := prec - obs]
      
      ################ Composite analysis and correlations of bias ##################
      
      
      for(lt in lead_times)
      {
        print(paste0('lead time = ',lt))
        
        
        ############################
        ## Nino 3.4 ##
        
        #print('N34')
        
        EN_dt = fread(paste0(data_dir,'N34.csv'))
        EN_dt = time_lag_monthly_dt(EN_dt,lag = lag + lt)
        
        ### composite analysis ###
        
        CA_dt = composite_analysis(var_dt = fc_dt[lead_time == lt],
                                   TC_dt = EN_dt,
                                   TC_name = 'N34',
                                   by_cols = c('month','lon','lat','lead_time'),
                                   var_name = 'bias')
          
        mm = CA_dt[,unique(month)] # target month
      
        plot_dir_temp = paste0(plot_dir,'/CompositeAnalysis/',system,'/lead_time',lt,'/lag',lag,'/N34/')
        dir.create(plot_dir_temp,recursive = T,showWarnings = F)
        
        rr = c(-4,4)
        pp = ggplot_dt(CA_dt[month == mm & lead_time == lt],data_col = 'x_plus',mn = paste0('Nino 3.4 x_plus, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_plus',limits = rr))
        
        png(file = paste0(plot_dir_temp,'x_plus_m',mm,'.png'))
        print(pp)
        dev.off()
        
        pp = ggplot_dt(CA_dt[month == mm & lead_time == lt],data_col = 'x_minus',mn = paste0('Nino 3.4 x_minus, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_minus',limits = rr))
        
        png(file = paste0(plot_dir_temp,'x_minus_m',mm,'.png'))
        print(pp)
        dev.off()
      
        ### correlation plots ###
        
        plot_dir_temp = paste0(plot_dir,'/correlations/',system,'/lead_time',lt,'/lag',lag,'/N34/')
        dir.create(plot_dir_temp,recursive = T,showWarnings = F)
        
        rr = c(-0.7,0.7)
        dt_temp = copy(fc_dt[lead_time == lt])
        dt_temp = merge(dt_temp,EN_dt,by = c('year','month'))
        dt_temp = mask_precip(dt_temp,value_col = 'obs')
        cor_dt = dt_temp[!(mask),cor(bias,N34),by = .(month,lon,lat)]
        setnames(cor_dt,'V1','correlation')
        
        pp = ggplot_dt(cor_dt[month == mm],data_col = 'correlation',mn = paste0('Nino 3.4 correlation, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "blue", mid = "white", high = "red",name = 'correlation',limits = rr))
        
        png(file = paste0(plot_dir_temp,'correlation_m',mm,'.png'))
        print(pp)
        dev.off()
        
        ####################################
        
        #### IOD ####
        
        #print('IOD')
        
        IOD_dt = fread(paste0(data_dir,'IOD.csv'))
        IOD_dt = time_lag_monthly_dt(IOD_dt,lag = lag + lt)
        
        ### composite analysis ###
        
        CA_dt = composite_analysis(var_dt = fc_dt[lead_time == lt],
                                   TC_dt = IOD_dt,
                                   TC_name = 'IOD',
                                   by_cols = c('month','lon','lat','lead_time'),
                                   var_name = 'bias')
        
        mm = CA_dt[,unique(month)] # target month
        
        plot_dir_temp = paste0(plot_dir,'/CompositeAnalysis/',system,'/lead_time',lt,'/lag',lag,'/IOD/')
        dir.create(plot_dir_temp,recursive = T,showWarnings = F)
        
        rr = c(-4,4)
        pp = ggplot_dt(CA_dt[month == mm & lead_time == lt],data_col = 'x_plus',mn = paste0('IOD x_plus, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_plus',limits = rr))
        
        png(file = paste0(plot_dir_temp,'x_plus_m',mm,'.png'))
          print(pp)
        dev.off()
        
        pp = ggplot_dt(CA_dt[month == mm & lead_time == lt],data_col = 'x_minus',mn = paste0('IOD x_minus, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_minus',limits = rr))
        
        png(file = paste0(plot_dir_temp,'x_minus_m',mm,'.png'))
          print(pp)
        dev.off()
        
        ### correlation plots ###
        
        plot_dir_temp = paste0(plot_dir,'/correlations/',system,'/lead_time',lt,'/lag',lag,'/IOD/')
        dir.create(plot_dir_temp,recursive = T,showWarnings = F)
        
        rr = c(-0.7,0.7)
        dt_temp = copy(fc_dt[lead_time == lt])
        dt_temp = merge(dt_temp,IOD_dt,by = c('year','month'))
        dt_temp = mask_precip(dt_temp,value_col = 'obs')
        cor_dt = dt_temp[!(mask),cor(bias,IOD),by = .(month,lon,lat)]
        setnames(cor_dt,'V1','correlation')
        
        pp = ggplot_dt(cor_dt[month == mm],data_col = 'correlation',mn = paste0('IOD correlation, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "blue", mid = "white", high = "red",name = 'correlation',limits = rr))
        
        png(file = paste0(plot_dir_temp,'correlation_m',mm,'.png'))
          print(pp)
        dev.off()
        
        #########################################
        ### Somali Jet ###
        
        #print('Somali Jet')
        
        SJ_dt = fread(paste0(data_dir,'SJ_obs.csv'))
        SJ_dt = time_lag_monthly_dt(SJ_dt,lag = lag + lt)
        
        ### composite analysis ###
        
        CA_dt = composite_analysis(var_dt = fc_dt[lead_time == lt],
                                   TC_dt = SJ_dt,
                                   TC_name = 'SJI',
                                   by_cols = c('month','lon','lat','lead_time'),
                                   var_name = 'bias')
        
        mm = CA_dt[,unique(month)] # target month
        
        plot_dir_temp = paste0(plot_dir,'/CompositeAnalysis/',system,'/lead_time',lt,'/lag',lag,'/SJ/')
        dir.create(plot_dir_temp,recursive = T,showWarnings = F)
        
        rr = c(-4,4)
        pp = ggplot_dt(CA_dt[month == mm & lead_time == lt],data_col = 'x_plus',mn = paste0('SJI x_plus, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_plus',limits = rr))
        
        png(file = paste0(plot_dir_temp,'x_plus_m',mm,'.png'))
        print(pp)
        dev.off()
        
        pp = ggplot_dt(CA_dt[month == mm & lead_time == lt],data_col = 'x_minus',mn = paste0('SJI x_minus, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "red", mid = "white", high = "blue",name = 'x_minus',limits = rr))
        
        png(file = paste0(plot_dir_temp,'x_minus_m',mm,'.png'))
        print(pp)
        dev.off()
        
        ### correlation plots ###
        
        plot_dir_temp = paste0(plot_dir,'/correlations/',system,'/lead_time',lt,'/lag',lag,'/SJ/')
        dir.create(plot_dir_temp,recursive = T,showWarnings = F)
        
        rr = c(-0.7,0.7)
        dt_temp = copy(fc_dt[lead_time == lt])
        dt_temp = merge(dt_temp,SJ_dt,by = c('year','month'))
        dt_temp = mask_precip(dt_temp,value_col = 'obs')
        cor_dt = dt_temp[!(mask),cor(bias,SJI),by = .(month,lon,lat)]
        setnames(cor_dt,'V1','correlation')
        
        pp = ggplot_dt(cor_dt[month == mm],data_col = 'correlation',mn = paste0('SJI correlation, month ',mm),rr=rr,
                       colorscale = scale_fill_gradient2(low = "blue", mid = "white", high = "red",name = 'correlation',limits = rr))
        
        png(file = paste0(plot_dir_temp,'correlation_m',mm,'.png'))
        print(pp)
        dev.off()
        
        
      }
    }
  }
}
      #########################################