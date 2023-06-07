## Predictability study 

rm(list = ls())


library(data.table)
library(ggplot2)

setwd('/nr/samba/user/claudio/pkg/PostProcessing/')

devtools::load_all()


### which systems do you want to have a look at? ###

fc_dir = '/nr/samba/user/claudio/bigdisk/CONFER/Systems_monthly_compiled/'
data_dir = '/nr/project/stat/CONFER/Data/'

fc_dt = fread(paste0(data_dir,'N34_IOD_forecasts_raw.csv'))
systems = unique(pred_dt[,system])

# observations:

N34_dt = fread(paste0(data_dir,'N34.csv'))
IOD_dt = fread(paste0(data_dir,'IOD.csv'))

### process forecasts ###

#Kelvin to Celsius
fc_dt[,c('N34','IOD_west','IOD_east') := lapply(.SD,FUN = function(x){x - 273.15}),.SDcols = c('N34','IOD_west','IOD_east')]

# standardize everything:
fc_dt[,c('N34','IOD_west','IOD_east') := lapply(.SD,FUN = function(x){(x - mean(x))/sd(x)}),.SDcols = c('N34','IOD_west','IOD_east'),by = .(system,forecast_month,target_month)]

fc_dt[,lead_time:= target_month - forecast_month]
fc_dt[lead_time <=0,lead_time := lead_time + 12]
setnames(fc_dt,'target_month','month')

fc_dt[,year:= forecast_year*(month > lead_time) + (forecast_year + 1)*(month <= lead_time)]


fc_dt[,IOD := IOD_west - IOD_east]


### process observations ###
obs_dt = merge(N34_dt,IOD_dt,by = c('year','month'))


# standardize everything:
obs_dt[,c('N34','IOD_west','IOD_east') := lapply(.SD,FUN = function(x){(x - mean(x))/sd(x)}),.SDcols = c('N34','IOD_west','IOD_east'),by = .(month)]


obs_dt[,IOD := IOD_west - IOD_east]


### get scores ###
for(variable in c('IOD','N34'))
{
  
  mse_dt = MSE_ens_fc(fc_dt,fc_cn = variable,
             obs_dt = obs_dt,obs_cn = variable,
             by_cns = c('month','lead_time','system'))
  
  msess_dt = MSESS_ens_fc(fc_dt,fc_cn = variable,
                      obs_dt = obs_dt,obs_cn = variable,
                      by_cns = c('month','lead_time','system'))
  
  crps_dt = CRPS_ens_fc(fc_dt,fc_cn = variable,
             obs_dt = obs_dt,obs_cn = variable,
             by_cns = c('month','lead_time','system'))
  
  crpss_dt = CRPSS_ens_fc(fc_dt,fc_cn = variable,
                        obs_dt = obs_dt,obs_cn = variable,
                        by_cns = c('month','lead_time','system'))
  
  
  
  crps_dt[,lead_time := as.factor(lead_time)]
  crpss_dt[,lead_time := as.factor(lead_time)]
  mse_dt[,lead_time := as.factor(lead_time)]
  msess_dt[,lead_time := as.factor(lead_time)]
  
  ### plotting ###
  plot_dir = '/nr/project/stat/CONFER/plots/paperTCs/'
  
  theme_set(theme_classic(base_size = 14))
  
  for(lt in 1:5)
  {
    
    ## MSESS ##
    ylims = c(min(c(msess_dt[,MSESS],0)),1)
    
    pp = ggplot(msess_dt[lead_time == lt]) + 
      geom_line(aes(x = month,y = MSESS,linetype = system, color = system)) + 
      geom_point(aes(x = month,y = MSESS,shape = system, color = system),size = 2) + 
      scale_x_continuous(breaks = 1:12,labels = 1:12) + 
      theme(panel.grid.major.x = element_line(color = 'gray',size = 0.2)) + 
      scale_y_continuous(limits = ylims) + 
      ggtitle(paste0(variable,' prediction, ',lt,ifelse(lt == 1,yes = ' month',no = ' months'),' ahead'))
    
    if(min(ylims) < 0 ) pp = pp + geom_hline(yintercept = 0,linetype = 'dashed',size = 0.4,color = 'gray')
    
    pdf(paste0(plot_dir,variable,'_MSESS_lt',lt,'.pdf'))  
      print(pp)
    dev.off()
    
    ## CRPSS ##
    
    ylims = c(min(c(crpss_dt[,CRPSS],0)),1)
    
    pp = ggplot(crpss_dt[lead_time == lt]) + 
      geom_line(aes(x = month,y = CRPSS,linetype = system, color = system)) + 
      geom_point(aes(x = month,y = CRPSS,shape = system, color = system),size = 2) + 
      scale_x_continuous(breaks = 1:12,labels = 1:12) + 
      theme(panel.grid.major.x = element_line(color = 'gray',size = 0.2)) + 
      scale_y_continuous(limits = ylims) + 
      ggtitle(paste0(variable,' prediction, ',lt,ifelse(lt == 1,yes = ' month',no = ' months'),' ahead'))
    
    if(min(ylims) < 0 ) pp = pp + geom_hline(yintercept = 0,linetype = 'dashed',size = 0.4,color = 'gray')
    
    pdf(paste0(plot_dir,variable,'_CRPSS_lt',lt,'.pdf'))  
    print(pp)
    dev.off()
    
      
  }
}  

#######################################################