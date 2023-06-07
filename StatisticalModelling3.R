
# tryout of several statistical models for ON precip, based on IOD and N34. We fit gridpointwise linear regression models on N34, IOD, and both, 
# as well as EOF-regression, i.e. regression on factor loadings. The results are derived based on SSTs from August,...,November.
# The script plots regression coefficients as well as diverse skill maps, and MSE bootstrap results.

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


new_plot_dir = paste0(plot_dir,'IOD_N34_regression_coefs/')

dir.create(new_plot_dir,showWarnings = F)

for(mon in c('Aug','Sep','Oct','Nov'))
{
  print(mon)
  temp = glm_dt(obs_dt,as.formula(paste0('prec~IOD_',mon)),bycols = c('lon','lat'))
  pp = ggplot_dt(temp,paste0('IOD_',mon,'_coef'),rr= c(-1.5,1.5),centering = 0)
  
  png(file = paste0(new_plot_dir,'IOD_',mon,'_coef.png'))
    print(pp)
  dev.off()
  
  temp = glm_dt(obs_dt,as.formula(paste0('prec~N34_',mon)),bycols = c('lon','lat'))
  pp = ggplot_dt(temp,paste0('N34_',mon,'_coef'),rr= c(-1.5,1.5),centering = 0)
  png(file = paste0(new_plot_dir,'N34_',mon,'_coef.png'))
    print(pp)
  dev.off()
}


# joint regression

for(mon in c('Aug','Sep','Oct','Nov'))
{
  print(mon)
  temp = glm_dt(obs_dt,as.formula(paste0('prec~IOD_',mon,' + N34_',mon)),bycols = c('lon','lat'))
  pp = ggplot_dt(temp,paste0('IOD_',mon,'_coef'),rr= c(-1.5,1.5),centering = 0)
  
  png(file = paste0(new_plot_dir,'jr_IOD_',mon,'_coef.png'))
  print(pp)
  dev.off()
  
  pp = ggplot_dt(temp,paste0('N34_',mon,'_coef'),rr= c(-1.5,1.5),centering = 0)
  png(file = paste0(new_plot_dir,'jr_N34_',mon,'_coef.png'))
  print(pp)
  dev.off()
}

###############################################
### assess skill of gridpointwise lr-models ###
###############################################

val_years = 1993:2016

skill_dt = data.table()
bt_dt = data.table()

for(mon in c('Aug','Sep','Oct','Nov'))
{
  print(mon)
  # IOD
  
  dt_lm_IOD = loyo_dt_lm(obs_dt,as.formula(paste0('prec~IOD_',mon)),bycols = c('lon','lat'),mc_cores = 8)
  dt_lm_IOD = merge(dt_lm_IOD,obs_dt[,.(lon,lat,year,prec)],by = c('year','lon','lat'))[year %in% val_years]
  
  dt_lm_IOD[,MSE := (prediction - prec)^2]
  
  bt_temp = bootstrap_scores_dt(dt_lm_IOD,score_col = 'MSE')
  bt_temp[,model:= paste0(mon,'_IOD')]
  
  bt_dt = rbindlist(list(bt_dt,bt_temp))
  skill_dt = rbindlist(list(skill_dt,dt_lm_IOD[,.(year,lon,lat,prediction,prec,MSE)][,model := paste0('lm_',mon,'_IOD')]))
  
  # N34
  
  dt_lm_N34 = loyo_dt_lm(obs_dt,as.formula(paste0('prec~N34_',mon)),bycols = c('lon','lat'),mc_cores = 8)
  dt_lm_N34 = merge(dt_lm_N34,obs_dt[,.(lon,lat,year,prec)],by = c('year','lon','lat'))[year %in% val_years]
  
  dt_lm_N34[,MSE := (prediction - prec)^2]
  
  bt_temp = bootstrap_scores_dt(dt_lm_N34,score_col = 'MSE')
  bt_temp[,model:= paste0(mon,'_N34')]
  
  bt_dt = rbindlist(list(bt_dt,bt_temp))
  skill_dt = rbindlist(list(skill_dt,dt_lm_N34[,.(year,lon,lat,prediction,prec,MSE)][,model := paste0('lm_',mon,'_N34')]))
  
  # both IOD and N34
  
  dt_lm_both = loyo_dt_lm(obs_dt,as.formula(paste0('prec~IOD_',mon,' + N34_',mon)),bycols = c('lon','lat'),mc_cores = 8)
  dt_lm_both = merge(dt_lm_both,obs_dt[,.(lon,lat,year,prec)],by = c('year','lon','lat'))[year %in% val_years]
  
  dt_lm_both[,MSE := (prediction - prec)^2]
  
  bt_temp = bootstrap_scores_dt(dt_lm_both,score_col = 'MSE')
  bt_temp[,model:= paste0(mon,'_both')]
  
  bt_dt = rbindlist(list(bt_dt,bt_temp))
  skill_dt = rbindlist(list(skill_dt,dt_lm_both[,.(year,lon,lat,prediction,prec,MSE)][,model := paste0('lm_',mon,'_both')]))
}


# append climatology

obs_dt = loyo(obs_dt,FUN = mean,'prec',bycols = c('lon','lat'))

setnames(obs_dt,'prec_new','clim')

temp = copy(obs_dt[,.(year,lon,lat,prec,clim)][,MSE := (clim - prec)^2][year %in% val_years])
bt_temp = bootstrap_scores_dt(temp,score_col = 'MSE')
bt_temp[,model:= 'clim']
bt_dt = rbindlist(list(bt_dt,bt_temp))

pp = ggplot(bt_dt) + geom_boxplot(aes(x = model,y = bootstrap_samples,color = model))
pp
 

###########################################
### assess skill of eof-based lr-models ###
###########################################

val_years = 1993:2016

bt_dt_eof = data.table()

for(mon in c('Aug','Sep','Oct','Nov'))
{
  print(mon)
  for(nv0 in c(38)) # 38 is the max number: this is the rank of the data matrix. Turns out that reducing dimensions leads to worse results here, I checked.
  {
    print(nv0)
  # IOD
  
  dt_eof_IOD = loyo_eof_regression(obs_dt,as.formula(paste0('prec~IOD_',mon)),nv = nv0,mc_cores = 8)
  dt_eof_IOD = merge(dt_eof_IOD,obs_dt[,.(lon,lat,year,prec)],by = c('year','lon','lat'))[year %in% val_years]
  
  dt_eof_IOD[,MSE := (prediction - prec)^2]
  
  bt_temp = bootstrap_scores_dt(dt_eof_IOD,score_col = 'MSE')
  bt_temp[,model:= paste0(mon,'_IOD')][,nv := nv0]
  
  bt_dt_eof = rbindlist(list(bt_dt_eof,bt_temp))
  skill_dt = rbindlist(list(skill_dt,dt_eof_IOD[,.(year,lon,lat,prediction,prec,MSE)][,model := paste0('lm_',mon,'_IOD_nv',nv0)]))
  
  # N34
  
  dt_eof_N34 = loyo_eof_regression(obs_dt,as.formula(paste0('prec~N34_',mon)),nv = nv0,mc_cores = 8)
  dt_eof_N34 = merge(dt_eof_N34,obs_dt[,.(lon,lat,year,prec)],by = c('year','lon','lat'))[year %in% val_years]
  
  dt_eof_N34[,MSE := (prediction - prec)^2]
  
  bt_temp = bootstrap_scores_dt(dt_eof_N34,score_col = 'MSE')
  bt_temp[,model:= paste0(mon,'_N34')][,nv := nv0]
  
  bt_dt_eof = rbindlist(list(bt_dt_eof,bt_temp))
  skill_dt = rbindlist(list(skill_dt,dt_eof_N34[,.(year,lon,lat,prediction,prec,MSE)][,model := paste0('lm_',mon,'_N34_nv',nv0)]))
  
  # both IOD and N34
  
  dt_eof_both = loyo_eof_regression(obs_dt,as.formula(paste0('prec~IOD_',mon,' + N34_',mon)),nv = nv0,mc_cores = 8)
  dt_eof_both = merge(dt_eof_both,obs_dt[,.(lon,lat,year,prec)],by = c('year','lon','lat'))[year %in% val_years]
  
  dt_eof_both[,MSE := (prediction - prec)^2]
  
  bt_temp = bootstrap_scores_dt(dt_eof_both,score_col = 'MSE')
  bt_temp[,model:= paste0(mon,'_both')][,nv := nv0]
  
  bt_dt_eof = rbindlist(list(bt_dt_eof,bt_temp))
  skill_dt = rbindlist(list(skill_dt,dt_eof_both[,.(year,lon,lat,prediction,prec,MSE)][,model := paste0('lm_',mon,'_both_nv',nv0)]))
  }
}




# append climatology

temp = copy(obs_dt[,.(year,lon,lat,prec,clim)][,MSE := (clim - prec)^2][year %in% val_years])
bt_temp = bootstrap_scores_dt(temp,score_col = 'MSE')
bt_temp[,model:= 'clim']
bt_temp[,nv:= NA]

bt_dt_eof = rbindlist(list(bt_dt_eof,bt_temp))


# compare both 

bt_dt_joint = rbindlist(list(bt_dt[,model := paste0(model,'_lr')][, nv := NA],bt_dt_eof[model != 'clim',model := paste0(model,'_eof')]))

month_as_number = function(model_names)
{
  return( 8*grepl('Aug',model_names,fixed = T) +
          9*grepl('Sep',model_names,fixed = T) +
          10 * grepl('Oct',model_names,fixed = T) +
          11 * grepl('Nov',model_names,fixed = T))
}

bt_dt_joint[,sst_month := month_as_number(model)]

#################
##### plots #####
#################

### start out with boxplots ###

### get August-initialized plot ###

old_model_names = unique(bt_dt_joint[sst_month %in% c(0,8),model])
# "Aug_IOD_lr"   "Aug_N34_lr"   "Aug_both_lr"  "clim_lr"      "Aug_IOD_eof"  "Aug_N34_eof"  "Aug_both_eof"
new_model_names = c('lr_IOD','lr_N34','lr_both','clim','eofr_IOD','eofr_N34','eofr_both')
model_order = c(4,1,2,3,5,6,7)

rr = c(0.75,1.55)

plot_dt = bt_dt_joint[sst_month %in% c(0,8),][,mod_ind := model_order[match(model,old_model_names)]]
plot_dt[,model:=new_model_names[match(model,old_model_names)]]

setkey(plot_dt,mod_ind)

pp = ggplot(plot_dt) + geom_boxplot(aes(x = model,y = bootstrap_samples,color = model))  + ylab('bootstrapped MSE') + ggtitle('August SST')
pp
pdf(paste0(plot_dir,'Skill_model_comparison/MSE_boxplots_August_Statistical.pdf'),width = 14)
pp
dev.off()

# November-initialized

old_model_names = unique(bt_dt_joint[sst_month %in% c(0,11),model])
# "Nov_IOD_lr"   "Nov_N34_lr"   "Nov_both_lr"  "clim_lr"      "Nov_IOD_eof"  "Nov_N34_eof"  "Nov_both_eof"
new_model_names = c('lr_IOD','lr_N34','lr_both','clim','eofr_IOD','eofr_N34','eofr_both')
model_order = c(4,1,2,3,5,6,7)

rr = c(0.75,1.55)

plot_dt = bt_dt_joint[sst_month %in% c(0,11),][,mod_ind := model_order[match(model,old_model_names)]]
plot_dt[,model:=new_model_names[match(model,old_model_names)]]

setkey(plot_dt,mod_ind)

pp = ggplot(plot_dt) + geom_boxplot(aes(x = model,y = bootstrap_samples,color = model)) + ylim(rr) + ylab('bootstrapped MSE') + ggtitle('November SST')


pdf(paste0(plot_dir,'Skill_model_comparison/MSE_boxplots_November_Statistical.pdf'),width = 14)
pp
dev.off()

#########################################

# skill-score maps:

skill_dt = mask_precip(skill_dt)[(!mask)][,mask:=NULL]

msess = MSESS_dt(skill_dt,fc_col = 'prediction',obs_col = 'prec',by_cols = c('lon','lat','model'))

pp = ggplot_dt(msess[model == 'lm_Aug_IOD_nv38'],'MSESS',centering = 0,rr = c(-1,1)) + ggtitle('EOFR on August IOD')

png(paste0(plot_dir,'Skill_model_comparison/MSESS_maps/eof_Aug_IOD_nv30.png'))
  print(pp)
dev.off()

pp = ggplot_dt(msess[model == 'lm_Aug_IOD'],'MSESS',centering = 0,rr = c(-1,1)) + ggtitle('LR on August IOD')
png(paste0(plot_dir,'Skill_model_comparison/MSESS_maps/lm_Aug_IOD.png'))
  print(pp)
dev.off()



pp = ggplot_dt(msess[model == 'lm_Aug_N34'],'MSESS',centering = 0,rr = c(-1,1)) + ggtitle('LR on August N3.4')
png(paste0(plot_dir,'Skill_model_comparison/MSESS_maps/lm_Aug_N34.png'))
print(pp)
dev.off()



pp = ggplot_dt(msess[model == 'lm_Aug_both'],'MSESS',centering = 0,rr = c(-1,1)) + ggtitle('LR on August IOD and N3.4')
png(paste0(plot_dir,'Skill_model_comparison/MSESS_maps/lm_Aug_both.png'))
print(pp)
dev.off()


pp = ggplot_dt(msess[model == 'lm_Nov_both'],'MSESS',centering = 0,rr = c(-1,1)) + ggtitle('LR on November IOD and N3.4')
png(paste0(plot_dir,'Skill_model_comparison/MSESS_maps/lm_Nov_both.png'))
print(pp)
dev.off()

#### save data ####

save_dir = paste0(data_dir,'TCpaper/')

fwrite(skill_dt,file = paste0(save_dir,'spatial_skill_statistical_models.csv'))
fwrite(bt_dt_joint,file = paste0(save_dir,'bootstrapped_MSEs_statistical_models.csv'))
