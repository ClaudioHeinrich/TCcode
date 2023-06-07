
# function for running a matrix linear regression

matrix_linear_regression = function(dt, formula, 
                                    bycols = intersect(c('lon','lat'),colnames(dt)),
                                    loyo = T, 
                                    alongcols = 'year',
                                    intercept = T)
{
  
  dependent_name = all.vars(formula)[1]
  covariate_names = all.vars(formula)[-1]
  
  # wrapper, useful for leva-one-year-out:
  wrapper_fct = function(dt_temp)
  {
    setkeyv(dt_temp,c(alongcols,bycols))
    
    # bring dependent variables into the right matrix shape:
    cast_formula = as.formula(paste0(paste(alongcols,collapse = ' + '), ' ~ ', paste(bycols,collapse = ' + ')))
    y_dt = dcast(dt_temp,cast_formula,value.var = dependent_name)  
    y_mat = as.matrix(y_dt[,(length(alongcols) + 1) : ncol(y_dt)])
    
    # bring covariates into the right matrix shape:
    x_dt = dt_temp[]dcast(dt_temp,cast_formula,value.var = covariate_names)  
    x_mat = as.matrix(x_dt[,(length(alongcols) + 1) : ncol(x_dt)])
    
    
    if(intercept) mod = lm(y_mat ~ x_mat)
    if(!intercept) mod = lm(y_mat ~ 0 + x_mat)
    
    reg_coefs = as.data.table(mod$coefficients)  
  }
  
  
  
  
  
  if(intercept) reg_coef_dt[,coef := c('intercept','beta1')]
  if(!intercept) reg_coef_dt[,coef := c('beta1')]
  
  reg_coef_dt = melt(reg_coef_dt,id.vars = 'coef')
  
  reg_coef_dt[,variable:= as.character(variable)]
  
  dt_temp = data.table(matrix(unlist(strsplit(x = reg_coef_dt[,variable],split = '_')),ncol = 3,byrow = T))
  setnames(dt_temp,c('season','lon','lat'))
  dt_temp[,c('lon','lat') := lapply(.SD,as.numeric),.SDcols = c('lon','lat')]
  
  reg_coef_dt = data.table(reg_coef_dt,dt_temp)
  reg_coef_dt[,variable:=NULL]
  
  return(reg_coef_dt)
}