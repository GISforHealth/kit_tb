## #############################################################################
## 
## RUN VARIANCE INFLATION FACTOR ANALYSIS ON COVARIATES
## Author: Kirsten Wiens, adapted by Nat Henry
## Created: Sept 24, 2019
## 
##
## Script designed to facilitate covariate selection using variance inflation 
##   factor (VIF) analysis
## The purpose of the VIF analysis is to test for multicollinearity in covariates
##
## (1) load covariate data cropped to MBG input data by region
## (2) perform VIF analysis and drop covariates above a threshold by region
## (3) take results and summarize in a table that is human-readable
## (4) take results and summarize in MBG input format
##
## References: 
## Faraway's Linear Models in R, Chapter 4: Diagnostics
## https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/
##
## #############################################################################

## SET INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Path to LBD Core repository
core_repo <- '/share/code/geospatial/nathenry/lbd_core/'

# Path to extracted covariates data file (in CSV format). This file should only
#  contain the covariates to be vetted, extracted at locations where data exists.
covs_file <- '/share/geospatial/mbg/tb/kit/model_assets/extracted_covs_full_20190826.csv'

# Path to covariates config file in repository
cov_config_fp <- '/share/code/geospatial/nathenry/kit_tb/covs_kit.csv'

# Path to directory where VIF outputs will be saved.
mbg_input_dir <- save_dir <- '/share/geospatial/mbg/tb/kit/model_assets/vif/'
dir.create(save_dir, showWarnings=FALSE)

# Vector of any covariate names that should be kept regardless of VIF results
covs_to_keep <- 'wp_poverty'

# List of thresholds to use for VIF analyses
thresholds <- seq(1.5,4,by=.5)

# Regions (just one for TB)
regions <- 'pak'
indicator_group <- 'tb'
indicator <- 'kit'
run_date <- "2019_08_26_new_custom_covs"

## SETUP AND LOAD DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Load all MBG functions and libraries
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)

# Load covariate matrix
all_cov_data <- list(fread(covs_file))
names(all_cov_data) <- regions


## FUNCTION FOR RUNNING STEPPWISE VIF SELECTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---------------------------------------------------------------
# Stepwise VIF selection function

# Start function
stepwise_vif_selection <- function(
  thresh, 
  covariate_data,
  reg,
  not_to_drop = NULL,
  trace = TRUE
){
  message(sprintf("\n\n******\nRUNNING VIF WITH THRESHOLD %s\n******\n",thresh))
  # ----------------------------------------------------------------------------
  # Data checks
  
  # make sure input data is a data table
  if(not('data.table' %in% class(covariate_data))){
    covariate_data<-as.data.table(covariate_data)
  }
  
  # remove any covariate values that do not vary accross the region
  var_names <- names(covariate_data)
  for (val in var_names) {
    if (!is.na(val) & length(unique(covariate_data[[val]])) < 2) {
      message(paste0(
        '\nRemoving covariate "', val, '" because it is constant across the region.'
      ))
      covariate_data[, (val) := NULL]
    }
  }
  # ----------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------
  # VIF selection
  
  # load vif analysis package
  library(fmsb, lib.loc='/home/j/temp/nathenry/u5m/visualization/inputs/packages/')
  
  # get initial vif value for all comparisons of variables
  vif_init<-NULL
  
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = covariate_data))))
  }
  
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    
    if(trace==T){ # print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      message(paste0(
        '\nAll variables have VIF < ', thresh,', max VIF ',round(vif_max,2),'\n\n'
      ))
    }
    cat(var_names)
    
  } else {
    
    # backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(covariate_data)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = covariate_data))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      
      # if there are any covs that we don't want dropped, set their vif to 0
      if (!is.null(not_to_drop)) {
        for (k in not_to_drop) {
          if (trace==T) message(paste0('Making sure we keep ', k))
          vif_vals[which(vif_vals[,1] == k),2] <- 0
        }
      }
      
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ # print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        message(paste0('removed: ',vif_vals[max_row,1], ' ', vif_max,'\n\n'))
        flush.console()
      }
      
      keep_cols <- var_names[-which(var_names == vif_vals[max_row,1])]
      covariate_data<-covariate_data[, keep_cols, with = FALSE]
      
    }
    
    message(paste0('\nCovariates selected for ',indicator,' in ',reg,
                   ' with a threshold of ',thresh,':'))
    cat(names(covariate_data))
    
  }
  
  # vector of selected covariates
  selected_covs <- c(names(covariate_data))
  # ----------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------
  # Add selection results to table
  
  # add 0/1 for whether or not covariate was selected to a data table
  d1 <- data.table(cov_names)
  d2 <- data.table(selected_covs)
  d2[, included := TRUE]
  dt <- merge(d1, d2, by.x = 'cov_names', by.y = 'selected_covs', all = T)
  dt[is.na(included), included := FALSE]
  dt <- dcast(melt(dt, id.vars = 'cov_names'), variable ~ cov_names)
  dt[, variable := NULL]
  
  # add threshold, region, number of covariates selected, and run date that corresponds to data
  dt[, region := reg]
  dt[, vif_threshold := thresh]
  dt[, num_covs_selected := length(selected_covs)]
  dt[, data_run_date := run_date]
  
  # bind to original data table
  selection_results <- rbindlist(list(selection_results, dt), use.names = TRUE)
  # ----------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------
  # Save mbg input files
  
  # create mbg input files
  selected_covs <- data.table(covariate = selected_covs)
  cov_config <- merge(covs, selected_covs, by = 'covariate')
  
  # save mbg input files
  dir.create(paste0(mbg_input_dir, 'vif_', i, '/'), showWarnings = FALSE)
  write.csv(cov_config, paste0(mbg_input_dir, 'vif_', i, '/covs_', indicator_group, '_', region, '.csv'))
  # ----------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------
  # End function
  
  # return a vector of covariates selected
  return(selection_results)
}
# ---------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Stepwise covariate removal using VIF

sr_all <- vector('list',length=length(thresholds))

# loop over thresholds
for (idx in 1:length(thresholds)) {
  # report threshold
  i <- thresholds[idx]
  message(paste0('\nThreshold: ', i))

  # load input covariate names
  covs <- fread(cov_config_fp)
  covs <- covs[(include == TRUE) | (include=='T'),]
  cov_names <- covs[, covariate]
  
  # set up table
  col_names <- c(cov_names, 'region', 'vif_threshold', 'num_covs_selected', 'data_run_date')
  selection_results <- setNames(
    data.table(matrix(nrow = 0, ncol = length(col_names))),
    col_names
  )
  
  # run vif selection analysis
  for (region in regions) {
    message(paste0('\n\n', region))
    
    # run stepwise vif selection function
    selection_results <- stepwise_vif_selection(thresh = i, 
                                                covariate_data = all_cov_data[[region]],
                                                reg = region,
                                                trace = T,
                                                not_to_drop = covs_to_keep)
    sr_all[[idx]] <- selection_results
  }
  
  # save selection results
  write.csv(selection_results, paste0(save_dir, 'selection_results_', indicator_group, '_vif_', i, '.csv'))
  
  # reset for next threshold
  message('\nSelection complete. Results saved.')
  
}
# ------------------------------------------------------------------------------

sr_all <- rbindlist(sr_all)
write.csv(sr_all, paste0(save_dir, 'selection_results_', indicator_group, '_vif_ALL_THRESHOLDS.csv'))
