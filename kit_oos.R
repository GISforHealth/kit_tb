################################################################################
## Generate out-of-sample predictive metrics for KIT TB project               ##
## Author: Nathaniel Henry                                                    ##
## Created: August 2019                                                       ##
##                                                                            ##
## NOTE: This code pulls from many functions in the LBDCore package, which is ##
##   still in development. Key functions from this package will be added to   ##
##   this repository as a demonstration. This package is on track to be       ##
##   publicly released as an open-source tool.                                ##
##                                                                            ##
################################################################################

rd <- run_date 

## #############################################################################
## MBG SETUP
## #############################################################################
# save some other important drive locations
root           <- '<<FILEPATH>>'
sharedir       <- '<<FILEPATH>>'

## Load libraries and MBG project functions.
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)
# In the future, we can also load from files named "XX_functions.R" in our KIT 
#  repository
# load_mbg_functions(ig_repo)

## Read config file and save all parameters in memory
config <- set_up_config(
  repo=ig_repo,
  core_repo=core_repo,
  indicator_group=indicator_group,
  indicator=indicator,
  config_file=paste0(ig_repo,config_name,'.csv'),
  covs_file=paste0(ig_repo,covs_config_name,'.csv'),
  post_est_only=FALSE,
  run_date=run_date,
  push_to_global_env=TRUE,
  run_tests=FALSE
)
slots      <- as.numeric(slots)
inla_cores <- as.numeric(inla_cores)

## Create directory structure for this model run, including an output for each agebin
create_dirs(indicator_group = indicator_group,
            indicator       = indicator)
dir.create(sprintf('%s/output/%s', sharedir, run_date), showWarnings=FALSE)
oos_dir <- sprintf('%s/output/%s/oos/',sharedir,run_date)
dir.create(oos_dir, showWarnings=FALSE)

reg         <- regions
age         <- 0
holdout     <- 0
samples <- as.numeric(samples)

## grab arguments
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
outputdir <- '<<FILEPATH>>'
dir.create(outputdir, showWarnings = FALSE)

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

# print out session info so we have it on record
sessionInfo()

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
check_config()

# We need to be in the singularity image, and specifically the LBD one if using TMB
if(!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if(as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Print the core_repo hash and check it
message("Printing git hash for 'core_repo' and checking against LBD Core Code master repo")
record_git_status(core_repo = core_repo, check_core_repo = TRUE)

## Make sure this inla patch is implemented if running on geos
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

## cores to use
cores_to_use <- Sys.getenv("SGE_HGR_fthread")
message(paste("Model set to use", cores_to_use, "cores"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load input data once to get a full list of holdouts
load('<<FILEPATH>>')
all_holdouts <- sort(na.omit(unique(df$cluster_number)))

for(holdout_i in all_holdouts){
  message(sprintf("\n######\n*** PREDICTING FOR HOLDOUT %i ***\n######\n",holdout_i))

  ## reload data an prepare for MBG
  load('<<FILEPATH>>')

  df <- df[cluster_number != holdout_i, ]

  ## convert stackers to transform space, if desired
  ## NOTE: we do this here to ensure that the stacker rasters are saved in prevalence/untransformed space
  ## this is useful for diagnostics and other code that was built expecting the untransformed rasters
  if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)){

    ## transform the rasters
    for (ii in child_model_names) {
      
      ## Preserve variable names in the raster first
      tmp_rastvar <- names(cov_list[[ii]])
      
      ## Logit
      cov_list[[ii]] <- logit(cov_list[[ii]])
      
      ## Reassign names
      names(cov_list[[ii]]) <- tmp_rastvar
      rm(tmp_rastvar)
    }

    ## transform the stacker values that are in df
    stacker_cols <- grep(paste0("(", paste(child_model_names, collapse="|"), ")(.*_pred)"), names(df), value=T)
    df[, (stacker_cols) := lapply(.SD, logit), .SDcols = stacker_cols]

  }

  input_data <- build_mbg_data_stack_tmb(
    d = df,
    fes = all_fixed_effects,
    indic = indicator, 
    mesh = mesh_s,
    cov_constraints = covariate_constraint_vectorize(config)
  )

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG BY HOLDOUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## combine all the inputs, other than cs_df these are not used if you are using TMB
  input_data$clamper$int <- c(-200, 200)
  stacked_input  <- input_data$Data
  cs_df          <- input_data$cs_df

  ## Generate other inputs necessary
  outcome <- df[[indicator]] # N+_i - event obs in cluster
  N       <- df$N                  # N_i - total obs in cluster
  weights <- df$weight

  ## catch in case there is no weight column
  if(is.null(weights)){
    weights = rep(1,nrow(df))
  }

  ## Set the number of cores to be equal to input;
  ## If missing, then revert to cores_to_use value
  if(Sys.getenv("OMP_NUM_THREADS") != "") {
      setompthreads(Sys.getenv("OMP_NUM_THREADS"))
  } else {
      setompthreads(cores_to_use)
  }

  if(Sys.getenv("MKL_NUM_THREADS") != "") {
      setmklthreads(Sys.getenv("MKL_NUM_THREADS"))
  } else {
      setmklthreads(cores_to_use)
  }

  ## Fit MBG model
  # run the model
  system.time(
    model_fit <- suppressMessages(suppressWarnings(fit_mbg_tmb(
      lbdcorerepo = core_repo,
      cpp_template = 'mbg_tmb_model',
      tmb_input_stack = input_data,
      control_list = list(trace=0, maxit=1000, itnmax=1000, eval.max=1000, iter.max=1000, abstol=-1E10, reltol=1E-8),
      optimizer = 'optimx',
      ADmap_list = NULL,
      newton_steps = 0,
      sparse_ordering = as.logical(sparse_ordering)
    )))
  )

  # clamping
  clamp_covs <- TRUE

  max_chunk <- 50
  samples   <- as.numeric(samples)

  if(is.character(year_list_predict)) year_list_predict <- eval(parse(text=year_list_predict))

  chunks <- rep(max_chunk, samples %/% max_chunk)
  if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
  message(sprintf("***** Predicting for holdout %i *****",holdout_i))
  pm <- lapply(chunks, function(samp) {
    suppressMessages(suppressWarnings(predict_mbg_tmb(
      samples              = samp,
      seed                 = NULL,
      tmb_input_stack      = input_data,
      model_fit_object     = model_fit,
      fes                  = all_fixed_effects,
      sr                   = simple_raster,
      yl                   = year_list_predict,
      zl                   = z_list,
      covs_list            = cov_list,
      clamp_covs           = clamp_covs,
      cov_constraints = covariate_constraint_vectorize(config)
    )))
  })

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Make cell preds and a mean raster
  cell_pred <- do.call(cbind, pm)
  unraked_mean_rr <- insertRaster(
    simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period))
  )
  unraked_median_rr <- insertRaster(
    simple_raster,
    matrix(
      apply(cell_pred, 1, function(x) quantile(x, .5, na.rm=T)),
      ncol = max(period_map$period)
    )
  )

  writeRaster(
    unraked_mean_rr,
    file=sprintf(
      '%s/%s_mean_unraked_%s_%s_holdout%s',oos_dir,indicator,min(year_list_predict),
      max(year_list_predict),holdout_i
    ),
    format = "GTiff", overwrite=TRUE
  )
  writeRaster(
    unraked_median_rr,
    file=sprintf(
      '%s/%s_median_unraked_%s_%s_holdout%s',oos_dir,indicator,min(year_list_predict),
      max(year_list_predict),holdout_i
    ),
    format = "GTiff", overwrite=TRUE
  )

  message(sprintf("*** DONE PREDICTING AND SAVING FOR HOLDOUT %i", holdout_i))
}

message("\n*** DONE WITH HOLDOUT FITTING ***\n")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~ GENERATE OOS VALIDITY AND SAVE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load GBD 2017 admin0 data for 2011 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source('<<FILEPATH>>')
## Get GBD location ID for Pakistan
ad0_code <- get_adm0_codes(reg)
gbd_loc_ids <- load_adm0_lookup_table()[gadm_geoid %in% ad0_code, loc_id]

## Load GBD 2017 admin0 data
over_15_ags <- get_age_metadata(age_group_set_id=12)[
  age_group_years_start >= 15, age_group_id
]
## Set cause_ids associated with TB all-forms (excluding HIV-TB)
tb_cause_ids <- c(934,946,947,948,949,950)

# Pull TB by age and sex
gbd2017_ad0_tb_agesex_specific <- get_outputs(
  'cause',
  cause_id = tb_cause_ids,
  age_group_id = over_15_ags,
  location_id = gbd_loc_ids,
  sex_id = 1:2,
  year_id = 2011,
  measure_id = 5, # Prevalence
  metric_id = 3, # Rate
  gbd_round_id = 5
)[, .(val=sum(val)), by=.(age_group_id, sex_id, year_id, location_id)]

# Pull population by age and sex
gbd2017_pops <- get_population(
  age_group_id = over_15_ags,
  location_id = gbd_loc_ids,
  sex_id = 1:2,
  year_id = 2011, 
  gbd_round_id = 5
)[, .(age_group_id, sex_id, year_id, location_id, population)]
# Run population-weighted aggregation
gbd2017_ad0 <- merge(
  x=gbd2017_ad0_tb_agesex_specific, y=gbd2017_pops,
  by=c('age_group_id','sex_id','location_id','year_id')
)[,
  .(gbd_est = weighted.mean(val, w=population)), by=.(location_id,year_id)
]$gbd_est

## Pull population raster in 2011 for weighting
pop_11 <- readRDS(
  sprintf('%s/model_assets/pop_15pl_%s.RDS',sharedir,reg)
)[[which(year_list_predict==2011)]]

## Load the dataset one last time
load('<<FILEPATH>>')

data_by_holdout <- lapply(all_holdouts, function(hh) df[cluster_number==hh,])
names(data_by_holdout) <- paste0('holdout_',all_holdouts)
data_by_holdout_filled <- vector('list', length=length(all_holdouts))
names(data_by_holdout_filled) <- paste0('holdout_',all_holdouts)

## Combine out-of-sample predictions with true holdout data
for(holdout_i in all_holdouts){
  message(sprintf('Prep holdout %s...',holdout_i))
  # Pull data for holdout
  data_sub <- data_by_holdout[[paste0('holdout_',holdout_i)]]
  # Pull rasters in 2011
  mean_r_11 <- raster::brick(sprintf(
    '%s/%s_mean_unraked_%s_%s_holdout%s.tif',oos_dir,indicator,min(year_list_predict),
    max(year_list_predict),holdout_i
  ))[[which(year_list_predict==2011)]]
  median_r_11 <- raster::brick(sprintf(
    '%s/%s_median_unraked_%s_%s_holdout%s.tif',oos_dir,indicator,min(year_list_predict),
    max(year_list_predict),holdout_i
  ))[[which(year_list_predict==2011)]]
  ## Create raking factors
  stopifnot(all(dim(mean_r_11) == dim(pop_11)))
  stopifnot(all(dim(median_r_11) == dim(pop_11)))
  rf_dt <- na.omit(data.table(
    mean_est = as.vector(mean_r_11),
    median_est = as.vector(median_r_11),
    pop = as.vector(pop_11)
  ))[, 
    .(lbd_est_mean = weighted.mean(mean_est, w=pop), 
      lbd_est_median=weighted.mean(median_est, w=pop))
  ]
  rf_dt$gbd_est <- gbd2017_ad0
  mean_rf <- rf_dt[, gbd_est / lbd_est_mean ]
  median_rf <- rf_dt[, gbd_est / lbd_est_median ]
  message(sprintf('  -- Mean RF=%.3f -- Median RF=%.3f',mean_rf,median_rf))
  # Convert to spatialPoints
  sub_pts <- sp::SpatialPoints(cbind(data_sub$longitude, data_sub$latitude))
  # Extract mean and median data at points
  data_sub$oos_est_mean <- raster::extract(x=mean_r_11, y=sub_pts)
  data_sub$oos_est_median <- raster::extract(x=median_r_11, y=sub_pts)
  data_sub[, oos_est_mean_raked := oos_est_mean * mean_rf ]
  data_sub[, oos_est_median_raked := oos_est_median * median_rf ]
  # Add to output list
  data_by_holdout_filled[[paste0('holdout_',holdout_i)]] <- data_sub
}

## Create combined out-of-sample predictions set
# Combine
oos_predict_full <- rbindlist(data_by_holdout_filled)

# Get data from in-sample model
is_mean_unraked_r <- raster::brick(
  sprintf('%s/output/%s/kit_mean_unraked_2010_2018.tif',sharedir,rd)
)[[which(year_list_predict==2011)]]
is_median_unraked_r <- raster::brick(
  sprintf('%s/output/%s/kit_median_unraked_2010_2018.tif',sharedir,rd)
)[[which(year_list_predict==2011)]]
is_mean_raked_r <- raster::brick(
  sprintf('%s/output/%s/kit_mean_raked_2010_2018.tif',sharedir,rd)
)[[which(year_list_predict==2011)]]
is_median_raked_r <- raster::brick(
  sprintf('%s/output/%s/kit_median_raked_2010_2018.tif',sharedir,rd)
)[[which(year_list_predict==2011)]]

oos_predict_full$is_est_mean <- raster::extract(
  is_mean_unraked_r, y=SpatialPoints(cbind(oos_predict_full$longitude, oos_predict_full$latitude))
)
oos_predict_full$is_est_median <- raster::extract(
  is_median_unraked_r, y=SpatialPoints(cbind(oos_predict_full$longitude, oos_predict_full$latitude))
)
oos_predict_full$is_est_mean_raked <- raster::extract(
  is_mean_raked_r, y=SpatialPoints(cbind(oos_predict_full$longitude, oos_predict_full$latitude))
)
oos_predict_full$is_est_median_raked <- raster::extract(
  is_median_raked_r, y=SpatialPoints(cbind(oos_predict_full$longitude, oos_predict_full$latitude))
)
oos_predict_full[, data_est := kit / N ]

# Save to file
write.csv(
  oos_predict_full, row.names=FALSE,
  file=sprintf('%s/out_of_sample_predicted_full.csv', oos_dir)
)


## Generate estimates of predictive validity
## Calculate Mean Squared Error
weighted_mse <- function(data, est, weight){
  full_dt <- na.omit(data.table(
    data = data,
    est = est,
    weight = weight
  ))
  total_weight <- sum(full_dt$weight)
  full_dt[, mse_component := (est - data)**2 * (weight/total_weight)]
  return(sum(full_dt$mse_component))
}

## Calculate Relative Squared Error
weighted_rse <- function(data, est, weight){
  full_dt <- na.omit(data.table(
    data = data,
    est = est,
    weight = weight
  ))
  ybar <- full_dt[, weighted.mean(data, w=weight)]
  total_weight <- sum(full_dt$weight)
  full_dt[, rse_num := (est - data)**2 * (weight/total_weight)]
  full_dt[, rse_denom := (data - ybar)**2 * (weight/total_weight)]
  return(sum(full_dt$rse_num)/sum(full_dt$rse_denom))
}

## Fill output data.table with MSE and RSE
oospf <- copy(oos_predict_full)
oospf <- oospf[,
  .(oos_est_mean = weighted.mean(oos_est_mean, w=weight),
    oos_est_median = weighted.mean(oos_est_median, w=weight),
    oos_est_mean_raked = weighted.mean(oos_est_mean_raked, w=weight),
    oos_est_median_raked = weighted.mean(oos_est_median_raked, w=weight),
    is_est_mean = weighted.mean(is_est_mean, w=weight),
    is_est_median = weighted.mean(is_est_median, w=weight),
    is_est_mean_raked = weighted.mean(is_est_mean_raked, w=weight),
    is_est_median_raked = weighted.mean(is_est_median_raked, w=weight),
    data_est = weighted.mean(data_est, w=weight),
    kit = sum(kit*weight),
    N = sum(N*weight)
  ),
  by=.(cluster_number)
]

oos_metrics <- data.table(
  summ_metric = c("mean","median","mean (raked)","median (raked)"),
  mse = c(
    weighted_mse(data=oospf[,kit/N], est=oospf$oos_est_mean, weight=oospf$N),
    weighted_mse(data=oospf[,kit/N], est=oospf$oos_est_median, weight=oospf$N),
    weighted_mse(data=oospf[,kit/N], est=oospf$oos_est_mean_raked, weight=oospf$N),
    weighted_mse(data=oospf[,kit/N], est=oospf$oos_est_median_raked, weight=oospf$N)
  ),
  rse = c(
    weighted_rse(data=oospf[,kit/N], est=oospf$oos_est_mean, weight=oospf$N),
    weighted_rse(data=oospf[,kit/N], est=oospf$oos_est_median, weight=oospf$N),
    weighted_rse(data=oospf[,kit/N], est=oospf$oos_est_mean_raked, weight=oospf$N),
    weighted_rse(data=oospf[,kit/N], est=oospf$oos_est_median_raked, weight=oospf$N)
  )
)
oos_metrics[, root_mse := sqrt(mse)]
oos_metrics[, root_rse := sqrt(rse)]

is_metrics <- data.table(
  summ_metric = c("mean","median","mean (raked)","median (raked)"),
  mse = c(
    weighted_mse(data=oospf[,kit/N], est=oospf$is_est_mean, weight=oospf$N),
    weighted_mse(data=oospf[,kit/N], est=oospf$is_est_median, weight=oospf$N),
    weighted_mse(data=oospf[,kit/N], est=oospf$is_est_mean_raked, weight=oospf$N),
    weighted_mse(data=oospf[,kit/N], est=oospf$is_est_median_raked, weight=oospf$N)
  ),
  rse = c(
    weighted_rse(data=oospf[,kit/N], est=oospf$is_est_mean, weight=oospf$N),
    weighted_rse(data=oospf[,kit/N], est=oospf$is_est_median, weight=oospf$N),
    weighted_rse(data=oospf[,kit/N], est=oospf$is_est_mean_raked, weight=oospf$N),
    weighted_rse(data=oospf[,kit/N], est=oospf$is_est_median_raked, weight=oospf$N)
  )
)
is_metrics[, root_mse := sqrt(mse)]
is_metrics[, root_rse := sqrt(rse)]

write.csv(oos_metrics, row.names=FALSE, file=sprintf('%s/oos_metrics.csv',oos_dir))


## Make diagnostic plots



ad0_sf <- sf::st_read('<<FILEPATH>>')

## Plot IS and OOS maps
# Function to make a diagnostic map
make_diag_map <- function(df, is_col, oos_col, fp){
  map_df <- copy(df)
  setnames(map_df, is_col,'is_col')
  setnames(map_df, oos_col, 'oos_col')
  map_df[, data_est := data_est * 1E5 ]
  map_df[, is_col := is_col * 1E5 ]
  map_df[, oos_col := oos_col * 1E5 ]
  map_df[data_est > 800, data_est := 800 ]
  map_df[is_col > 800, is_col := 800 ]
  map_df[oos_col > 800, oos_col := 800 ]

  map_df[, is_data_diff := is_col - data_est ]
  map_df[ is_data_diff < -400, is_data_diff := -400 ]
  map_df[ is_data_diff > 400, is_data_diff := 400 ]
  map_df[, oos_data_diff := oos_col - data_est ]
  map_df[ oos_data_diff < -400, oos_data_diff := -400 ]
  map_df[ oos_data_diff > 400, oos_data_diff := 400 ]

  col_breaks <- seq(0,800,by=200)
  col_breaks_diff <- seq(-400, 400, by=100)
  fig_data <- ggplot(data=map_df) + 
    geom_sf(data=ad0_sf, fill='white', color='black', lwd=.5) +
    geom_point(aes(x=longitude, y=latitude, color=data_est, size=weight)) +
    scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(name='Spectral', n=11)), breaks=col_breaks, limits=range(col_breaks)) +
    labs(title='Raw Data', color='Prevalence\n per 100k', x='', y='')
  fig_is <- ggplot(data=map_df) + 
    geom_sf(data=ad0_sf, fill='white', color='black', lwd=.5) +
    geom_point(aes(x=longitude, y=latitude, color=is_col, size=weight)) +
    scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(name='Spectral', n=11)), breaks=col_breaks, limits=range(col_breaks)) +
    labs(title='In-sample estimates', color='Prevalence\n per 100k', x='', y='')
  fig_oos <- ggplot(data=map_df) + 
    geom_sf(data=ad0_sf, fill='white', color='black', lwd=.5) +
    geom_point(aes(x=longitude, y=latitude, color=oos_col, size=weight)) +
    scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(name='Spectral', n=11)), breaks=col_breaks, limits=range(col_breaks)) +
    labs(title='Out-of-sample estimates', color='Prevalence\n per 100k', x='', y='')
  fig_is_diff <- ggplot(data=map_df) + 
    geom_sf(data=ad0_sf, fill='white', color='black', lwd=.5) +
    geom_point(aes(x=longitude, y=latitude, color=is_data_diff, size=weight)) +
    scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(name='PuOr', n=11)), breaks=col_breaks_diff, limits=range(col_breaks_diff)) +
    labs(title='Diff: IS - Data', color='Diff', x='', y='')
  fig_oos_diff <- ggplot(data=map_df) + 
    geom_sf(data=ad0_sf, fill='white', color='black', lwd=.5) +
    geom_point(aes(x=longitude, y=latitude, color=oos_data_diff, size=weight)) +
    scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(name='PuOr', n=11)), breaks=col_breaks_diff, limits=range(col_breaks_diff)) +
    labs(title='Diff: OOS - Data', color='Diff', x='', y='')

  pdf(fp, height=10, width=25)
  grid.arrange(
    ggplotGrob(fig_data), ggplotGrob(fig_is), ggplotGrob(fig_oos), 
    grid.rect(gp=gpar(col="white")), ggplotGrob(fig_is_diff), ggplotGrob(fig_oos_diff), ncol=3)
  dev.off()
}

make_diag_map(
  oos_predict_full,
  is_col='is_est_median_raked',
  oos_col='oos_est_median_raked',
  fp=sprintf('%s/plot_oos_median_raked_differences.pdf',oos_dir)
)
make_diag_map(
  oos_predict_full,
  is_col='is_est_mean_raked',
  oos_col='oos_est_mean_raked',
  fp=sprintf('%s/plot_oos_mean_raked_differences.pdf',oos_dir)
)

oos_predict_full[is.na(is_est_mean), unique(cluster_number)]
