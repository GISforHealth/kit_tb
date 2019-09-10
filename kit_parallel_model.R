################################################################################
## Generic parallel script for running model-based geostatistics              ##
##                                                                            ##
## Created by Roy Burstein, Nick Graetz, Aaron Osgood-Zimmerman, Jon Mosser   ##
## Adapted by Nat Henry for KIT TB competition                                ##
## Adapted: August 10, 2019                                                   ##
##                                                                            ##
## NOTE: This code pulls from many functions in the LBDCore package, which is ##
##   still in development. Key functions from this package will be added to   ##
##   this repository as a demonstration. This package is on track to be       ##
##   publicly released as an open-source tool.                                ##
##                                                                            ##
################################################################################


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ## make a pathaddin that get used widely
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PID <- Sys.getpid()
tic("Entire script") # Start master timer

## skip a large chunk if requested in config
if(as.logical(skiptoinla) == FALSE){

  message('You have chosen to not skip directly to inla.')

  ## some set up
  if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
  if (class(year_list_predict) == "character") year_list_predict <- eval(parse(text=year_list_predict))
  if (class(z_list)    == "character") z_list    <- eval(parse(text=z_list))

  ## Load simple polygon template to model over
  gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(
    gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
    custom_shapefile_path = '<<FILEPATH>>'
  )
  subset_shape        <- simple_polygon_list[[1]]
  simple_polygon      <- simple_polygon_list[[2]]

  ## Load list of raster inputs (pop and simple)
  raster_list        <- build_simple_raster_pop(subset_shape, link_table=NULL)
  simple_raster      <- raster_list[['simple_raster']]
  pop_raster         <- raster_list[['pop_raster']]

  ## Load input data based on stratification and holdout, OR pull in data as 
  ## normal and run with the whole dataset if holdout == 0. For holdouts, we 
  ## have depreciated val level, so each val level must be recorded in a different run date
  message('Holdout == 0 so loading in full dataset using load_input_data()')
  df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                        simple      = simple_polygon,
                        agebin      = age,
                        removeyemen = FALSE,
                        pathaddin   = pathaddin,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share),
                        yl          = year_list)
  if(holdout!=0) {
    message(paste0('Holdout != 0 so loading holdout data only from holdout ',holdout))
    nids_to_drop <- oos_fold_table[fold==holdout,]$nid
    df <- subset(df,! nid %in% nids_to_drop)
    write.csv(df,sprintf('%s/input_data_%s.csv',outputdir,pathaddin))
  }
  
  ## if testing, we only keep 1000 or so observations
  if(test == 1){
    test_pct <- as.numeric(test_pct)

    message(paste0('Test option was set on and the test_pct argument was found at ',test_pct,'% \n
                 ... keeping only ', round(nrow(df)*(test_pct/100),0),' random rows of data.'))
    df <- df[sample(nrow(df),round(nrow(df)*(test_pct/100),0)),]

    message('Also, making it so we only take 100 draws')
    samples <- 100
  }

  ## if there is another weight column, multiply it with weight now
  if(exists('other_weight')) if(other_weight!='') {
                               message(paste0('Multiplying weight and ',other_weight))
                               df[['weight']] <- df[['weight']]*df[[other_weight]]
                             }

  ## Some built in data checks that cause known problems later on
  if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
    stop('You have binomial data where k > N. Check your data before proceeding')
  if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
    stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')

  ## Save distribution of data for this region
  png(paste0(outputdir, reg, '.png'))
  if(indicator_family=='binomial') hist(df[, get(indicator)]/df$N) else hist(df[, get(indicator)])
  dev.off()

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Define modeling space. In years only for now.
  period_map <- make_period_map(
    modeling_periods = year_list_predict
  )

  ## Make placeholders for covariates
  cov_layers <- gbd_cov_layers <- NULL

  ## Pull all covariate bricks/layers
  # Use a separate process to pull the poverty covariate
  custom_covars_list <- c(
    'wp_poverty', 'acled_protests', 'acled_violence', 'tb_facil_density'
  )
  custom_covars <- custom_covars_list[
    custom_covars_list %in% fixed_effects_config$covariate
  ]
  fixed_effects_no_custom <- fixed_effects_config[!(covariate %in% custom_covars_list),]
  if (nrow(fixed_effects_no_custom) > 0) {
    message('Grabbing raster covariate layers')
    loader <- MbgStandardCovariateLoader$new(start_year = min(year_list_predict),
                                             end_year = max(year_list_predict),
                                             interval = as.numeric(interval_mo),
                                             covariate_config = fixed_effects_no_custom)
    cov_layers <- loader$get_covariates(simple_polygon)
  }
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Pull custom covariates, if specified
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  custom_templ_rr <- simple_raster
  kit_asset_dir <- '<<FILEPATH>>'
  ## Poverty 
  if('wp_poverty' %in% custom_covars){
    message("Loading wp_poverty which is a synoptic **custom covariate**")
    # Load poverty raster and extend to the template
    pov_r <- raster::raster(paste0(kit_asset_dir,'/wp_poverty_pak_resampled_filled.tif'))
    pov_r <- extend(pov_r, y=custom_templ_rr)
    pov_r <- crop(pov_r, y=custom_templ_rr)
    if(any(dim(pov_r) != dim(custom_templ_rr))) stop("Issue with poverty raster dimensions.")
    # Convert to raster brick
    pov_rr <- raster::brick(lapply(year_list_predict, function(x) pov_r))
    # Append to the covariate list
    cov_layers <- c(cov_layers, list(wp_poverty=pov_rr))
  }
  ## ACLED: Protests
  if('acled_protests' %in% custom_covars){
    message("Loading acled_protests which is an annual **custom covariate**")
    # Load, extend, subset
    acled_p <- raster::brick(paste0(kit_asset_dir,'/acled_protest_covar.tif'))
    acled_p <- crop(extend(acled_p, y=custom_templ_rr), y=custom_templ_rr)
    if(any(dim(acled_p)[1:2] != dim(custom_templ_rr)[1:2])) stop("Issue with acled_p.")
    # Append to covariate list
    cov_layers <- c(cov_layers, list(acled_protests=acled_p))
  }
  ## ACLED: Violence
  if('acled_violence' %in% custom_covars){
    message("Loading acled_violence which is an annual **custom covariate**")
    # Load, extend, subset
    acled_v <- raster::brick(paste0(kit_asset_dir,'/acled_violence_covar.tif'))
    acled_v <- crop(extend(acled_v, y=custom_templ_rr), y=custom_templ_rr)
    if(any(dim(acled_v)[1:2] != dim(custom_templ_rr)[1:2])) stop("Issue with acled_v.")
    # Append to covariate list
    cov_layers <- c(cov_layers, list(acled_violence=acled_v))
  }
  ## TB facility density
  if('tb_facil_density' %in% custom_covars){
    message("Loading tb_facil_density which is a synoptic **custom covariate**")
    # Load, extend, convert to raster brick
    tb_density <- raster::brick(paste0(kit_asset_dir,'/facils_per_100k_covar.tif'))
    tb_density <- crop(extend(tb_density, y=custom_templ_rr), y=custom_templ_rr)
    if(any(dim(tb_density)[1:2] != dim(custom_templ_rr)[1:2])) stop("Issue with tb_density.")
    tb_rr <- raster::brick(lapply(year_list_predict, function(x) tb_density))
    # Append to covariate list
    cov_layers <- c(cov_layers, list(tb_facil_density=tb_rr))
  }

  ## Pull country level gbd covariates
  if (nchar(gbd_fixed_effects) > 0) {
    message('Grabbing GBD covariates')

    effects <- trim(strsplit(gbd_fixed_effects, "\\+")[[1]])
    measures <- trim(strsplit(gbd_fixed_effects_measures, "\\+")[[1]])
    gbd_cov_layers <- load_gbd_covariates(covs     = effects,
                                          measures = measures,
                                          year_ids = year_list,
                                          age_ids  = gbd_fixed_effects_age,
                                          template = cov_layers[[1]][[1]],
                                          simple_polygon = simple_polygon,
                                          interval_mo = interval_mo)
  }

  ## Combine all covariates
  all_cov_layers <- c(cov_layers, gbd_cov_layers)

  ## regenerate all fixed effects equation from the cov layers
  all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")

  ## Make stacker-specific formulas where applicable
  all_fixed_effects_brt <- all_fixed_effects


  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Stacking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  tic("Stacking - all") ## Start stacking master timer

  ## Figure out which models we're going to use
  child_model_names <- stacked_fixed_effects        %>%
                          gsub(" ", "", .)          %>%
                          strsplit(., "+", fixed=T) %>%
                          unlist
  message(paste0('Child stackers included are: ',paste(child_model_names,collapse=' // ')))

  the_covs <- format_covariates(all_fixed_effects)

  ## copy the dataset to avoid unintended namespace conflicts
  the_data <- copy(df)

  ## shuffle the data into six folds
  the_data <- the_data[sample(nrow(the_data)),]
  the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=as.numeric(n_stack_folds),labels=FALSE)]

  ## add a row id column
  the_data[, a_rowid := seq(1:nrow(the_data))]

  ## extract covariates to the points and subset data where its missing covariate values
  cs_covs <- extract_covariates(the_data,
                                all_cov_layers,
                                id_col              = "a_rowid",
                                return_only_results = TRUE,
                                reconcile_timevarying = TRUE,
                                centre_scale        = TRUE,
                                period_var          = 'year',
                                period_map          = period_map)

  # A check to see if any of the variables do not vary across the data. 
  #  This could break model later so we check and update some objects
  covchecklist <- check_for_cov_issues(
    check_pixelcount = FALSE,
    check_pixelcount_thresh = ifelse(
      exists("pixelcount_thresh"), as.numeric(pixelcount_thresh), 0.95
    )
  )
  for(n in names(covchecklist)){
    assign(n, covchecklist[[n]])
  }

  # plot covariates as a simple diagnostic here
  pdf(sprintf('%s/raw_covariates_%s.pdf',outputdir,pathaddin), height=12, width=12)
  for(covname in names(all_cov_layers)){
    plot(all_cov_layers[[covname]],main=covname,maxpixel=1e6)
  }
  dev.off()

  ## Check for data where covariate extraction failed
  rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
  if (rows_missing_covs > 0) {
    pct_missing_covs <- round((rows_missing_covs/nrow(the_data))*100, 2)
    warning(paste0(rows_missing_covs, " out of ", nrow(the_data), " rows of data ",
                   "(", pct_missing_covs, "%) do not have corresponding ",
                   "covariate values and will be dropped from child models..."))
    if (rows_missing_covs/nrow(the_data) > 0.1) {
      stop(paste0("Something has gone quite wrong: more than 10% of your data does not have ",
                  "corresponding covariates.  You should investigate this before proceeding."))
    }
  }

  the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)

  ## store the centre scaling mapping
  covs_cs_df  <-  cs_covs[[2]]

  ## this will drop rows with NA covariate values
  the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))

  ## stop if this na omit demolished the whole dataset
  if(nrow(the_data) == 0) stop(
    'You have an empty df, make sure one of your covariates was not NA everywhere.'
  )

  ## seperated out into a different script
  if(as.logical(use_stacking_covs)){
    message('Fitting Stackers')
    
    # Run the child stacker models 
    child_model_run <- run_child_stackers(models = child_model_names, input_data = the_data)
    
    # Bind the list of predictions into a data frame
    child_mods_df <- do.call(cbind, lapply(child_model_run, function(x) x[[1]]))
    
    ## combine the children models with the_data
    the_data  <- cbind(the_data, child_mods_df)
    
    ## Rename the child model objects into a named list
    child_model_objs <- setNames(lapply(child_model_run, function(x) x[[2]]), child_model_names)
    
    
    
    ## return the stacked rasters
    stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                          period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                          child_models     = child_model_objs,
                                          indicator_family = indicator_family,
                                          centre_scale_df  = covs_cs_df)

    ## plot stackers
    pdf(paste0(outputdir, 'stacker_rasters', pathaddin, '.pdf'))
    for(i in 1:length(stacked_rasters))
      plot(stacked_rasters[[i]],main=names(stacked_rasters[[i]]),maxpixel=ncell(stacked_rasters[[i]]))
    dev.off()

    message('Stacking is complete')
  } ## if(use_stacking_covs)


  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## set the fixed effects to use in INLA based on config args
  if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- ''
  }
  if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- stacked_fixed_effects
  }
  if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- fixed_effects ## from config
  }
  if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
  }
  if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, sep = " + ")
  }
  if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }
  if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }
  if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }

  ## copy things back over to df
  df <- copy(the_data)

  ## remove the covariate columns so that there are no name conflicts when they get added back in
  df <- df[,paste0(the_covs) := rep(NULL, length(the_covs))]

  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
  df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]

  ## create a full raster list to carry though to the shiny/next steps
  if(as.logical(use_stacking_covs)){
    cov_list      <- c(unlist(stacked_rasters),unlist(all_cov_layers))
    child_mod_ras <- cov_list[child_model_names]
  }else{
    cov_list <- unlist(all_cov_layers)
    child_model_names <- ''
  }

  toc(log = T) ## End stacking master timer

  ## make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

  ## Build spatial mesh over modeling area
  message(paste0("Creating spatial mesh, max edge parameter: ", mesh_s_max_edge))
  mesh_s_max_edge <- eval(parse(text = mesh_s_max_edge))
  mesh_s_offset <- eval(parse(text = mesh_s_offset))
  mesh_s <- inla.mesh.2d(
    boundary = inla.sp2segment(simple_polygon),
    loc = cbind(df$longitude, df$latitude),
    max.edge = mesh_s_max_edge,
    offset = mesh_s_offset,
    max.n.strict = 40,
    cutoff = .5
  )
  # Plot spatial mesh
  pdf(sprintf('%s/spatial_mesh_with_data.pdf',outputdir), height=8, width=8)
  plot(mesh_s, asp = 1)
  lines(subset_shape)
  points(df$longitude, df$latitude, col = df$year)
  dev.off()


  ## Build temporal mesh (standard for now)
  if (length(unique(year_list)) == 1) {
    mesh_t <- NULL
  } else { 
    mesh_t <- build_time_mesh(periods = eval(parse(text = mesh_t_knots)))
  }


  ## For raw covs, don't want to center-scale (as that will happen in `build_mbg_data_stack()`)
  ##
  ## For predict_mbg, non-center-scaled covs are pulled from cov_list (either stackes or raw) and
  ## center-scaled within the function.  So both fit and predict take place on center-scaled covs
  if (as.logical(use_raw_covs) == TRUE) {
    centre_scale_covs <- FALSE
  } else {
    centre_scale_covs <- TRUE
  }

  ## Save all inputs for MBG model into correct location on /share
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 mesh_s            = mesh_s,
                 mesh_t            = mesh_t,
                 cov_list          = cov_list,
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = centre_scale_covs)

} else { ## END !SKIPTOINLA
  message(paste0('You have chosen to skip directly to INLA. Picking up objects from run_date ',skiptoinla_from_rundate))
  message('Now copying saved MBG inputs from that chosen run_date.')

  file.copy(from = '<<FILEPATH>>', to = '<<FILEPATH>>')
}

## reload data an prepare for MBG
load('<<FILEPATH>>')

# Bound GBM to 0-1 if desired
if (exists("gbm_bounded_0_1")) {
  if (as.logical(gbm_bounded_0_1) == T) {
    message("Truncating GBM values > 1 to 0.999")
    values(cov_list[["gbm"]])[values(cov_list[["gbm"]]) >= 1 & !is.na(values(cov_list[["gbm"]]))] <- 0.999
    gbm_cols <- grep(paste0("(gbm)(.*_pred)"), names(df), value=T)
    replace_one <- function(x) {
      x[x>=1 & !is.na(x)] <- 0.999
      return(x)
    }
    df[, (gbm_cols) := lapply(.SD, replace_one), .SDcols = gbm_cols]
  }
}

## convert stackers to transform space, if desired
## NOTE: we do this here to ensure that the stacker rasters are saved in prevalence/untransformed space
## this is useful for diagnostics and other code that was built expecting the untransformed rasters
if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)){
  message('Converting stackers to logit space')

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


tic("MBG - all") ## Start MBG master timer

input_data <- build_mbg_data_stack_tmb(
  d = df,
  fes = all_fixed_effects,
  indic = indicator, 
  mesh = mesh_s,
  cov_constraints = covariate_constraint_vectorize(config)
)

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

tic("MBG - fit model") ## Start MBG - model fit timer

## Set the number of cores to be equal to input;
## If missing, then revert to cores_to_use value
if(Sys.getenv("OMP_NUM_THREADS") != "") {
    setompthreads(Sys.getenv("OMP_NUM_THREADS"))
} else {
    print("Threading information not found; setting cores_to_use as the input OpenMP threads.")
    setompthreads(cores_to_use)
}

if(Sys.getenv("MKL_NUM_THREADS") != "") {
    setmklthreads(Sys.getenv("MKL_NUM_THREADS"))
} else {
    print("Threading information not found; setting cores_to_use as the input MKL threads.")
    setmklthreads(cores_to_use)
}

## Fit MBG model
message('Fitting model with TMB')
message(sprintf('%i Data points and %i mesh nodes',nrow(df),length(input_data$Parameters$epsilon_s)))

# save RDS file of input data for replication
saveRDS(object = input_data, ## save this here in case predict dies
        file = '<<FILEPATH>>')
# run the model
system.time(
  model_fit <- fit_mbg_tmb(
    lbdcorerepo = core_repo,
    cpp_template = 'mbg_tmb_model',
    tmb_input_stack = input_data,
    control_list = list(trace=1, maxit=1000, itnmax=1000, eval.max=1000, iter.max=1000, abstol=-1E10, reltol=1E-8),
    optimizer = 'optimx',
    ADmap_list = NULL,
    newton_steps = 0,
    sparse_ordering = as.logical(sparse_ordering)
  )
)

# clamping
clamp_covs <- TRUE

saveRDS(object = model_fit, ## save this here in case predict dies
        file = '<<FILEPATH>>')

tic("MBG - predict model") ## Start MBG - model predict timer

## Run predict_mbg on chunks of 50 samples (to avoid memory issues)
message('Making predictions in 50 draw chunks.')

max_chunk <- 50
samples   <- as.numeric(samples)

if(is.character(year_list_predict)) year_list_predict <- eval(parse(text=year_list_predict))

pred_transform <- ifelse(indicator_family=='poisson','exp','inverse-logit')

chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
pm <- lapply(chunks, function(samp) {
  predict_mbg_tmb(
    samples          = samp,
    seed             = NULL,
    tmb_input_stack  = input_data,
    model_fit_object = model_fit,
    fes              = all_fixed_effects,
    sr               = simple_raster,
    yl               = year_list_predict,
    zl               = z_list,
    covs_list        = cov_list,
    clamp_covs       = clamp_covs,
    transform        = pred_transform,
    cov_constraints  = covariate_constraint_vectorize(config)
  )
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Make cell preds and a mean raster
cell_pred <- do.call(cbind, pm)
mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))
toc(log = T) # Stop MBG - model predict timer

message('Wrapping up')
save_mbg_preds(config     = config,
               time_stamp = time_stamp,
               run_date   = run_date,
               mean_ras   = mean_ras,
               sd_ras     = NULL,
               res_fit    = model_fit,
               cell_pred  = cell_pred,
               df         = df,
               pathaddin  = pathaddin)


# plot the mean raster
pdf(paste0(outputdir,'mean_rasterXX', pathaddin, '.pdf'))
plot(mean_ras,maxpixel=1e6)
dev.off()

# plot the mean raster WITH DATA OVERLAYED
pts_for_plot <- SpatialPoints(cbind(the_data$longitude, the_data$latitude))
pdf(paste0(outputdir,'mean_rasterXX_with_data', pathaddin, '.pdf'))
plot(mean_ras[[2]],maxpixel=1e6)
points(pts_for_plot, col='blue')
dev.off()

pdf(paste0(outputdir,'mesh_with_data', pathaddin, '.pdf'))
plot(mesh_s)
points(pts_for_plot, col='blue')
dev.off()

## timer stuff
toc(log = T) # End master timer

## Format timer
ticlog   <- tic.log(format = F)
df_timer <- generate_time_log(ticlog)
df_timer[, region := reg]
df_timer[, holdout := holdout]
setcolorder(df_timer, c("region", "holdout", "step", "time"))

## Pull run time for this run
run_time_all <- df_timer[step == "Entire script", time]

## Write to a run summary csv file in the output directory
output_file <- paste0(outputdir, "run_summary_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}

# Write a an empty file to indicate done with this parallel script
write(NULL, file = paste0(outputdir, "/fin_", pathaddin))

## Make a nice fitted parameter summary table
param_table <- fitted_param_table_tmb(
  model_fit=model_fit,
  exp_fixed_effects     = TRUE,
  transform_hyperparams = TRUE,
  draws                 = samples,
  calculate_range       = TRUE,
  calculate_nom_var     = TRUE,
  cov_constraints = covariate_constraint_vectorize(config)
)
write.csv(
  param_table, file=paste0(outputdir,'/fitted_parameter_summary_table.csv'),
  row.names=FALSE
)

# Plot parameter table
param_fig <- ggplot(
      param_table[!(param_name %in% c('int', 'kappa','range','tau','nominal_variance')),],
      aes(x=param_name, y=median)
    ) +
  geom_hline(yintercept=1, color='black', linetype=3) + 
  ylab('Median effect size and 95% UI') +
  xlab("Fixed Effect Coefficients for Spatial Covariates and Hyperparameters") +
  geom_pointrange(aes(ymin=lower, ymax=upper), position=position_dodge(width=.5)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14))

viz_dir <- sprintf('%s/visualizations/',outputdir)
dir.create(viz_dir, showWarnings=FALSE)
png(sprintf('%s/coefficients_across_regions.png',viz_dir),height=600,width=1000)
print(param_fig)
dev.off()

## Write CSV
write.csv(df_timer,file = output_file, row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
