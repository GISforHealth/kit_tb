## #############################################################################
## 
## AUTHOR:  Nat Henry
## CREATED: August 7, 2019
## PURPOSE: Model-based geostatistics launch script for KIT TB Project
## 
## NOTE: This code pulls from many functions in the LBDCore package, which is
##   still in development. Key functions from this package will be added to
##   this repository as a demonstration. This package is on track to be
##   publicly released as an open-source tool.
## 
## #############################################################################


## #############################################################################
## DEFINE INPUTS
## #############################################################################

## Set repo location and indicator group
user <- Sys.info()['user'] # user running this code
core_repo <- sprintf('<<FILEPATH>>')
ig_repo <- sprintf('<<FILEPATH>>')
core_remote <- 'origin'
core_branch <- 'develop'
ig_remote <- 'origin'
ig_branch <- 'master'
indicator_group <- 'tb'
indicator <- 'kit'
region <- reg <- Regions <- regions <- 'pak'

## Create config and covariate config names
config_name <- 'config_kit'
covs_config_name <- 'covs_kit_newcovs_vif2'

# set to true if you will pull into the repo on share before running code
pullgit         <- FALSE  

## IF YOU WANT TO PRESET TO AN OLD RUN DATE DO IT HERE, ELSE SET IT NULL AND 
##  IT WILL SPIN UP A NEW RUN DATE(see run_date_notes.csv)
preset_run_date <- '2019_09_03_binom_for_submission'

## sort some directory stuff and pull newest code into share
if(pullgit) {
  system(sprintf('cd %s\ngit pull %s %s', core_repo, core_remote, core_branch))
  system(sprintf('cd %s\ngit pull %s %s', ig_repo,   ig_remote,   ig_branch))
}

## Create run date in correct format
if(is.null(preset_run_date)){
  run_date <- make_time_stamp(TRUE)
  message(sprintf('New run_date: %s', run_date))
} else {
  run_date <- preset_run_date
  message(sprintf('Preset run_date: %s', run_date))
}
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
dir.create(sprintf('%s/output/%s', sharedir, run_date))

## Copy config file and covariates config file to the run directory
message("Copying config and covariate config files to the run directory.")
file.copy(
  from = sprintf('%s/%s.csv', ig_repo, config_name), 
  to = sprintf('%s/output/%s/config.csv', sharedir, run_date), 
  overwrite = TRUE
)
file.copy(
  from = sprintf('%s/%s.csv', ig_repo, covs_config_name), 
  to = sprintf('%s/output/%s/covariates_config.csv', sharedir, run_date), 
  overwrite = TRUE
)


## #############################################################################
## RUN MBG PARALLEL SCRIPT
## #############################################################################

reg         <- regions
age         <- 0
holdout     <- 0
samples <- as.numeric(samples)

source(sprintf('%s/kit_parallel_model.R',ig_repo))

source(sprintf('%s/kit_postest.R',ig_repo))

source(sprintf('%s/kit_oos.R',ig_repo))

message('~~~FIN~~~')
