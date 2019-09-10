## #############################################################################
## 
## ASSEMBLE MODELING DATASET
## 
## Author: Nat Henry
## Created: August 9, 2019
## Purpose: Merge the polygon-resampled point locations with the survey 
##   microdata results in the Pakistan 2010-2011 TB Prevalence Survey
## 
## #############################################################################

## SET INPUTS

# Code filepaths
user <- 'nathenry'
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/',user)

# Year of the prevalence survey
svy_year <- 2011

## Data filepaths
# Survey microdata
kit_dir <- '/ihme/limited_use/LIMITED_USE/LU_GEOSPATIAL/KIT_TB_Project/'
# Input survey microdata
tb_prev_micro_fp <- paste0(
  '/ihme/limited_use/LIMITED_USE/PROJECT_FOLDERS/PAK/NTP_TB_HACKATHON_KIT/',
  '1 Prevalence Survey Cluster Estimates/prevalence_survey_case_based_2010_2011.xlsx'
)
# Point-polygon resampled cluster locations
ptpoly_file <- paste0(kit_dir, '/prep/data_prep/prevalence_survey_point_resampled.csv')
# Output filepaths
out_file <- paste0(kit_dir,'/prep/data_prep/kit_collapsed_resampled.csv')
out_file_mbg <- '/share/geospatial/mbg/input_data/kit_tb_test_dataset_20190809'

## SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load all MBG functions and libraries
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(
  t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)),
  'readxl'
)
mbg_setup(package_list=package_list, repos=core_repo)
source(sprintf('%s/data_central/sae_functions.R',core_repo))


## LOAD AND PREP MICRODATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load microdata file; keep only relevant fields
prev_micro <- readxl::read_excel(tb_prev_micro_fp) %>% as.data.table
keep_fields <- c('cluster','eligible','participant','age','smpostbcase','smear','bact')
prev_micro <- prev_micro[, ..keep_fields]
# Drop non-eligible and non-participating individuals
prev_micro <- prev_micro[(age > 15),]
# Sum total number sampled and total bacteriologically positive by area
prev_micro[, bact_pos := 0]
prev_micro[bact=='Positive', bact_pos := 1]
prev_micro[smpostbcase!='No SS+ TB case', bact_pos := 1]

# Aggregate
prev_agg <- prev_micro[, .(N = .N, kit = sum(bact_pos)), by=cluster]
setnames(prev_agg, 'cluster','cluster_number')


## MERGE ON CLUSTER-LEVEL DATA AND ADD METADATA FIELDS ~~~~~~~~~~~~~~~~~~~~~~~~~

# Load cluster file
clusters <- fread(ptpoly_file)
# Merge data
collapsed_resampled <- merge(
  x = prev_agg,
  y = clusters,
  by = 'cluster_number',
  all.x = TRUE
)[order(cluster_number,-weight)]

## Add metadata fields
collapsed_resampled[, geo_unit := as.character(cluster_number)]
collapsed_resampled[, year := 2011 ] # Keep all one year for now
collapsed_resampled[, svyyr := 2011 ]
collapsed_resampled[, source := 'COUNTRY_SPECIFIC']
collapsed_resampled[, period := 1]
collapsed_resampled[, point := 1 ]
collapsed_resampled[ pseudocluster == TRUE, point := 0]
collapsed_resampled[, nid := 9999999 ] # Dummy NID
collapsed_resampled[, country := 'PAK' ]

# Run a test to make sure that all weights sum to 1
test_data <- collapsed_resampled[, .(weight=sum(weight)), by=cluster_number]
message(paste0(
  "Range of aggregated weights: ",paste(range(test_data$weight), collapse='-')
))

## SAVE THAT DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(collapsed_resampled, file=out_file, row.names=FALSE)
saveRDS(collapsed_resampled, file=out_file_mbg)
