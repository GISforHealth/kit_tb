## #############################################################################
## 
## POINT RESAMPLE PAKISTAN ADMIN3 DATA
## 
## Author: Nat Henry
## Created: August 6, 2019
## Purpose: For locations that were matched to the Admin3 level but not the
##   admin4 level, point resample the data by taking the population-weighted
##   centroids of each admin4 and splitting admin3 observations into admin4
##   estimates proportional to the population in each admin4.
## 
## #############################################################################

## SET INPUTS

# Code filepaths
user <- 'nathenry'
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/',user)

# Year of the prevalence survey
svy_year <- 2011

# Data filepaths
kit_dir <- '/ihme/limited_use/LIMITED_USE/LU_GEOSPATIAL/KIT_TB_Project/'
geomatching_fp <- paste0(
  kit_dir,'/prep/raw_data/prevalence_survey_geomatched_20190806.csv'
)

# Shapefile filepaths and identifying codes
ad3_shp_fp <- paste0(
  kit_dir,'/shapefiles/HDX_Shp_Admin2/pak_adm_ocha_pco_gaul_20181218_SHP/',
  'pak_admbnda_adm3_ocha_pco_gaul_20181218.shp'
)
ad3_id <- 'ADM3_PCODE'
ad4_shp_fp <- paste0(
  kit_dir,'shapefiles/HDX_Shp_Admin3/Adminbdy Shapefile/Union_Council.shp'
)
ad4_id <- 'UC_C'

# Link table filepath (to be created and saved if the file does not exist)
ad4_link_fp <- paste0(
  kit_dir,'shapefiles/HDX_Shp_Admin3/Adminbdy Shapefile/Union_Council_link.RDS'
)

# Output file for geocoding
out_file <- paste0(kit_dir, '/prep/data_prep/prevalence_survey_point_resampled.csv')


## SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load all MBG functions and libraries
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)
source(sprintf('%s/data_central/sae_functions.R',core_repo))


## LOAD AND PREP SHAPEFILES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load admin3 shapefile
ad3_sp <- rgdal::readOGR(ad3_shp_fp)
ad3_sp$ad3_id <- ad3_sp@data[[ad3_id]]
ad3_sp@data <- ad3_sp@data[, c('ADM3_EN','ad3_id')]

## Load admin4 shapefile
ad4_sp <- rgdal::readOGR(ad4_shp_fp)

## Load admin4 link table
if(file.exists(ad4_link_fp)){
  # Load from file
  ad4_link_list <- readRDS(ad4_link_fp)
} else {
  # If a link table does not already exist, build one and save to file
  ad4_link_list <- build_link_table(
    shapefile_version=NULL,
    cores=get_max_forked_threads(nobjs=10),
    region=NULL,
    custom_shapefile_path=ad4_shp_fp,
    custom_shapefile_field=ad4_id
  )
  saveRDS(ad4_link_list, file=ad4_link_fp)
}
ad4_link <- ad4_link_list$link_table

## ADMIN3 POPULATION CENTROID DATASET
## Merge on this only in cases where no admin4s overlap an admin3
ad3_latlong <- frac_agg_covs(
  covs=c('lat','lon'), years=svy_year, measures=c('mean','mean'),
  releases=c('2019_06_10','2019_06_10'), shapefile_path=ad3_shp_fp,
  shapefile_field=ad3_id, core_repo=core_repo, agg_method="pop_weight",
  worldpop_age_sex_group="total", worldpop_age_sex_release="2017_04_27",
  cores=get_max_forked_threads(nobjs=10), link_table=NULL
)
ad3_latlong$year <- NULL
setnames(ad3_latlong, c(ad3_id,'lat','lon'), c('ad3_id','latitude','longitude'))
ad3_latlong$weight <- 1

## ADMIN4 POPULATION CENTROID DATASET
## Create a latitude-longitude SpatialPointsDataFrame containing population info
# Get population-weighted centroids for each UC
ad4_latlong <- frac_agg_covs(
  covs=c('lat','lon'), years=svy_year, measures=c('mean','mean'),
  releases=c('2019_06_10','2019_06_10'), shapefile_path=ad4_shp_fp,
  shapefile_field=ad4_id, core_repo=core_repo, agg_method="pop_weight",
  worldpop_age_sex_group="total", worldpop_age_sex_release="2017_04_27",
  cores=get_max_forked_threads(nobjs=10), link_table=ad4_link
)
# Get total population for each UC
ad4_population <- frac_agg_covs(
  covs=c('worldpop'), years=svy_year, measures=c('total'),
  releases=c('2017_04_27'), shapefile_path=ad4_shp_fp,
  shapefile_field=ad4_id, core_repo=core_repo, agg_method="sum",
  cores=get_max_forked_threads(nobjs=10), link_table=ad4_link
)
# Remove erroneous code, merge, and convert to SPDF
ad4_data_full <- merge(x=ad4_latlong, y=ad4_population, by=c(ad4_id, 'year'))
setnames(ad4_data_full, c(ad4_id, 'total'), c('admin4_code','pop'))
ad4_data_full <- ad4_data_full[admin4_code > 0,]
ad4_points <- SpatialPointsDataFrame(
  coords=ad4_data_full[, .(lon, lat)],
  proj4string=CRS(proj4string(ad3_sp)),
  data=ad4_data_full[, .(admin4_code, pop)]
)

## Check overlap between admin3 codes and admin4 codes
admin_overlap <- as.data.table(cbind(
  ad4_data_full,
  sp::over(x=ad4_points, y=ad3_sp)
))
# Observations will be split based on the proportion of population in each 
#  admin4 contained within an admin3 location
admin_overlap[, total_pop_in_ad3 := sum(pop), by=ad3_id]
admin_overlap[, ad4_split_proportion := pop / total_pop_in_ad3]


## SPLIT GEOMATCHING DATASET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load and prep geomatching dataset
match_dt <- fread(geomatching_fp)
# Column formatting
names(match_dt) <- tolower(gsub(' ','_',names(match_dt)))
keep_fields <- c('cluster_number','adm3_pcode','confidence','lat','long')
match_dt <- match_dt[, ..keep_fields]
setnames(match_dt, c('adm3_pcode','lat','long'), c('ad3_id','latitude','longitude'))
# Split into two datasets containing "true" points and geomatched points
# NOTE: This could also be done on the basis of geomatching confidence
match_true_points <- match_dt[!is.na(latitude) & !is.na(longitude), ]
match_need_split <- match_dt[is.na(latitude) | is.na(longitude), ]
# Add metadata for true points
match_true_points[, weight := 1 ]
match_true_points[, shapefile := NA ]
match_true_points[, location_code := NA ]
match_true_points[, pseudocluster := FALSE ]

## Prepare admin overlap data table for merge
ad4_for_split <- admin_overlap[, 
  .(ad3_id, admin4_code, lat, lon, ad4_split_proportion)
]
setnames(
  ad4_for_split, 
  c('lat','lon','admin4_code','ad4_split_proportion'), 
  c('latitude','longitude','location_code','weight')
)

## Merge geomatching dataset onto the splitting dataset
match_need_split[, c('latitude', 'longitude') := NULL]
match_split <- merge(
  x=match_need_split, y=ad4_for_split, by='ad3_id', all.x=TRUE
)
# Split out any unsuccessfull matches
need_admin3_centroid <- match_split[is.na(weight)]
match_split <- match_split[!is.na(weight)]
# Add metadata
match_split[, shapefile := 'HDX_Shp_Admin4']
match_split[, pseudocluster := TRUE ]

## Merge on admin3 centroid in case admin4 match was unsuccessfull
message(sprintf(
  "%i locations were missing admin4 overlap: %s", nrow(need_admin3_centroid),
  paste0(need_admin3_centroid$ad3_id, collapse=', ')
))
need_admin3_centroid[, c('location_code','latitude','longitude','weight'):=NULL]
admin3_centroid <- merge(
  x=need_admin3_centroid, y=ad3_latlong, by=c('ad3_id'), all.x=TRUE
)
# If the admin3 centroid approach does not work, stop the program
admin3_merge_fail <- admin3_centroid[is.na(weight),]
if(nrow(admin3_merge_fail)>0){
  stop(sprintf(
    "COULD NOT FIND ADMIN3: %s",paste0(admin3_merge_fail$ad3_id, collapse=', ')
  ))
}
# Add admin3 merge metadata
admin3_centroid[, shapefile := 'HDX_Shp_Admin3']
admin3_centroid[, location_code := ad3_id ]
admin3_centroid[, pseudocluster := TRUE ]
admin3_centroid[, weight := 1 ] # Only 'split' into a single point

## MERGE ALL 3 SPLIT TYPES TOGETHER AND SAVE
geocoding_combined <- rbindlist(
  list(match_true_points, match_split, admin3_centroid),
  use.names=TRUE
)[order(cluster_number)]
# Drop admin3 code (duplicative)
geocoding_combined[, ad3_id := NULL ]
# Save
write.csv(geocoding_combined, file=out_file, row.names=FALSE)
