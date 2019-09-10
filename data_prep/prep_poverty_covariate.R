## #############################################################################
## 
## PREP POVERTY COVARIATE
## Author: Nat Henry
## Created: August 21, 2019
## Purpose: Prep poverty covariate
## 
## #############################################################################

## SET INPUTS 

core_repo <- '<<FILEPATH>>'
asset_dir <- '<<FILEPATH>>'
raw_wp_dir <- '<<FILEPATH>>'
modeling_shapefile_version <- '2019_05_06'

## SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load all MBG functions and libraries
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)


## LOAD AND FORMAT INPUT DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load (or create) Worldpop total 1k raster
wp_cache_file <- sprintf('%s/wp_pakistan_1k_2007.tif',asset_dir)
if(file.exists(wp_cache_file)){
  # Load from cache
  full_wp_raster <- raster::raster(wp_cache_file)
} else {
  # Create raster
  wp_sub_dirs <- grep('//[m|f]_', list.dirs(raw_wp_dir), value=TRUE)
  full_wp_raster <- Reduce('+', lapply(wp_sub_dirs, function(wp_dir){
    # Retrieve age group
    ag <- gsub('/','',gsub(raw_wp_dir,'',wp_dir))
    message(sprintf("Pulling raster for %s...",ag))
    return(raster::raster(sprintf('%s/pak_%s_2007_agg.tif',wp_dir,ag)))
  }))
  # Save to cache file
  writeRaster(full_wp_raster, file=wp_cache_file, format="GTiff", overwrite=TRUE)
}

## Load Pakistan poverty data
wp_pov_pak <- raster::raster(sprintf('%s/pak_1k_poverty/pak07povmpi.tif',asset_dir))
# Check that the dimensions match the worldpop population raster
if(any(dim(wp_pov_pak) != dim(full_wp_raster))) stop("Raster dimension mismatch")

## Load simple raster based on Pakistan admin2 shapefile
gaul_list           <- get_adm0_codes('pak')
subset_shape <- load_simple_polygon(
  gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
  custom_shapefile_path = sprintf('%s/pak_ad2.shp',asset_dir)
)[[1]]
simple_raster <- build_simple_raster_pop(
  subset_shape, link_table=NULL
)$simple_raster

## Plot raster BEFORE resample
pdf(sprintf('%s/wp_poverty_pak_original.pdf',asset_dir), height=10, width=10)
raster::plot(wp_pov_pak, col=viridis::viridis(n=100), main='Original (1k x 1k)')
lines(subset_shape)
dev.off()

## Create a new cropped simple raster matching 1k raster extents
cropped_sr <- crop(simple_raster, y=wp_pov_pak)
cropped_sr <- extend(simple_raster, y=wp_pov_pak)
## Add unique identifiers for each non-NA cell
idx_vals <- cellIdx(cropped_sr)
index_sr <- insertRaster(cropped_sr, new_vals=matrix(idx_vals,ncol=1))


## Create population-weighted aggregates for poverty ~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Resample index_sr to match dimensions of poverty raster
index_sr_resample <- raster::resample(
  x=index_sr, y=wp_pov_pak, method='ngb', # Nearest neighbor
)
# Create data.table using index values, then merge on cell_idx to ensure that
#  all values were included
pov_data_full <- data.table(
  pov = as.vector(wp_pov_pak),
  pop = as.vector(full_wp_raster),
  id = as.vector(index_sr_resample)
)
pov_data_full <- merge(
  x=data.table(id=idx_vals),
  y=pov_data_full,
  by='id',
  all.x=TRUE
)
pov_data_agg <- pov_data_full[, 
  .(pov=weighted.mean(pov,w=pop,na.rm=T), pop=sum(pop, na.rm=T)), by=id
][order(id)]

## Insert back into the simple_raster and resize to the correct value
if(any(pov_data_agg$id != idx_vals)) stop('Issue with simple_raster ordering')
pov_resampled <- insertRaster(cropped_sr, new_vals=matrix(pov_data_agg$pov,ncol=1))
pop_resampled <- insertRaster(cropped_sr, new_vals=matrix(pov_data_agg$pop,ncol=1))
# Extend
pov_resampled <- extend(pov_resampled, simple_raster)
pop_resampled <- extend(pop_resampled, simple_raster)

## SAVE IT!
writeRaster(
  pov_resampled, format="GTiff", overwrite=TRUE,
  file=sprintf('%s/wp_poverty_pak_resampled.tif',asset_dir)
)

## Plot resampled poverty raster for vetting purposes
pdf(sprintf('%s/wp_poverty_pak_resampled.pdf',asset_dir), height=10, width=10)
raster::plot(pov_resampled, col=viridis::viridis(n=100), main='Resampled')
lines(subset_shape)
dev.off()


## FILL MISSING RASTER AREAS WITH MEAN POVERTY BY ADMIN2 ~~~~~~~~~~~~~~~~~~~~~~~

# Rasterize subset_shape to extent of simple_raster
r_template_full <- rasterize(subset_shape, y=simple_raster, field='ADM2_CODE')
# Combine admin2 codes with pixel-level estimates of poverty and population
ad2_fill_dt <- data.table(
  pov = as.vector(pov_resampled),
  pop = as.vector(pop_resampled),
  ADM2_CODE = as.vector(r_template_full)
)[!is.na(ADM2_CODE)]
ad2_fill_dt[, sort_order := .I ]
# Get mean poverty by admin2 (ensuring that no admin2 remains blank)
pov_by_ad2 <- ad2_fill_dt[, .(pov=weighted.mean(pov,w=pop,na.rm=T)),by=ADM2_CODE]
pov_by_ad2 <- merge(
  x = data.table(ADM2_CODE=subset_shape$ADM2_CODE),
  y = pov_by_ad2,
  by='ADM2_CODE',
  all.x=TRUE
)
## ** Fill India-administered Kashmir with data from a separate source **
# Source: "Table 162, Number and Percentage of Population Below Poverty Line". 
#   Reserve Bank of India, Government of India. 2013. Retrieved April 20, 2014.
pov_by_ad2[ADM2_CODE==0, pov:=.1035]
# Ensure that poverty is filled in all admin2 units
if(any(is.na(pov_by_ad2$pov))){
  stop("Some districts have not been a assigned poverty value.")
}

# Merge data back onto the admin fill data.table, sorting to ensure consistent
#  cell order
ad2_fill_dt <- merge(
  x = ad2_fill_dt,
  y = pov_by_ad2,
  by = 'ADM2_CODE',
  all.x = TRUE,
  suffixes = c('_grid_cell','_adm2')
)[order(sort_order)]
ad2_fill_dt[, poverty_filled := pov_grid_cell ]
ad2_fill_dt[is.na(poverty_filled), poverty_filled := pov_adm2 ]

# Insert filled data back into raster
pov_resampled_filled <- insertRaster(
  r_template_full, new_vals=matrix(ad2_fill_dt$poverty_filled, ncol=1)
)
# Save that goodness!
writeRaster(
  pov_resampled_filled, format="GTiff", overwrite=TRUE,
  file=sprintf('%s/wp_poverty_pak_resampled_filled.tif',asset_dir)
)
# Plot final raster for vetting purposes
pdf(sprintf('%s/wp_poverty_pak_resampled_filled.pdf',asset_dir), height=10, width=10)
raster::plot(
  pov_resampled_filled, col=viridis::viridis(n=100), main='Resampled and Filled'
)
lines(subset_shape)
dev.off()

