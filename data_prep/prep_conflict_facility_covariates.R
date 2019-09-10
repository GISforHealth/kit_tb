## #############################################################################
## 
## PREPARE NEW COVARIATES
## 
## 
## 
## 
## #############################################################################


## SET INPUTS

core_repo <- '/share/code/geospatial/nathenry/lbd_core/'

asset_dir <- '/share/geospatial/mbg/tb/kit/model_assets/'
lu_dir <- '/ihme/limited_use/LIMITED_USE/LU_GEOSPATIAL/KIT_TB_Project/'
facil_fp <- paste0(lu_dir, 'facility/pak_dis_facility_dens.csv')
conflict_fp <- paste0(lu_dir, 'conflict/all_years_Adm2_conflict_violence.csv')
years <- 2010:2018

## SETUP AND LOAD DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load all MBG functions and libraries
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)

# Load population raster (which will also be a template raster)
pop_rr <- readRDS(paste0(asset_dir,'pop_15pl_pak.RDS'))

# Load admin2 shapefile and rasterize to template
ad2_shp <- sf::st_read(paste0(asset_dir,'pak_ad2.shp'))

# Load link table and corresponding simple raster
link_list <- readRDS(paste0(asset_dir,'/pak_link_list.RDS'))
subset_shape <- load_simple_polygon(
  gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
  custom_shapefile_path = '/share/geospatial/mbg/tb/kit/model_assets/pak_ad2.shp'
)[[1]]
simple_raster <- build_simple_raster_pop(subset_shape, link_table=NULL)$simple_raster

# Load raw covariate datasets
facil <- fread(facil_fp)
conflict <- fread(conflict_fp)

## PREP FACILITY DENSITY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get total population by district
if(any(dim(pop_rr)[1:2] != dim(simple_raster)[1:2])) stop("Issue with pop dimensions")
pop_dt <- data.table(
  dummy = rep(as.vector(simple_raster), length(years)),
  pop = as.vector(pop_rr),
  year = rep(years, each=prod(dim(simple_raster)[1:2]))
)[!is.na(dummy), .(year, pop)]
# Associate each pixel with a cellIdx
pop_dt$pixel_id <- rep(cellIdx(simple_raster), length(years))
# Merge on the link table and fractionally sum population
link_table_ad2 <- link_list$ad2$link_table[, .(pixel_id, ADM2_CODE, area_fraction)]
pop_dt <- merge(
  x = pop_dt,
  y = link_table_ad2,
  by = 'pixel_id',
  allow.cartesian = TRUE
)
pop_dt[, pop_adj := pop*area_fraction ]
pop_agg <- pop_dt[, .(pop = sum(pop_adj)), by=.(ADM2_CODE, year)]

## Get number of TB facilities per 100,000 people by district
facil <- facil[, .(num_facils = sum(Facility_density)), by=ADM2_CODE]
facil_pop <- merge(
  x = facil,
  y = pop_agg[year==2018, .(pop, ADM2_CODE)],
  by = 'ADM2_CODE'
)
facil_pop[, facils_per_100k := num_facils / pop * 1E5 ]

## Merge district data, number of facilities, and population in 2018
facil_sf <- merge(
  x = ad2_shp,
  y = facil_pop,
  by = 'ADM2_CODE',
  all.x = TRUE
)
facil_sf[is.na(facil_sf$facils_per_100k), 'facils_per_100k'] <- 0
# NA out India-controlled Pakistan (no data)
facil_sf[facil_sf$ADM2_NAME=='India-administered Kashmir', 'facils_per_100k'] <- NA

## Rasterize and save to file
facil_r <- fasterize::fasterize(
  sf=facil_sf,
  raster=simple_raster,
  field = 'facils_per_100k'
)
writeRaster(
  facil_r, format='GTiff', overwrite=TRUE,
  file=paste0(asset_dir,'/facils_per_100k_covar.tif')
)
pdf(paste0(asset_dir,'/facils_per_100k_covar.pdf'), height=10, width=10)
plot(facil_r)
dev.off()

## PREP CONFLICT COVARS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Merge binary conflict data with shapefile and rasterize by year
create_conflict_cov <- function(yr, field){
  conflict_sub <- copy(conflict[year==yr,])
  setnames(conflict_sub, field, 'to_map')
  ad2_data <- merge(
    x = ad2_shp,
    y = conflict_sub,
    by.x = 'ADM2_CODE',
    by.y = 'adm_uid',
    all.x = TRUE
  )
  ad2_data[is.na(ad2_data$to_map), 'to_map'] <- 0
  # NA out India-controlled Pakistan
  ad2_data[ad2_data$ADM2_NAME=='India-administered Kashmir', 'to_map'] <- NA
  # Rasterize
  ad2_r <- fasterize::fasterize(
    sf = ad2_data,
    raster = simple_raster,
    field = 'to_map'
  )
  return(ad2_r)
}

# Apply rasterizing function to protest field
protest_rr <- raster::brick(
  lapply(2010:2018, function(yr) create_conflict_cov(yr=yr, field='protest'))
)
names(protest_rr) <- paste0('protest.',2010:2018)

# Apply rasterizing function to violence field
violence_rr <- raster::brick(
  lapply(2010:2018, function(yr) create_conflict_cov(yr=yr, field='violence'))
)
names(violence_rr) <- paste0('violence.',2010:2018)

## Save rasters
writeRaster(
  protest_rr, format="GTiff", overwrite=TRUE,
  file=paste0(asset_dir, '/acled_protest_covar.tif')
)
pdf(paste0(asset_dir,'/acled_protest_covar.pdf'), height=15, width=15)
plot(protest_rr)
dev.off()

writeRaster(
  violence_rr, format="GTiff", overwrite=TRUE,
  file=paste0(asset_dir, '/acled_violence_covar.tif')
)
pdf(paste0(asset_dir,'/acled_violence_covar.pdf'), height=15, width=15)
plot(violence_rr)
dev.off()
