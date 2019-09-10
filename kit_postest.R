## #############################################################################
## 
## Postestimation script for KIT TB
## Author: Nat Henry
## Created: August 16, 2019
## Purpose: Raking, aggregation, and visualization script for KIT TB
## 
## NOTE: This code pulls from many functions in the LBDCore package, which is
##   still in development. Key functions from this package will be added to
##   this repository as a demonstration. This package is on track to be
##   publicly released as an open-source tool.
##
## #############################################################################

## NOTE: Most globals are passed in from launch script

# Custom worldpop age aggregation
wp_age_groups <- c(
  'a1549t','a5054m','a5054f','a5559m','a5559f','a6064m','a6064f','a65plf','a65plm'
)
pop_release <- '2017_04_27'
years <- 2010:2018


## SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message(sprintf("**** RUNNING POSTESTIMATION FOR %s ****",rd))
## Set input filepaths
root     <- "/home/j/"
sharedir <- sprintf('<<FILEPATH>>')
rd_dir   <- sprintf('%s/output/%s/',sharedir,rd)

## Load libraries and MBG project functions.
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)

## Load GBD libraries with relevant information
source('<<FILEPATH>>')


## LOAD DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("** Loading data **")
# Get admin0 code for region
ad0_code <- get_adm0_codes(reg)
gbd_loc_ids <- load_adm0_lookup_table()[gadm_geoid %in% ad0_code, loc_id]
# GBD 2019 subnationals
gbd_subnats <- get_location_metadata(
  gbd_round_id=6, location_set_id=35
)[parent_id %in% gbd_loc_ids, location_id]

## Load simple raster
simple_raster_fp <- '<<FILEPATH>>'
attach(simple_raster_fp) # Create a temporary environment from this Rdata object
simple_raster <- simple_raster # Add an object from the temporary env to your working env
detach()

## Load UNRAKED cell pred
unraked_cp <- get(load('<<FILEPATH>>'))

## Load UNRAKED mean raster
unraked_mean_ras <- raster::brick('<<FILEPATH>>')

## Load population, pulling from cache if possible
pop_cache_fp <- sprintf('%s/model_assets/pop_15pl_%s.RDS',sharedir,reg)
if(file.exists(pop_cache_fp)){
  pop_15pl_rr <- readRDS(pop_cache_fp)
} else {
  message("Constructing population and cacheing...")  
  pop_15pl_rr <- Reduce("+", lapply(wp_age_groups, function(ag){
    message(sprintf("Loading %s...",ag))
    covariate_config <- data.table(
      covariate='worldpop', measure=ag, release=pop_release
    )
    loader <- MbgStandardCovariateLoader$new(start_year = min(years),
                                             end_year = max(years),
                                             interval = 12,
                                             covariate_config = covariate_config)
    return(loader$get_covariates(simple_raster)[[1]])
  }))
  saveRDS(pop_15pl_rr, file=pop_cache_fp)
}

## Load admin0, admin1, admin2 for Pakistan
ad0_sf <- sf::st_read(sprintf('%s/model_assets/pak_ad0.shp',sharedir))
ad1_sf <- sf::st_read(sprintf('%s/model_assets/pak_ad1.shp',sharedir))
ad2_sf <- sf::st_read(sprintf('%s/model_assets/pak_ad2.shp',sharedir))

## Load link table from cache if possible; otherwise build
link_list_fp <- sprintf('%s/model_assets/%s_link_list.RDS',sharedir,reg)
if(file.exists(link_list_fp)){
  link_list_full <- readRDS(link_list_fp)
} else {
  # Create the link list from scratch and cache
  link_list_full <- vector('list', length=3)
  for(lev in 0:2){
    link_list_full[[sprintf('ad%s',lev)]] <- suppressMessages(build_link_table(
      cores=get_max_forked_threads(nobjs=10),
      custom_shapefile_path=sprintf('%s/model_assets/pak_ad%s.shp',sharedir,lev),
      custom_shapefile_field=sprintf('ADM%s_CODE',lev)
    ))
  }
  # Cache
  saveRDS(link_list_full, file=link_list_fp)
}

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
  sex_id = 3,
  year_id = years[years < 2018],
  measure_id = 5, # Prevalence
  metric_id = 3, # Rate
  gbd_round_id = 5
)[, .(val=sum(val, na.rm=T)), by=.(age_group_id, sex_id, year_id, location_id)]

# Pull population by age and sex
gbd2017_pops <- get_population(
  age_group_id = over_15_ags,
  location_id = gbd_loc_ids,
  sex_id = 3,
  year_id = years[years < 2018],
  gbd_round_id = 5
)[, .(age_group_id, sex_id, year_id, location_id, population)]
# Run population-weighted aggregation
gbd2017_ad0 <- merge(
  x=gbd2017_ad0_tb_agesex_specific, y=gbd2017_pops,
  by=c('age_group_id','sex_id','location_id','year_id')
)[,
  .(gbd_est = weighted.mean(val, w=population)), by=.(location_id,year_id)
]
# Impute a GBD 2018 value based on annualized rate of change for observed years
gbd_2018 <- merge(
  x = gbd2017_ad0[year_id==years[1], .(gbd_est, location_id)],
  y = gbd2017_ad0[year_id==rev(years)[2], .(gbd_est, location_id)],
  by = 'location_id',
  suffixes = c("_fy", "_ly")
)
gbd_2018[, aroc := (gbd_est_ly/gbd_est_fy)**(1/diff(range(years)))]
gbd_2018[, gbd_est := gbd_est_ly * aroc ]
gbd_2018[, year_id := 2018 ]
gbd_2018[, c('gbd_est_fy','gbd_est_ly','aroc') := NULL ]
# Merge back onto GBD 2017 estimates
gbd2017_ad0 <- rbindlist(list(gbd2017_ad0, gbd_2018), fill=FALSE, use.names=TRUE)


## CREATE SUMMARY RASTERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("** Creating summary rasters **")
## Summarize
unraked_mean_rr <- insertRaster(
  raster=simple_raster,
  new_vals=matrix(rowMeans(unraked_cp, na.rm=TRUE), ncol=length(years))
)
unraked_lower_rr <- insertRaster(
  raster=simple_raster,
  new_vals=matrix(
    apply(unraked_cp,1,function(x) quantile(x, 0.025, na.rm=TRUE)), 
    ncol=length(years)
  )
)
unraked_median_rr <- insertRaster(
  raster=simple_raster,
  new_vals=matrix(
    apply(unraked_cp,1,function(x) quantile(x, 0.5, na.rm=TRUE)), 
    ncol=length(years)
  )
)
unraked_upper_rr <- insertRaster(
  raster=simple_raster,
  new_vals=matrix(
    apply(unraked_cp,1,function(x) quantile(x, 0.975, na.rm=TRUE)), 
    ncol=length(years)
  )
)
## Save
writeRaster(
  unraked_mean_rr,
  file=sprintf('%s/%s_mean_unraked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)
writeRaster(
  unraked_lower_rr,
  file=sprintf('%s/%s_lower_unraked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)
writeRaster(
  unraked_median_rr,
  file=sprintf('%s/%s_median_unraked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)
writeRaster(
  unraked_upper_rr,
  file=sprintf('%s/%s_upper_unraked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)


## CUSTOM FRACTIONAL AGGREGATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("** Fractionally aggregating to a custom shapefile **")
## Convert and prep cell_pred object
# Convert to data.table
cp_dt <- as.data.table(unraked_cp)
# Check number of draws
ndraw <- ncol(cp_dt)
draw_cols <- paste0('V',1:ndraw)
names(cp_dt) <- draw_cols
# Extract population at simple_raster and add to cell pred
if(any(dim(pop_15pl_rr)[1:2] != dim(simple_raster)[1:2])) stop("Issue with pop dimensions")
pop_dt <- data.table(
  dummy = rep(as.vector(simple_raster), length(years)),
  pop = as.vector(pop_15pl_rr),
  year = rep(years, each=prod(dim(simple_raster)[1:2]))
)[!is.na(dummy), .(year, pop)]
if(nrow(pop_dt) != nrow(cp_dt)) stop("Issue with pop dimensions")
cp_dt <- cbind(pop_dt, cp_dt)
# Associate each pixel with a cellIdx
cp_dt$pixel_id <- rep(cellIdx(simple_raster), length(years))

## Merge on link table and prep
link_table_ad2 <- link_list_full$ad2$link_table[,
  .(pixel_id, ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME, 
    area_fraction)
]
cp_dt <- merge(
  x = cp_dt,
  y = link_table_ad2,
  by = 'pixel_id',
  allow.cartesian=TRUE
)
cp_dt[, pop_adj := pop * area_fraction ]

## Aggregate by admin unit
# Admin2
ad2_agg <- cp_dt[, 
  lapply(.SD, function(x) weighted.mean(x, w=pop_adj,na.rm=T)), 
  by=.(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME, year),
  .SDcols=draw_cols
]
ad2_agg$pop <- cp_dt[, 
  .(pop=sum(pop_adj)), 
  by=.(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME, year)
]$pop
ad2_agg <- ad2_agg[ADM2_CODE != 0,]
# Admin1
ad1_agg <- cp_dt[, 
  lapply(.SD, function(x) weighted.mean(x, w=pop_adj,na.rm=T)), 
  by=.(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, year),
  .SDcols=draw_cols
]
ad1_agg$pop <- cp_dt[, 
  .(pop=sum(pop_adj)), 
  by=.(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, year)
]$pop
ad1_agg <- ad1_agg[ADM1_CODE != 0,]
# Admin0
ad0_agg <- cp_dt[, 
  lapply(.SD, function(x) weighted.mean(x, w=pop_adj,na.rm=T)), 
  by=.(ADM0_CODE, ADM0_NAME, year),
  .SDcols=draw_cols
]
ad0_agg$pop <- cp_dt[,
  .(pop=sum(pop_adj)), by=.(ADM0_CODE, ADM0_NAME, year)
]$pop


## Summarize by admin unit
summarize_admin_draws <- function(adm_draws, num_draws){
  draw_cols <- paste0('V',1:num_draws)
  draw_dt <- adm_draws[, ..draw_cols]
  id_cols <- names(adm_draws)[!(names(adm_draws) %in% draw_cols)]
  id_dt <- adm_draws[, ..id_cols]
  id_dt$mean <- rowMeans(draw_dt, na.rm=TRUE)
  id_dt$lower <- apply(draw_dt, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
  id_dt$median <- apply(draw_dt, 1, function(x) quantile(x, 0.5, na.rm=TRUE))
  id_dt$upper <- apply(draw_dt, 1, function(x) quantile(x, 0.975, na.rm=TRUE))
  return(id_dt)
}
ad2_summ <- summarize_admin_draws(ad2_agg, ndraw)
ad1_summ <- summarize_admin_draws(ad1_agg, ndraw)
ad0_summ <- summarize_admin_draws(ad0_agg, ndraw)

## Save outputs to shared directory
wcsv <- function(dt, ff) write.csv(dt, file=ff, row.names=FALSE)
wcsv(ad2_agg, sprintf('%s/%s_unraked_ad2_draws.csv',rd_dir,indicator))
wcsv(ad1_agg, sprintf('%s/%s_unraked_ad1_draws.csv',rd_dir,indicator))
wcsv(ad0_agg, sprintf('%s/%s_unraked_ad0_draws.csv',rd_dir,indicator))
wcsv(ad2_summ, sprintf('%s/%s_unraked_ad2_fullsummary.csv',rd_dir,indicator))
wcsv(ad1_summ, sprintf('%s/%s_unraked_ad1_fullsummary.csv',rd_dir,indicator))
wcsv(ad0_summ, sprintf('%s/%s_unraked_ad_fullsummary.csv',rd_dir,indicator))


## RAKE INDICATOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("** Raking indicators **")
## Determine raking factors by year
rf_dt <- merge(
  x=ad0_summ[, .(year, mean)],
  y=gbd2017_ad0[, .(year_id, gbd_est)],
  by.x='year', by.y='year_id'
)
setnames(rf_dt, 'mean', 'lbd_est')
rf_dt[, rf := gbd_est / lbd_est ]

write.csv(rf_dt, row.names=FALSE, file=sprintf('%s/%s_rfs_%s.csv',rd_dir,indicator,reg))

## Apply raking factors to cell pred and save
raker <- do.call('rbind', lapply(years, function(yy){
  matrix(rf_dt[year==yy, rf], ncol=ndraw, nrow=length(cellIdx(simple_raster)))
}))
raked_cp <- unraked_cp * raker
save(raked_cp, file=sprintf(
  '%s/output/%s/%s_cell_draws_raked_eb_bin0_%s_0.RData',sharedir,rd,indicator,reg
))

## Apply raking factors to aggregates and save
apply_rfs <- function(lbd_est, rfs, cols){
  lbd_est <- merge(x=lbd_est, y=rfs[, .(year, rf)], by='year')
  lbd_est <- lbd_est[, (cols) := lapply(.SD, '*', rf), .SDcols=cols]
  return(lbd_est)
}
ad2_agg_raked <- apply_rfs(ad2_agg, rf_dt, cols=draw_cols)
ad1_agg_raked <- apply_rfs(ad1_agg, rf_dt, cols=draw_cols)
ad0_agg_raked <- apply_rfs(ad0_agg, rf_dt, cols=draw_cols)
ad2_summ_raked <- apply_rfs(ad2_summ, rf_dt, cols=c('mean','lower','median','upper'))
ad1_summ_raked <- apply_rfs(ad1_summ, rf_dt, cols=c('mean','lower','median','upper'))
ad0_summ_raked <- apply_rfs(ad0_summ, rf_dt, cols=c('mean','lower','median','upper'))
# Save
wcsv(ad2_agg_raked, sprintf('%s/%s_raked_ad2_draws.csv',rd_dir,indicator))
wcsv(ad1_agg_raked, sprintf('%s/%s_raked_ad1_draws.csv',rd_dir,indicator))
wcsv(ad0_agg_raked, sprintf('%s/%s_raked_ad0_draws.csv',rd_dir,indicator))
wcsv(ad2_summ_raked, sprintf('%s/%s_raked_ad2_fullsummary.csv',rd_dir,indicator))
wcsv(ad1_summ_raked, sprintf('%s/%s_raked_ad1_fullsummary.csv',rd_dir,indicator))
wcsv(ad0_summ_raked, sprintf('%s/%s_raked_ad0_fullsummary.csv',rd_dir,indicator))

## Apply raking factors to summary rasters and save
raked_mean_rr <- brick(lapply(1:length(years), function(i){
  unraked_mean_rr[[i]] * rf_dt[year==years[i],rf]
}))
raked_lower_rr <- brick(lapply(1:length(years), function(i){
  unraked_lower_rr[[i]] * rf_dt[year==years[i],rf]
}))
raked_median_rr <- brick(lapply(1:length(years), function(i){
  unraked_median_rr[[i]] * rf_dt[year==years[i],rf]
}))
raked_upper_rr <- brick(lapply(1:length(years), function(i){
  unraked_upper_rr[[i]] * rf_dt[year==years[i],rf]
}))
writeRaster(
  raked_mean_rr,
  file=sprintf('%s/%s_mean_raked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)
writeRaster(
  raked_lower_rr,
  file=sprintf('%s/%s_lower_raked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)
writeRaster(
  raked_median_rr,
  file=sprintf('%s/%s_median_raked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)
writeRaster(
  raked_upper_rr,
  file=sprintf('%s/%s_upper_raked_%s_%s',rd_dir,indicator,min(years),max(years)),
  format = "GTiff", overwrite=TRUE
)


## VISUALIZATIONS AND VETTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("** Running visualizations **")
## Make maps
viz_dir <- sprintf('%s/visualizations/',rd_dir)
dir.create(viz_dir, showWarnings=FALSE)
col_scale <- rev(RColorBrewer::brewer.pal("Spectral", n=11))
col_breaks <- c(0, 100, 200, 400, 600)
col_max <- max(col_breaks)

make_tb_map <- function(sf_obj, data_dt, merge_col, map_col, yr, title){
  # Prep data
  data_for_merge <- data_dt[year==yr, ]
  keep_cols <- c(merge_col, map_col)
  data_for_merge <- data_for_merge[, ..keep_cols]
  setnames(data_for_merge, map_col, 'map_col')
  # Define color scale
  data_for_merge[, map_col := map_col * 1E5]
  data_for_merge[map_col > col_max, map_col := col_max ]
  # Merge on data
  sf_for_mapping <- merge(x=sf_obj, y=data_for_merge, by=merge_col, all.x=TRUE)
  sf_for_mapping$map_col <- as.numeric(sf_for_mapping$map_col)
  # Make plot
  fig <- ggplot() + 
    geom_sf(data=sf_for_mapping, aes(fill=map_col), size=.5, color='#444444') + 
    scale_fill_gradientn(
      colors=col_scale, limits=range(col_breaks), breaks=col_breaks
    ) + 
    labs(
      title=title,
      fill='TB Prevalence\nper 100k'
    ) + 
    theme_minimal()
  return(fig)
}

## PLOT MEANS
for(lev in 0:2){
  pdf(sprintf('%s/adm%i_maps.pdf',viz_dir,lev), height=10, width=10)
  in_sf <- get(sprintf('ad%i_sf',lev))
  in_data <- get(sprintf('ad%i_summ_raked',lev))
  for(year in years){
    fig <- make_tb_map(
      sf_obj=in_sf, data_dt=in_data, merge_col=sprintf('ADM%i_CODE',lev),
      map_col='mean', yr=year,
      title=sprintf('TB prevalence at the Admin %i level in %i',lev,year)
    )
    print(fig)
  }
  dev.off()
}
## PLOT MEDIANS
for(lev in 0:2){
  pdf(sprintf('%s/adm%i_median_maps.pdf',viz_dir,lev), height=10, width=10)
  in_sf <- get(sprintf('ad%i_sf',lev))
  in_data <- get(sprintf('ad%i_summ_raked',lev))
  for(year in years){
    fig <- make_tb_map(
      sf_obj=in_sf, data_dt=in_data, merge_col=sprintf('ADM%i_CODE',lev),
      map_col='median', yr=year,
      title=sprintf('TB prevalence at the Admin %i level in %i',lev,year)
    )
    print(fig)
  }
  dev.off()
}


## Make admin2 trends plot WITH DATA
prepped_data <- readRDS('<<FILEPATH>>')
ad2_sp <- as(ad2_sf, 'Spatial')
prepped_meta <- as.data.table(sp::over(
  SpatialPoints(
    cbind(prepped_data$longitude, prepped_data$latitude),
    proj4string=CRS(proj4string(ad2_sp))
  ),
  ad2_sp
))[, .(ADM2_CODE, ADM1_CODE)]

prepped_data <- cbind(prepped_data, prepped_meta)
prepped_data[, year := 2011 ]
prepped_data[, resampled := as.character(as.integer(pseudocluster))]
prepped_data[resampled == '0', resampled := "No"]
prepped_data[resampled == '1', resampled := "Yes"]
prepped_data <- prepped_data[!is.na(resampled)]

## Merge data and admin time trends
ad1_summ_plotting <- merge(
  x = ad1_summ,
  y = prepped_data[, .(ADM1_CODE, year, kit, N, weight, resampled)],
  by = c('ADM1_CODE','year'),
  all.x = TRUE,
  sort= FALSE
)[order(ADM1_CODE, year)]
ad1_summ_plotting[, mean_plot := mean * 1E5 ]
ad1_summ_plotting[, lower_plot := lower * 1E5 ]
ad1_summ_plotting[, upper_plot := upper * 1E5 ]
ad1_summ_plotting[, data_rate := kit/N * 1E5 ]
ad1_summ_plotting[, N_adj := N * weight ]

ad2_summ_plotting <- merge(
  x = ad2_summ,
  y = prepped_data[, .(ADM2_CODE, year, kit, N, weight, resampled)],
  by = c('ADM2_CODE','year'),
  all.x = TRUE,
  sort= FALSE
)
ad2_summ[, mean_plot := mean * 1E5 ]
ad2_summ_plotting[, mean_plot := mean * 1E5 ]
ad2_summ_plotting[, lower_plot := lower * 1E5 ]
ad2_summ_plotting[, upper_plot := upper * 1E5 ]
ad2_summ_plotting[, data_rate := kit/N * 1E5 ]
ad2_summ_plotting[, N_adj := N * weight ]

admin_trends_plot_ad1 <- ggplot(ad1_summ_plotting) + 
  facet_wrap("ADM1_NAME", ncol=4) +
  geom_ribbon(
    aes(x=year, ymin=lower_plot, ymax=upper_plot),
    fill = '#01CFD8', # LBD in blue
    color = '#222222',
    alpha = 0.5,
    lwd = .5
  ) + 
  geom_line(
    aes(x=year, y=mean_plot),
    color = '#222222',
    lwd=1
  ) +
  geom_line(
    data=ad2_summ,
    aes(x=year, y=mean_plot, group=ADM2_CODE),
    color = '#D53E4F',
    lwd = .25,
    alpha = .5
  ) +
  geom_point(
    data=ad1_summ_plotting,
    aes(x=year, y=data_rate, size=N_adj, color=resampled),
    alpha=.5
  ) + 
  theme_bw() + 
  lims(x=range(years)) + 
  coord_cartesian(ylim=c(0, 1500)) +
  scale_color_manual(values=c('#ff33cc','#669900')) + 
  labs(
    title='Province and District-level Time Trends',
    x="Year",
    y="TB Prevalence per 100k",
    size='Adjusted\nSample\nSize',
    color="Resampled?"
  )

pdf(sprintf('%s/admin1_time_series.pdf',viz_dir), height=7, width=12)
print(admin_trends_plot_ad1)
dev.off()



## NEW PLOT BY ADMIN2 UNIT
admin_trends_plot_ad2 <- ggplot(ad2_summ_plotting) + 
  facet_wrap("ADM2_NAME", ncol=6) +
  geom_ribbon(
    aes(x=year, ymin=lower_plot, ymax=upper_plot),
    fill = '#01CFD8', # LBD in blue
    color = '#222222',
    alpha = 0.5,
    lwd = .5
  ) + 
  geom_line(
    aes(x=year, y=mean_plot),
    color = '#222222',
    lwd=1
  ) +
  geom_point(
    data=ad2_summ_plotting,
    aes(x=year, y=data_rate, size=N_adj, color=resampled),
    alpha=.5
  ) + 
  theme_bw() + 
  lims(x=range(years)) + 
  coord_cartesian(ylim=c(0, 1500)) +
  scale_color_manual(values=c('#ff33cc','#669900')) + 
  labs(
    title='District-level Time Trends',
    x="Year",
    y="TB Prevalence per 100k",
    size='Adjusted\nSample\nSize',
    color="Resampled?"
  )

pdf(sprintf('%s/admin2_time_series.pdf',viz_dir), height=30, width=18)
print(admin_trends_plot_ad2)
dev.off()


pdf(sprintf('%s/admin_time_series_by_ad2.pdf',viz_dir), height=7, width=12)
for(ad2_name in unique(ad2_summ_plotting$ADM2_NAME)){

  fig <- ggplot(data=ad2_summ_plotting[ADM2_NAME == ad2_name,], aes(x=year)) + 
    geom_ribbon(
      aes(x=year, ymin=lower_plot, ymax=upper_plot),
      fill = '#01CFD8', # LBD in blue
      color = '#222222',
      alpha = 0.5,
      lwd = .5
    ) + 
    geom_line(
      aes(x=year, y=mean_plot),
      color = '#222222',
      lwd=1
    ) + 
    geom_point(
      aes(y=data_rate, size=N),
      color = 'blue'
    ) + 
    labs(title=ad2_name)
  print(fig)
}
dev.off()


## NEW PLOT BY PIXEL
data_spatialpts <- SpatialPoints(
  cbind(prepped_data$longitude, prepped_data$latitude),
  proj4string=CRS(proj4string(ad2_sp))
)
prepped_data[, rs_cluster_i := .I ]
px_level_list <- vector('list', length=length(years))
for(i in 1:length(years)){
  yr <- years[i]
  px_level_list[[i]] <- data.table(
    rs_cluster_i = prepped_data$rs_cluster_i,
    px_mean = raster::extract(x=raked_mean_rr[[i]], y=data_spatialpts),
    px_median = raster::extract(x=raked_median_rr[[i]], y=data_spatialpts),
    px_lower = raster::extract(x=raked_lower_rr[[i]], y=data_spatialpts),
    px_upper = raster::extract(x=raked_upper_rr[[i]], y=data_spatialpts),
    year = yr
  )    
  if(yr==2011){
    px_level_list[[i]]$cluster_number <- prepped_data$cluster_number
    px_level_list[[i]]$N <- prepped_data$N
    px_level_list[[i]]$kit <- prepped_data$kit
    px_level_list[[i]]$weight <- prepped_data$weight
    px_level_list[[i]]$resampled <- as.integer(prepped_data$pseudocluster)
  }
}
px_data_full <- rbindlist(px_level_list, use.names=TRUE, fill=TRUE)
px_data_full[, outcome := kit / N ]
px_data_full[, N_adj := N * weight]

# Multiply everything by 100k
px_data_full[, px_mean := px_mean * 1E5 ]
px_data_full[, px_lower := px_lower * 1E5 ]
px_data_full[, px_upper := px_upper * 1E5 ]
px_data_full[, outcome := outcome * 1E5 ]
px_data_full[, resampled := as.character(resampled)]
px_data_full[resampled == '0', resampled := "No"]
px_data_full[resampled == '1', resampled := "Yes"]

all_clusters <- sort(unique(prepped_data$rs_cluster_i))

pdf(sprintf('%s/admin_time_series_by_pixel_onepage.pdf',viz_dir), height=60, width=24)
fig <- ggplot(data=px_data_full, aes(x=year)) + 
  facet_wrap(~rs_cluster_i, ncol=8) + 
  geom_ribbon(
    aes(x=year, ymin=px_lower, ymax=px_upper),
    fill = '#01CFD8', # LBD in blue
    color = '#222222',
    alpha=.5,
    lwd=.5
  ) +
  geom_line(
    aes(x=year, y=px_mean), color='#222222', lwd=1
  ) +
  geom_point( 
    aes(x=year, y=outcome, size=N_adj, color=resampled)
  ) + 
  labs(
    title="Time trend plots: Pixels overlapping data",
    x='Year',
    y='TB Prevalence per 100k',
    size='Weighted\nSample\nSize',
    color='Resampled?'
  ) + 
  scale_color_manual(values=c('#ff33cc','#669900')) + 
  coord_cartesian(ylim=c(0, 1500)) +
  theme_bw()
print(fig)
dev.off()

px_11 <- px_data_full[year==2011,]
cluster_agg <- unique(px_11[, .(cluster_number, N, outcome)])
px_11_agg <- px_11[, .(
    px_mean = weighted.mean(px_mean, w=weight), 
    px_median = weighted.mean(px_median, w=weight) * 1E5
  ),
  by = .(cluster_number)
]

px_11_agg <- merge(x=px_11_agg, y=cluster_agg, by='cluster_number')

test_outcome <- na.omit(px_11_agg[, .(px_mean, outcome, N)])
outcome_r2 <- boot::corr(
  d=cbind(test_outcome$px_mean,test_outcome$outcome),
  w=test_outcome$N
)**2
test_outcome_median <- na.omit(px_11_agg[, .(px_median, outcome, N)])
outcome_r2_median <- boot::corr(
  d=cbind(test_outcome_median$px_median,test_outcome_median$outcome),
  w=test_outcome$N
)**2



## Test correlation between pixels and data in 2011
correlation_fig <- ggplot(data=px_11_agg) + 
  geom_abline(intercept=0, slope=1, linetype=2, color='#666666') +
  geom_point(aes(x=outcome, y=px_mean, size=N), color='blue', alpha=.5) +
  theme_bw() + 
  labs(
    title='Correlation between data and modeled estimates',
    subtitle=sprintf('Observed in 2011. R squared: %s',round(outcome_r2,4)),
    x='Raw data estimate (per 100k)',
    y='Modeled value at pixel (mean, per 100k)',
    size="Sample Size"
  )
pdf(sprintf('%s/data_estimate_correlation_px.pdf',viz_dir), height=7, width=7)
print(correlation_fig)
dev.off()

correlation_fig <- ggplot(data=px_11_agg) + 
  geom_abline(intercept=0, slope=1, linetype=2, color='#666666') +
  geom_point(aes(x=outcome, y=px_median, size=N), color='blue', alpha=.5) +
  theme_bw() + 
  labs(
    title='Correlation between data and modeled estimates',
    subtitle=sprintf('Observed in 2011. R squared: %s',round(outcome_r2_median,4)),
    x='Raw data estimate (per 100k)',
    y='Modeled value at pixel (median, per 100k)',
    size="Sample Size"
  )
pdf(sprintf('%s/data_estimate_correlation_px_median.pdf',viz_dir), height=7, width=7)
print(correlation_fig)
dev.off()


## END OF POSTESTIMATION SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("\n**** All done with postestimation!! ****\n")
