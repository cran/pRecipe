## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  eval = TRUE,
  fig.width = 7,
  warning = FALSE,
  message = FALSE
)
library(pRecipe)
data("era5_global_ts")
data("era5_ce_ts")
data("era5_cze_ts")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('pRecipe')
#  library(pRecipe)

## ---- eval=FALSE--------------------------------------------------------------
#  download_data(data_name = "era5")
#  show_info("era5_tp_mm_global_195901_202112_025_monthly.nc")

## ---- eval=FALSE--------------------------------------------------------------
#  subset_spacetime("era5_tp_mm_global_195901_202112_025_monthly.nc",
#                   years = c(1981, 2020), bbox = c(2,28,42,58), autosave = TRUE)
#  show_info("era5_tp_mm_subset_198101_202012_025_monthly.nc")

## ---- eval=FALSE--------------------------------------------------------------
#  crop_data(data_file = "era5_tp_mm_subset_198101_202012_025_monthly.nc",
#            shp_path = "CZE_adm0.shp", autosave = TRUE)
#  show_info("era5_tp_mm_cropped_198101_202012_025_monthly.nc")

## ---- eval=FALSE--------------------------------------------------------------
#  era5_global_ts <- make_ts("era5_tp_mm_global_195901_202112_025_monthly.nc")
#  head(era5_global_ts, 12)

## ---- eval=FALSE--------------------------------------------------------------
#  era5_ce_ts <- make_ts("era5_tp_mm_subset_198101_202012_025_monthly.nc")
#  head(era5_ce_ts, 12)

## ---- eval=FALSE--------------------------------------------------------------
#  era5_cze_ts <-make_ts("era5_tp_mm_cropped_198101_202012_025_monthly.nc")
#  head(era5_cze_ts, 12)

## ---- eval=FALSE--------------------------------------------------------------
#  global <- raster::brick("era5_tp_mm_global_195901_202112_025_monthly.nc")
#  plot_map(global[[1]])

## ---- eval=FALSE--------------------------------------------------------------
#  central_europe <- raster::brick("era5_tp_mm_subset_198101_202012_025_monthly.nc")
#  plot_map(central_europe[[1]])

## ---- eval=FALSE--------------------------------------------------------------
#  czechia <- raster::brick("era5_tp_mm_cropped_198101_202012_025_monthly.nc")
#  plot_map(czechia[[1]])

## -----------------------------------------------------------------------------
plot_line(era5_global_ts)
plot_line(era5_ce_ts)
plot_line(era5_cze_ts)

## -----------------------------------------------------------------------------
plot_heatmap(era5_global_ts)
plot_heatmap(era5_ce_ts)
plot_heatmap(era5_cze_ts)

## -----------------------------------------------------------------------------
plot_box(era5_global_ts)
plot_box(era5_ce_ts)
plot_box(era5_cze_ts)

## -----------------------------------------------------------------------------
plot_density(era5_global_ts)
plot_density(era5_ce_ts)
plot_density(era5_cze_ts)

## ---- eval=FALSE--------------------------------------------------------------
#  plot_summary(era5_global_ts)
#  #plot_summary(era5_ce_ts)
#  #plot_summary(era5_cze_ts)

