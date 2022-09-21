## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  eval = FALSE
)

## -----------------------------------------------------------------------------
#  install.packages('pRecipe')
#  library(pRecipe)

## -----------------------------------------------------------------------------
#  download_data(name = "era5", project_folder = ".")
#  show_info("./data/database/era5_tp_mm_global_195901_202112_025_monthly.nc")

## -----------------------------------------------------------------------------
#  subset_spacetime("era5", start_year = 1981, end_year = 2020,
#                   bbox = c(2,28,42,58),
#                   database_path = "./data/database/")
#  show_info("./data/processed/era5_tp_mm_subset_1981_2020_025_monthly.nc")

## -----------------------------------------------------------------------------
#  crop_data(nc_path = "./data/processed/era5_tp_mm_subset_1981_2020_025_monthly.nc",
#            shp_path = "./shapefiles/CZE_adm0.shp",
#            save_nc = TRUE)
#  show_info("./data/processed/era5_tp_mm_subset_1981_2020_025_monthly_cropped.nc")

## -----------------------------------------------------------------------------
#  make_ts("era5", database_path = "./data/database/")
#  tp_ts <- readr::read_csv("./data/processed/era5_ts.csv")
#  head(tp_ts, 12)

## -----------------------------------------------------------------------------
#  era5_subset <- raster::brick("./data/processed/era5_tp.mm_subset_1981_2020_025_monthly.nc")
#  make_ts(era5_subset)
#  subset_ts <- readr::read_csv("./data/processed/era5_tp.mm_subset_1981_2020_025_monthly_ts.csv")
#  head(subset_ts, 12)

## -----------------------------------------------------------------------------
#  era5_cropped <- raster::brick("./data/processed/era5_tp.mm_subset_1981_2020_025_monthly_cropped.nc")
#  make_ts(era5_cropped)
#  cropped_ts <- readr::read_csv("./data/processed/era5_tp.mm_subset_1981_2020_025_monthly_cropped_ts.csv")
#  head(cropped_ts, 12)

## -----------------------------------------------------------------------------
#  global <- raster::raster("./data/database/era5_tp_mm_global_195901_202112_025_monthly.nc")
#  plot_map(global[[1]])

## -----------------------------------------------------------------------------
#  central_europe <- raster::brick("./data/processed/era5_tp_mm_subset_1981_2020_025_monthly.nc")
#  plot_map(central_europe[[1]])

## -----------------------------------------------------------------------------
#  czechia <- raster::brick("./data/processed/era5_tp_mm_subset_1981_2020_025_monthly_cropped.nc")
#  plot_map(czechia[[1]])

## -----------------------------------------------------------------------------
#  plot_line(tp_ts)
#  #plot_line(subset_ts)
#  #plot_line(cropped_ts)

## -----------------------------------------------------------------------------
#  plot_heatmap(tp_ts)
#  #plot_heatmap(subset_ts)
#  #plot_heatmap(cropped_ts)

## -----------------------------------------------------------------------------
#  plot_box(tp_ts)
#  #plot_box(subset_ts)
#  #plot_box(cropped_ts)

## -----------------------------------------------------------------------------
#  plot_density(tp_ts)
#  #plot_density(subset_ts)
#  #plot_density(cropped_ts)

## -----------------------------------------------------------------------------
#  plot_summary(tp_ts)
#  #plot_summary(subset_ts)
#  #plot_summary(cropped_ts)

