## ----start, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  eval = TRUE,
  fig.width = 7,
  warning = FALSE,
  message = FALSE
)
library(pRecipe)
library(kableExtra)

## ----gauge, echo=FALSE, results = 'asis'--------------------------------------
tibble::tribble(
  ~"Dataset", ~"Spatial Coverage", ~"Highest Temporal Resolution Available", ~"Record Length", ~"Get Data", ~"Reference",
"CPC-Global", "Land", "Daily", "1979/01-2023/09", "[Download](https://psl.noaa.gov/data/gridded/data.cpc.globalprecip.html)", "@xie_cpc_2010",
"CRU TS v4.08", "Land", "Monthly", "1901/01-2023/12", "[Download](https://crudata.uea.ac.uk/cru/data/hrg/)", "@harris_version_2020",
"EM-Earth", "Land", "Daily", "1950/01-2019/12", "[Download](https://www.frdr-dfdr.ca/repo/dataset/8d30ab02-f2bd-4d05-ae43-11f4a387e5ad)", "@tang_em-earth_2022",
"GHCN v2", "Land", "Monthly", "1900/01-2015/05", "[Download](https://psl.noaa.gov/data/gridded/data.ghcngridded.html)", "@peterson_overview_1997",
"GPCC v2022", "Land", "Daily", "1891/01-2020/10", "[Download](https://psl.noaa.gov/data/gridded/data.gpcc.html)", "@schneider_gpcc_2011",
"PREC/L", "Land", "Monthly", "1948/01-2024/10", "[Download](https://psl.noaa.gov/data/gridded/data.precl.html)", "@chen_global_2002"
) |>
  kbl(align = 'lccccr') |>
  kable_styling("striped") |>
  unclass() |> cat()

## ----satellite, echo=FALSE, results = 'asis'----------------------------------
tibble::tribble(
  ~"Dataset", ~"Spatial Coverage", ~"Highest Temporal Resolution Available", ~"Record Length", ~"Get Data", ~"Reference",
"CHIRPS v2.0", "Land 50°SN", "Daily", "1981/01-2023/08", "[Download](https://www.chc.ucsb.edu/data/chirps)", "@funk_climate_2015",
"CMAP", "Global", "Monthly", "1979/01-2024/10", "[Download](https://psl.noaa.gov/data/gridded/data.cmap.html)", "@xie_global_1997",
"CMORPH-CDR", "Global 60°SN", "Daily", "1998/01-2023/04", "[Download](https://www.ncei.noaa.gov/data/cmorph-high-resolution-global-precipitation-estimates/)", "@joyce_cmorph_2004",
"GPCP v3.2", "Global", "Daily", "1979/01-2021/09", "[Download](https://psl.noaa.gov/data/gridded/data.gpcp.html)", "@adler_global_2018",
"GPM IMERGM Final v07", "Global", "Daily", "1998/01-2024/06", "[Download](https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDF_07/summary?keywords=GPM_3IMERGDF_07)", "@huffman_gpm_2019",
"GSMaP v8", "Global", "Daily", "1998/01-2023/06", "[Download](https://sharaku.eorc.jaxa.jp/GSMaP/)", "@kubota_global_2020",
"MSWEP v2.8", "Global", "Daily", "1979/01-2024/11", "[Download](https://www.gloh2o.org/mswep/)", "@beck_mswep_2019",
"PERSIANN-CDR", "Global 60°SN", "Daily", "1983/01-2023/12", "[Download](https://chrsdata.eng.uci.edu/)", "@ashouri_persiann-cdr_2015"
) |>
  kbl(align = 'lccccr') |>
  kable_styling("striped") |>
  unclass() |> cat()

## ----reanalysis, echo=FALSE, results = 'asis'---------------------------------
tibble::tribble(
  ~"Dataset", ~"Spatial Coverage", ~"Highest Temporal Resolution Available", ~"Record Length", ~"Get Data", ~"Reference",
"20CR v3", "Global", "Daily", "1836/01-2015/12", "[Download](https://psl.noaa.gov/data/gridded/data.20thC_ReanV3.html)", "@slivinski_towards_2019",
"ERA-20C", "Global", "Daily", "1900/01-2010/12", "[Download](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-20th-century)", "@poli_era-20c_2016",
"ERA5", "Global", "Monthly", "1959/01-2021/12", "[Download](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)", "@hersbach_era5_2020",
"ERA5-Land", "Land", "Monthly", "1959/01-2021/12", "[Download](https://www.ecmwf.int/en/era5-land)", "@munoz_era5_2021",
"JRA-55", "Global", "Daily", "1958/01-2023/09", "[Download](https://rda.ucar.edu/datasets/ds628.1/dataaccess/)", "@kobayashi_jra-55_2015",
"MERRA-2", "Global", "Daily", "1980/01-2024/10", "[Download](https://disc.gsfc.nasa.gov/datasets?page=1&project=MERRA-2)", "@gelaro_modern-era_2017",
"NCEP/NCAR R1", "Global", "Daily", "1948/01-2023/12", "[Download](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.html)", "@kalnay_ncepncar_1996",
"NCEP/DOE R2", "Global", "Daily", "1979/01-2023/12", "[Download](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html)", "@kanamitsu_ncepdoe_2002"
) |>
  kbl(align = 'lccccr') |>
  kable_styling("striped") |>
  unclass() |> cat()

## ----models, echo=FALSE, results = 'asis'-------------------------------------
tibble::tribble(
  ~"Dataset", ~"Spatial Coverage", ~"Highest Temporal Resolution Available", ~"Record Length", ~"Get Data", ~"Reference",
"FLDAS", "Land", "Monthly", "1982/01-2024/10", "[Download](https://ldas.gsfc.nasa.gov/fldas/fldas-data-download)", "@mcnally_land_2017",
"GLDAS CLSM v2.0", "Land", "Daily", "1948/01-2014/12", "[Download](https://ldas.gsfc.nasa.gov/gldas/gldas-get-data)", "@rodell_global_2004",
"GLDAS NOAH v2.0", "Land", "Monthly", "1948/01-2014/12", "[Download](https://ldas.gsfc.nasa.gov/gldas/gldas-get-data)", "@rodell_global_2004",
"GLDAS VIC v2.0", "Land", "Monthly", "1948/01-2014/12", "[Download](https://ldas.gsfc.nasa.gov/gldas/gldas-get-data)", "@rodell_global_2004",
"TerraClimate", "Land", "Monthly", "1958/01-2023/12", "[Download](https://www.climatologylab.org/terraclimate.html)", "@abatzoglou_terraclimate_2018"
) |>
  kbl(align = 'lccccr') |>
  kable_styling("striped") |>
  unclass() |> cat()

## ----install, eval = FALSE----------------------------------------------------
# install.packages('pRecipe')
# library(pRecipe)

## ----download, eval = FALSE---------------------------------------------------
# download_data(dataset = 'gpm-imerg')
# gpm_global <- raster::brick('gpm-imerg-v7_tp_mm_global_199801_202406_025_monthly.nc')
# infoNC(gpm_global)

## ----subset, eval = FALSE-----------------------------------------------------
# gpm_subset <- subset_data(gpm_global, box = c(-96, -30, -56, 24), yrs = c(2001, 2020))
# infoNC(gpm_subset)

## ----crop, eval = FALSE-------------------------------------------------------
# gpm_bol <- crop_data(gpm_subset, "gadm41_BOL_0.shp")
# infoNC(gpm_bol)

## ----global_ts, eval = FALSE--------------------------------------------------
# gpm_global_ts <- fldmean(gpm_global)
# head(gpm_global_ts, 12)

## ----subset_ts, eval = FALSE--------------------------------------------------
# gpm_subset_ts <- fldmean(gpm_subset)
# head(gpm_subset_ts, 12)

## ----bol_ts, eval = FALSE-----------------------------------------------------
# gpm_bol_ts <- fldmean(gpm_bol)
# head(gpm_bol_ts, 12)

## ----map_global, eval = FALSE-------------------------------------------------
# plot_map(gpm_global)

## ----map_subset, eval = FALSE-------------------------------------------------
# plot_map(gpm_subset)

## ----map_bo, eval = FALSE-----------------------------------------------------
# plot_map(gpm_bol)

## ----lines, eval = FALSE------------------------------------------------------
# plot_line(gpm_global_ts)

## ----lines_sa, eval = FALSE---------------------------------------------------
# plot_line(gpm_subset_ts)

## ----lines_bo, eval = FALSE---------------------------------------------------
# plot_line(gpm_bol_ts)

## ----hearmaps, eval = FALSE---------------------------------------------------
# plot_heatmap(gpm_global_ts)

## ----hearmaps_sa, eval = FALSE------------------------------------------------
# plot_heatmap(gpm_subset_ts)

## ----hearmaps_bo, eval = FALSE------------------------------------------------
# plot_heatmap(gpm_bol_ts)

## ----boxplots, eval = FALSE---------------------------------------------------
# plot_box(gpm_global_ts)

## ----boxplots_sa, eval = FALSE------------------------------------------------
# plot_box(gpm_subset_ts)

## ----boxplots_bo, eval = FALSE------------------------------------------------
# plot_box(gpm_bol_ts)

## ----histograms, eval = FALSE-------------------------------------------------
# plot_density(gpm_global_ts)

## ----histograms_sa, eval = FALSE----------------------------------------------
# plot_density(gpm_subset_ts)

## ----histograms_bo, eval = FALSE----------------------------------------------
# plot_density(gpm_bol_ts)

## ----summary, eval=FALSE------------------------------------------------------
# plot_summary(gpm_global_ts)
# #plot_summary(gpm_subset_ts)
# #plot_summary(gpm_cz_ts)

