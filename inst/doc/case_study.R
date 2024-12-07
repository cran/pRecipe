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
library(cowplot)
library(data.table)
library(foreach)
library(ggplot2)
library(ggpubr)
library(grid)
library(raster)
library(scales)

## ----eval = FALSE-------------------------------------------------------------
# library(cowplot)
# library(data.table)
# library(foreach)
# library(ggpubr)
# library(ggplot2)
# library(grid)
# library(kableExtra)
# library(pRecipe)
# library(raster)
# library(scales)

## ----eval = FALSE-------------------------------------------------------------
# download_data()
# prec_obs <- fread("chmi.csv")
# prec_obs$dataset <- "chmi"
# GAUGE_BASED <- c("cpc-global", "cru-ts-v4-08", "em-earth", "ghcn-v2",
#                  "gpcc-v2022", "precl")
# SATELLITE_BASED <- c("chirps-v2", "cmap", "cmorph-cdr", "gpcp-v3-2",
#                      "gpm-imerg-v7", "gsmap-v8", "mswep-v2-8", "persiann-cdr")
# REANALYSES <- c("20cr-v3", "era20c", "era5", "era5-land", "jra55", "merra-2",
#                 "ncep-doe", "ncep-ncar")
# FORCINGS <- c("fldas", "gldas-clsm-v2-0", "gldas-noah-v2-0", "gldas-vic-v2-0",
#               "terraclimate")

## ----eval = FALSE-------------------------------------------------------------
# prec_names <- list.files(full.names = TRUE)
# prec_data <- foreach (dataset_count = 1:length(prec_names)) %do% {
#   dummie_data <- subset_data(prec_names[dataset_count], box = c(11, 19, 46, 54),
#                              yrs = c(1981, 2020)) %>%
#     crop_data("gadm41_CZE_0.shp")
#   dummie_name <- sub(".*/([^_/]*)_.*", "\\1", prec_names[dataset_count])
#   dummie_data@file@name <- dummie_name
#   return(dummie_data)
#   }

## ----eval = FALSE-------------------------------------------------------------
# prec_ts <- foreach (idx = 1:length(prec_data), .combine = rbind) %do% {
#   dummie_data <- prec_data[[idx]]
#   dummie_name <- dummie_data@file@name
#   dummie_data <- fldmean(dummie_data)
#   dummie_data$dataset <- dummie_name
#   return(dummie_data)
# }
# prec_ts[dataset %in% GAUGE_BASED, source := "Gauge-based"
#         ][dataset %in% FORCINGS, source := "Model forcing"
#           ][dataset %in% REANALYSES, source := "Reanalysis"
#             ][dataset %in% SATELLITE_BASED, source := "Satellite-based"]

## ----eval = FALSE-------------------------------------------------------------
# prec_ts[, asum := sum(value, na.rm = TRUE), by = .(year(date), dataset)
#         ][, amean := mean(value, na.rm = TRUE), by = .(year(date), dataset)
#           ][, amax := max(value, na.rm = TRUE), by = .(year(date), dataset)
#             ][, amin := min(value, na.rm = TRUE), by = .(year(date), dataset)
#               ][, amed := median(value, na.rm = TRUE), by = .(year(date), dataset)]
# prec_ts$dataset <- as.factor(prec_ts$dataset)

## ----eval = FALSE-------------------------------------------------------------
# prec_color <- setNames(hue_pal()(length(unique(prec_ts$dataset))),
#                        levels(prec_ts$dataset))
# p01 <- ggplot(prec_ts, aes(x = date)) +
#   geom_line(aes(y = value, color = dataset)) +
#   scale_color_manual(values = prec_color, guide = 'none') +
#   labs(y = NULL, title = 'Monthly') + theme_bw()
# p02 <- ggplot(prec_ts, aes(x = date)) +
#   geom_line(aes(y = asum, color = dataset)) +
#   scale_color_manual(values = prec_color, guide = 'none') +
#   labs(y = NULL, title = 'Annual Total') + theme_bw()
# p03 <- ggplot(prec_ts, aes(x = date)) +
#   geom_line(aes(y = amin, color = dataset)) +
#   scale_color_manual(values = prec_color, guide = 'none') +
#   labs(y = NULL, title = 'Annual Min') + theme_bw()
# p04 <- ggplot(prec_ts, aes(x = date)) +
#   geom_line(aes(y = amax, color = dataset)) +
#   scale_color_manual(values = prec_color, guide = 'none') +
#   labs(y = NULL, title = 'Annual Max') + theme_bw()
# p05 <- ggplot(prec_ts, aes(x = date)) +
#   geom_line(aes(y = amed, color = dataset)) +
#   scale_color_manual(values = prec_color, guide = 'none') +
#   labs(y = NULL, title = 'Annual Median') + theme_bw()
# p06 <- ggplot(prec_ts, aes(x = date)) +
#   geom_line(aes(y = amean, color = dataset)) +
#   scale_color_manual(values = prec_color, guide = 'none') +
#   labs(y = NULL, title = 'Annual Average') + theme_bw()
# 
# # Aux legends
# gauge_legend <- ggplot(prec_ts[type == 'Gauge-based']) +
#   geom_line(aes(x = date, y = value, color = dataset)) +
#   scale_color_manual(values = prec_color, name = 'Gauge-based') +
#   theme_bw() +
#   theme(legend.text = element_text(size = 20),
#         legend.title = element_text(size = 24))
# model_legend <- ggplot(prec_ts[type == 'Model forcing']) +
#   geom_line(aes(x = date, y = value, color = dataset)) +
#   scale_color_manual(values = prec_color, name = 'Model forcing') +
#   theme_bw() +
#   theme(legend.text = element_text(size = 20),
#         legend.title = element_text(size = 24))
# reanalysis_legend <- ggplot(prec_ts[type == 'Reanalysis']) +
#   geom_line(aes(x = date, y = value, color = dataset)) +
#   scale_color_manual(values = prec_color, name = 'Reanalysis') +
#   theme_bw() +
#   theme(legend.text = element_text(size = 20),
#         legend.title = element_text(size = 24))
# satellite_legend <- ggplot(prec_ts[type == 'Satellite-based']) +
#   geom_line(aes(x = date, y = value, color = dataset)) +
#   scale_color_manual(values = prec_color, name = 'Satellite-based') +
#   theme_bw() +
#   theme(legend.text = element_text(size = 20),
#         legend.title = element_text(size = 24))
# 
# plot_grid(plot_grid(p01, p02, p03, p04, p05, p06, nrow = 3, ncol = 2),
#           plot_grid(get_legend(gauge_legend), get_legend(model_legend),
#                     get_legend(reanalysis_legend),
#                     get_legend(satellite_legend), ncol = 1),
#           ncol = 2, rel_widths = c(5, 0.75)) %>%
#   annotate_figure(left = textGrob("Precipitation in [mm]", rot = 90,
#                                   vjust = 1, gp = gpar(fontsize = 24)))

## ----eval=FALSE---------------------------------------------------------------
# plot_taylor(prec_ts, prec_obs, cols = prec_color)

## ----eval=FALSE---------------------------------------------------------------
# prec_8 <- prec_ts[dataset %in% c('gpcc-v2022', 'cpc-global', 'era5-land',
#                               'merra-2', 'gldas-noah-v2-0', 'gldas-clsm-v2-0',
#                               'gpcp-v3-2', 'mswep-v2-8')]
# prec8_colors <- setNames(c("#999999", "#E69F00", "#56B4E9", "#009E73",
#                            "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
#                          unique(prec_8$name))
# plot_taylor(prec_8, prec_obs, groups = 'seasons', cols = prec8_colors)

## ----eval=FALSE---------------------------------------------------------------
# prec_4 <- prec_ts[dataset %in% c('gpcc-v2022', 'merra-2', 'mswep-v2-8',
#                                  'gldas-noah-v2-0')] %>%
#   .[year(date) <= 2000, period := "1981 - 2000"] %>%
#   .[year(date) > 2000, period := "2001 - 2020"]
# 
# plot_density(prec_4) +
#   facet_grid(dataset + source ~ period) +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 24),
#         strip.text = element_text(size = 24),
#         strip.background = element_rect(fill = "white",
#                                         color = "black", linewidth = 1))

## ----eval=FALSE---------------------------------------------------------------
# prec_grid <- lapply(prec_data, function(x) x@file@name %in% c('gpcc-v2022',
#                                                               'merra-2',
#                                                               'mswep-v2-8',
#                                                               'gldas-noah-v2-0'))
# prec_grid <- prec_data[which(prec_grid == TRUE)]
# prec_grid <- foreach (dataset_count = 1:4, .combine = rbind) %do% {
#   dummie <- prec_grid[[dataset_count]]
#   dummie_name <- dummie@file@name
#   dummie <- tabular(dummie) %>%
#     .[year(date) <= 2000, period := "1981 - 2000\nMedian"] %>%
#     .[year(date) > 2000, period := "2001 - 2020\nMedian"]
#   dummie <- dummie[, .(value = median(value, na.rm = TRUE)),
#                    by = .(lon, lat, period)]
#   dummie$dataset <- dummie_name
#   dummie
# }
# 
# prec_grid[dataset %in% GAUGE_BASED, source := "Gauge-based"
#           ][dataset %in% FORCINGS, source := "Model forcing"
#             ][dataset %in% REANALYSES, source := "Reanalysis"
#               ][dataset %in% SATELLITE_BASED, source := "Satellite-based"]
# 
# 
# plot_map(prec_grid, timestamp = FALSE) +
#   facet_grid(dataset + source ~ period) +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 24),
#         strip.text = element_text(size = 24),
#         strip.background = element_rect(fill = "white",
#                                         color = "black", linewidth = 2))

