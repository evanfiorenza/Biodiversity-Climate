#--------------------------------------------------------------------------------------

# Set the working directory
setwd('C:/Users/User/Documents/School/UCI/Projects/Biodiversity Assessment/Share')

# Load needed libraries or install them if you do not already have them on your computer
library(vroom)
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(viridis)
library(RColorBrewer)
library(scales)
library(patchwork)
library(glmmTMB)
library(vegan)
library(units)
library(ncdf4)
library(raster)
library(mgcv)
library(rasterVis)
library(scales)
library(exactextractr)

# Install and load BioTIMEr
# I had to download the zip from https://github.com/bioTIMEHub/BioTIMEr then install manually
#devtools::install_local('BioTIMEr-main.zip')
library(BioTIMEr)

# Set a seed for reproducible results
set.seed(41)

#--------------------------------------------------------------------------------------

# FUNCTIONS TO HELP WITH DATA PROCESSING

# Central meridian for Robinson reprojection mapping (center of map at this longitude in degrees)
lon_0 <- -100

# Function to reproject and center points so all are in same projected coordinate systems (for mapping)
centerPointsDF <- function(df) {
  
  return(df %>%
           st_cast('POINT') %>% 
           st_make_valid() %>%
           st_wrap_dateline() %>% 
           st_transform(crs = paste0('+proj=robin +lon_0=', lon_0))
  )
}

# Function to reproject and center polygons so all are in same projected coordinate systems (for mapping)
centerPolygonsDF <- function(df) {
  
  return(df %>%
           st_cast('MULTIPOLYGON') %>% 
           st_cast('POLYGON') %>% 
           st_make_valid() %>%
           st_wrap_dateline() %>% 
           st_break_antimeridian(lon_0 = lon_0) %>% 
           st_transform(crs = paste0('+proj=robin +lon_0=', lon_0))
  )
}

# Function to center based on WGS (for grid cell analysis)
centerPolygonsDFWGS <- function(df) {
  
  return(df %>%
           st_cast('MULTIPOLYGON') %>% 
           st_cast('POLYGON') %>% 
           st_make_valid() %>%
           st_wrap_dateline() %>% 
           st_transform(crs = 4326)
  )
}

#--------------------------------------------------------------------------------------

# HELPER OBJECTS

# Load global countries to use to filter BioTIME to NAM


worldmap <- rnaturalearth::ne_download(scale = 'large',
                                       type = "countries",
                                       category = "cultural",
                                       destdir = tempdir(),
                                       load = TRUE,
                                       returnclass = "sf")


countriesRobinson <- worldmap %>%
  st_transform(4326) %>%
  centerPolygonsDF()

# Load global countries to use to filter BioTIME to NAM

countriesWGS <- worldmap %>%
  st_transform(4326) %>%
  centerPolygonsDFWGS()

#--------------------------------------------------------------------------------------

# MAKE SCALES FOR DIFFERENT MAPS

# Global temperature viz
temperatureWorldViz <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')),
  breaks = c(-20, 0, 30),
  limits = c(-20, 30),
  labels = c(-20, 0, 30),
  oob = squish,
  na.value = NA)

# Temperature trend viz
temperatureTrendViz <- scale_color_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')),
  breaks = c(-0.1, 0, 0.1),
  limits = c(-0.1, 0.1),
  labels = c(-0.1, 0, 0.1),
  oob = squish,
  na.value = NA)

# Temperature trend viz to match Antao et al., 2020
temperateTemperatureTrendVizAntao <- scale_color_viridis_c(
  option = 'magma',
  breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6),
  limits = c(-0.4, 0.6),
  labels = c(-0.4, '', 0, '', 0.4, 0.6),
  oob = squish,
  na.value = NA)

# Standardized temperatures viz
temperatureStandardizedViz <- scale_color_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')),
  breaks = c(-1, 0, 1),
  limits = c(-1, 1),
  labels = c(-1, 0, 1),
  oob = squish,
  na.value = NA)

temperatureStandardizedVizFill <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')),
  breaks = c(-1, 0, 1),
  limits = c(-1, 1),
  labels = c(-1, 0, 1),
  oob = squish,
  na.value = NA)

# Global precipitation viz
precipitationWorldViz <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlGn')),
  breaks = c(0, 100, 200),
  limits = c(0, 200),
  labels = c(0, 100, 200),
  oob = squish,
  na.value = NA)

precipitationTrendViz <- scale_color_gradientn(
  colours = (RColorBrewer::brewer.pal(11, 'RdYlGn')),
  breaks = c(-1, 0, 1),
  limits = c(-1, 1),
  labels = c(-1, 0, 1),
  oob = squish,
  na.value = NA)

#--------------------------------------------------------------------------------------