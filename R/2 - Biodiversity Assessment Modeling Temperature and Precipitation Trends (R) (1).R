#--------------------------------------------------------------------------------------

# Sources

# Pulling in files from "Biodiversity Assessment.R"
# Antao et al., 2020 - https://www.nature.com/articles/s41559-020-1185-7#Abs1
# Repo - https://github.com/lauraantao/Temp_Biodiv_Change
# Blowes et al., 2019 - https://www.science.org/doi/10.1126/science.aaw1620#elettersSection
# Repo - https://github.com/sChange-workshop/BioGeo-BioDiv-Change/tree/master

#---------------------------------------------------------------------------------------

# Load all libraries and functions
source('C:/Users/User/Documents/School/UCI/Projects/Biodiversity Assessment/Share/0 - Biodiversity Assessment Source and Functions (R).R')

#--------------------------------------------------------------------------------------

# DYNAMIC VARIABLES

# Filter out values that started before this year
earliestYearSample <- 1965

#--------------------------------------------------------------------------------------

# TEMPERATURE TREND ANALYSIS

# Obtained temp records from HadCRUT4
# HadSST3- marine SST (1 degree resolution)
# CRUTEM4- air temps on land (0.5 degree resolution)

# Jones, P. D. et al. Hemispheric and large-scale land-surface air
# temperature variations: an extensive revision and an update to 2010. J. Geophys.
# Res. Atmos. 117, D05127 (2012)., Harris, I., Jones, P. D., Osborn, T. J. & Lister,
# D. H. Updated high-resolution grids of monthly climatic observations—the CRU
# TS3.10 Dataset. Int. J. Climatol. 34, 623–642 (2014).

# Extracted monthly mean temp for monitoring period (Yearstart:Yearend)
# Estimated mean temp trends (generalized additive models + package mgcv)
# mgcv= temporally autocorrelated error structure
# Compared models using AIC
# With and without “month”
# Extracted linear slope from best model
# Summary of trend of annual temp change

#--------------------------------------------------------------------------------------

# LOAD BIOTIME DATA FROM "1 - Biodiversity Assessment Filtering"

# Filtered assemblages df with updated start and end years
# Isolate cell and STUDY ID columns
biotimeRenamedGriddedFiltereAssemblages <- read.csv('Raw Data/Biotime_Renamed_Gridded_12_Filtered_New_Assemblage_Start_End_Year_Metadata.csv') %>%
  dplyr::select(-1, -STUDY_ID) %>% 
  separate(assemblageID, into = c('STUDY_ID', 'cell'), remove = FALSE) %>% 
  mutate(STUDY_ID = as.numeric(STUDY_ID), cell = as.numeric(cell))

# Filtered data and cells
biotimeRenamedGridded <- vroom('Raw Data/Biotime_Renamed_Gridded_12.csv')
biotimeRenamedGriddedFiltered <- vroom('Raw Data/Biotime_Renamed_Gridded_12_Filtered.csv')

# Must keep in WGS until after the extraction of data
biotimeRenamedGriddedCells <- st_read('Raw Data/Biotime_Renamed_Gridded_12.shp')
biotimeRenamedGriddedFilteredCells <- st_read('Raw Data/Biotime_Renamed_Gridded_12_Filtered.shp')

# Get a bounding box of the cells to use for all maps
bbox <- biotimeRenamedGriddedFilteredCells %>%
  centerPolygonsDF() %>% 
  st_bbox()

# Calculate the percent of total data and total grid cells retained
nrow(biotimeRenamedGriddedFiltered) / nrow(biotimeRenamedGridded) * 100
nrow(biotimeRenamedGriddedFilteredCells) / nrow(biotimeRenamedGriddedCells) * 100

#--------------------------------------------------------------------------------------

# LOAD AND EXPORT TEMPERATURE AND PRECIPITATION DATA

# Load netcdfs, convert to rasters, reproject, and export
# Start index should be for 1965-01-01
loadData <- function(file, startIndex, endIndex, label) {

  # # Open  brick, select bands
  # ncBrick <- brick(file)[[startIndex:endIndex]]
  # 
  # # Convert to raster, select bands
  # # Rename from brick to keep the band names as dates
  # ncRast <- setNames(terra::rast(file)[[startIndex:endIndex]], names(ncBrick))

  # For ERSST, skip brick, above ncRast steps and add c('sst')
  # Inherit names from previous raster
  ncRast <- setNames(terra::rast(file, 'sst')[[startIndex:endIndex]], names(ncRastMarineHadSST))
  
  # Write as raster, reproject to fix inconsistencies in nc projections
  writeRaster(terra::project(ncRast, '+init=EPSG:4326 +lon_0='), paste('Raw Data/NC_Raster_', label, '.tif', sep = ''), overwrite = TRUE)

  return(ncRast)
}

# Most recent year in BioTIME (need up through that year)
max(biotimeRenamedGriddedFiltered$YEAR)

# Sandbox to find starting index for 1965-01-01 and ending index for 2018-12-01
rast('Raw Data/cru_ts4.07.1901.2022.tmp.dat.nc')[[769]]
rast('Raw Data/cru_ts4.07.1901.2022.tmp.dat.nc')[[1416]]

# Terrestrial temperature - get just the data dates needed to reduce computation
# https://catalogue.ceda.ac.uk/uuid/5fda109ab71947b6b7724077bf7eb753
ncRastTerrestrial <- loadData('Raw Data/cru_ts4.07.1901.2022.tmp.dat.nc', 769, 1416, 'Terrestrial_Temperature')

# Marine temperature - get just the data dates needed to reduce computation
# https://www.metoffice.gov.uk/hadobs/hadsst4/Raw Data/download.html
ncRastMarineHadSST <- loadData('Raw Data/HadSST.4.0.1.0_median.nc', 1381, 2028, 'Marine_Temperature_HadSST')

# See notes in loadData function for specifically dealing with ERSST data as it comes in individual nc files
# Marine temperature - get just the data dates needed to reduce computation
# https://www.ncei.noaa.gov/pub/Raw Data/cmb/ersst/v5/netcdf/
ncRastMarineERSSTFiles <- list.files(recursive = TRUE, pattern = '.nc', ignore.case = TRUE)[3:650]
ncRastMarineERSST <- loadData(ncRastMarineERSSTFiles, 1, 648, 'Marine_Temperature_ERSST')

plot(ncRastMarineERSST[[1]])

# Monthly precipitation - get just the data dates needed to reduce computation
# https://psl.noaa.gov/Raw Data/gridded/data.gpcc.html
ncRastPrecipitation <- loadData('Raw Data/precip.mon.total.0.25x0.25.v2020.nc', 889, 1536, 'Precipitation')

#--------------------------------------------------------------------------------------

# EXTRACT DATA FROM NAM CELLS

# Function to get values from all bands for NAM cells
extractData <- function(raster, label) {

  # Extract mean values at cell, cbind with cell vales
  # After loading geoTIFF, raster runs faster than brick for extract
  extractedRaw <- cbind(
    data.frame(cell = biotimeRenamedGriddedFilteredCells %>% pull(cell)),
    exact_extract(raster, biotimeRenamedGriddedFilteredCells, fun = 'mean')
  )

  # Remove mean. from exact extract column names to match the others
  columnNames <- colnames(extractedRaw)
  formattedColumnNames <- columnNames %>% str_replace('mean.', '')

  # Rename columns
  extractedRaw <- extractedRaw %>% rename_at(vars(columnNames), ~ formattedColumnNames)

  # Convert to df by cell and measurements at each time increment
  extracted <- extractedRaw %>%

    # To get the date from the layer name...
    pivot_longer(cols = 2:ncol(extractedRaw)) %>%
    mutate(name = str_replace(name, 'X', '')) %>%
    mutate(name = str_replace(name, '\\.', '-')) %>%
    mutate(name = str_replace(name, '\\.', '-')) %>%
    mutate(date = name) %>%

    # Break date into parts so you can filter by year later
    separate(date, into = c('year', 'month', 'day')) %>%
    rename(date = name) %>%
    mutate(year = parse_integer(year)) %>%
    mutate(month = as.factor(month)) %>%
    mutate(date = as.Date(date)) %>%

    # Filter to just the years of choice
    filter(year >= earliestYearSample)

  # EXPORT FOR QUICKER LOADING
  write.csv(extracted, paste('Raw Data/Biotime_Renamed_Gridded_12_Filtered_', label, '_Time_Series.csv', sep = ''))

  return(extracted)
}

# Extract values for each variable time series, see how long each takes
extractedTerrestrialTemperature <- extractData(ncRastTerrestrial, 'Terrestrial_Temperature')
extractedMarineTemperatureHadSST <- extractData(ncRastMarineHadSST, 'Marine_Temperature_HadSST')
extractedMarineTemperatureERSST <- extractData(ncRastMarineERSST, 'Marine_Temperature_ERSST')
extractedPrecipitation <- extractData(ncRastPrecipitation, 'Precipitation')

#--------------------------------------------------------------------------------------

# FILTER TO JUST CELLS WITH ADEQUATE DATA COVERAGE

# Cell must have at least this percent of total monthly data to be kept
monthPercentThreshold <- 0.95

# Get the number of monthly layers by cell, filter those with less than threshold percent
filterMonthCoverage <- function(df, raster) {
  
  # Get counts per cells of monthly observations as a percent of total available
  filteredMonthCoverageIDs <- df %>%
    filter(!is.na(value)) %>% 
    group_by(cell) %>%
    count() %>% 
    mutate(percentMonths = (n / (dim(raster)[3]))) %>% 
    filter(percentMonths >= monthPercentThreshold) %>% 
    pull(cell)
  
  # Filter and join to have lat long data
  filteredMonthCoverage <- df %>% filter(cell %in% filteredMonthCoverageIDs)
  
  return(filteredMonthCoverage)
}

# Filter to cells and data with monthly coverage
extractedTerrestrialTemperatureFiltered <- filterMonthCoverage(extractedTerrestrialTemperature, ncRastTerrestrial)
extractedMarineTemperatureFilteredHadSST <- filterMonthCoverage(extractedMarineTemperatureHadSST, ncRastMarineHadSST)
extractedMarineTemperatureFilteredERSST <- filterMonthCoverage(extractedMarineTemperatureERSST, ncRastMarineERSST)
extractedPrecipitationFiltered <- filterMonthCoverage(extractedPrecipitation, ncRastPrecipitation)

#--------------------------------------------------------------------------------------

# MODEL TRENDS FOR EACH ASSEMBLAGE ACCORDING TO START AND END YEAR

trendData <- function(df, label) {

  # Create a blank data frame to append data after each model in the loop
  blank <- data.frame()

  # Get distinct cells from the extracted dataset
  dfCells <- df %>%
    distinct(cell) %>%
    pull(cell)

  # Get assemblages that took place in those cells
  dfAssemblages <<- biotimeRenamedGriddedFiltereAssemblages %>% filter(cell %in% dfCells)

  # Run temperature model for each assemblage and collate results
  for(row in 1:nrow(dfAssemblages)) {

    # Try because sometimes models don't work any throw error
    try({
      # Get the assemblage based on the row of all of the assemblages
      assemblage <- dfAssemblages %>% slice(row)

      # Get the ID of the cell as a number, filter sf to the cell (Should only be one)
      assemblageCellID <- assemblage %>% pull(cell) %>% as.numeric()
      assemblageCell <- biotimeRenamedGriddedFilteredCells %>% filter(cell %in% assemblageCellID)

      # Spot check for one
      print(nrow(assemblageCell))

      # Get the updated start and end year
      assemblageStartYear <- assemblage %>% pull(START_YEAR)
      assemblageEndYear <- assemblage %>% pull(END_YEAR)

      # Filter extracted values to just those within the cell and years of interest
      assemblageExtracted <- df %>%

        # Filter to the cell
        filter(cell == assemblageCellID) %>%

        # Filter to years of study
        filter(year >= assemblageStartYear) %>%
        filter(year <= assemblageEndYear)

      # Run models with and without the random effect of month using the study and point of interest
      modelRE <- glmmTMB(
        value ~ year + (1|month),
        #value ~ year + as.integer(month),
        data = assemblageExtracted %>% mutate(month = as.factor(month)),
        family = 'gaussian'
      )

      modelNoRE <- glmmTMB(
        value ~ year,
        data = assemblageExtracted,
        family = 'gaussian'
      )

      # Get the AICs for both models
      aicModelRE <- AIC(modelRE)
      aicModelNoRE <- AIC(modelNoRE)

      # If there is not an AIC, it will cause an error in the if else below
      # Set NA AICs to a very high number so they will always be higher
      if(is.na(aicModelRE)) {aicModelRE <- 1e13}
      if(is.na(aicModelNoRE)) {aicModelNoRE <- 1e13}

      # Chose the lower one as the final model
      if(aicModelRE > aicModelNoRE) {

        model <- modelNoRE
        modelType <- 'Model without RE'
        aicModel <- aicModelNoRE

      } else if(aicModelRE < aicModelNoRE) {

        model <- modelRE
        modelType <- 'Model with RE'
        aicModel <- aicModelRE

      } else if(aicModelRE == 1e13 & aicModelNoRE == 1e13) {

        model <- modelRE # Just take the model with month if both NA
        modelType <- 'Model with RE NA AIC'
        aicModel <- NA
      }

      # Get the stats from the chosen model
      modelSummary <- summary(model)

      modelInterceptEstimate <- modelSummary$coefficients[[1]][1,1]
      modelYearEstimate <- modelSummary$coefficients[[1]][2,1]
      modelYearSE <- modelSummary$coefficients[[1]][2,2]
      modelYearZ <- modelSummary$coefficients[[1]][2,3]
      ModelYearP <- modelSummary$coefficients[[1]][2,4]

      # Compile into one row
      assemblageModel <- assemblage %>%
        mutate(
          modelType = modelType,
          aicModel = aicModel,
          modelInterceptEstimate = modelInterceptEstimate,
          modelYearEstimate = modelYearEstimate,
          modelYearSE = modelYearSE,
          modelYearZ = modelYearZ,
          ModelYearP = ModelYearP
        )

      # Add to running data frame
      blank <- rbind(blank, assemblageModel)

      print(paste('Assemble ID #', row))
    })
  }

  # Export trend analysis
  write.csv(blank, paste('Raw Data/Biotime_Renamed_Gridded_12_Filtered_', label, '_Trends.csv', sep = ''))

  return(blank)
}

assemblageTerrestrialTemperatureTrends <- trendData(extractedTerrestrialTemperatureFiltered, 'Terrestrial_Temperature')
assemblageMarineTemperatureTrendsHadSST <- trendData(extractedMarineTemperatureFilteredHadSST, 'Marine_Temperature_HadSST')
assemblageMarineTemperatureTrendsERSST <- trendData(extractedMarineTemperatureFilteredERSST, 'Marine_Temperature_ERSST')
assemblagePrecipitationTrends <- trendData(extractedPrecipitationFiltered, 'Precipitation')

#--------------------------------------------------------------------------------------