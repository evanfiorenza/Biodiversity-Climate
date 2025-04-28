#--------------------------------------------------------------------------------------

# Sources (papers, Git hubs, data access, etc.)
# Antao et al., 2020 - https://www.nature.com/articles/s41559-020-1185-7
# Github for ^ https://github.com/lauraantao/Temp_Biodiv_Change
# Blowes et al., 2019 https://www.science.org/doi/epdf/10.1126/science.aaw1620
# BioTIME - https://biotime.st-andrews.ac.uk/
# BoTIMEr Gitgub - https://github.com/bioTIMEHub/BioTIMEr

#--------------------------------------------------------------------------------------

# Load all libraries and functions
source('R/0 - Biodiversity Assessment Source and Functions (R).R')

#--------------------------------------------------------------------------------------

# DYNAMIC VARIABLES

# "By default, the function separates studies into hexagonal cells of ~96 km2 (`res=12`), as this resolution
# was found to be the most appropriate when working on the whole BioTIME database. Yet sometimes, particularly
# when working with a subset of BioTIME, you may want to explore alternative grid resolutions."   
griddingResolution <- 12

# Filter out samples that started before this year
earliestYearSample <- 1965

# Minimum number of years per time series (5 per Antao et al., 2020)
minYearsTimeSeries <- 5

#--------------------------------------------------------------------------------------

# LOAD BIOTIME DATA, PROCESS AND FILTER TO JUST THE DATA WE WANT

# Use vroom to load quickly and join the metadata based on STUDY_ID (need REALM and START_YEAR for filtering)
# Can join the metadata here if you do not need to sub out any of the study ID points
biotime <- vroom('Raw Data/BioTIMEQuery_24_06_2021.csv')

# BioTIME and metadata from https://biotime.st-andrews.ac.uk/download.php
# Load the metadata CSV
biotimeMetadata <- vroom('Raw Data/BioTIMEMetadata_24_06_2021.csv')

# Get a unique list of genus and species by species ID
biotimeSpeciesID <- biotime %>% distinct(ID_SPECIES, GENUS_SPECIES, GENUS, SPECIES)

# Swap out STUDY_ID 169
# Augment columns to match that from the website version
study169Server <- vroom('Raw Data/BioTIME_study169_fromserver.csv') %>%

  # Rename abundance an biomass columns
  rename(sum.allrawdata.ABUNDANCE = ABUNDANCE) %>%
  rename(sum.allrawdata.BIOMASS = BIOMASS) %>%
  rowid_to_column('...1') %>%

  # Remove unnecessary columns
  dplyr::select(-DEPTH, -ID_ALL_RAW_DATA, -newID) %>%

  # Join with metadata
  left_join(biotimeSpeciesID, by = c('ID_SPECIES'))

# Reorder to get the same order as web version
study169Server <- study169Server[,colnames(biotime)]

# Remove original 169, substitute with server 169
# Join with metadata
biotime <- biotime %>%
  filter(STUDY_ID != 169) %>%
  rbind(study169Server) %>%
  left_join(biotimeMetadata, by = c('STUDY_ID'))

# Get a list of unique survey sites using latitude, longitude, and study ID
# Convert to a spatial layer and center
biotimeSites <- biotime %>%
  group_by(LATITUDE, LONGITUDE, STUDY_ID) %>%
  count() %>%
  st_as_sf(coords = c('LONGITUDE', 'LATITUDE'), crs = 4326, remove = FALSE) %>%
  centerPointsDF()

# Map out BioTIME sites
ggplot() +
  geom_sf(data = countriesRobinson) +
  geom_sf(data = biotimeSites, color = 'red', size = 0.1) +
  theme_classic()

# Save BioTIME sites map
ggsave('Figures/BioTIME_Sites.png', width = 4000, height = 2000, units = 'px')

#--------------------------------------------------------------------------------------

# TROUBLESHOOT MISMATCHING COLUMN NAMES BETWEEN DOWNLOADED BIOTIME DATASET AND REQUIRED COLUMNS FOR GRIDDING() FUNCTION

# This will produce an error message and tell you the needed columns
biotimeGridded <- gridding(meta = biotimeMetadata, btf = biotime, res = griddingResolution, resByData = FALSE)

# Error in gridding(meta = biotimeMetadata, btf = biotime, res = griddingResolution,  :
# Assertion on 'colnames(btf)' failed: Colnames must include the elements
# {'valid_name','STUDY_ID','DAY','MONTH','YEAR','LATITUDE','LONGITUDE','ABUNDANCE','BIOMASS','SAMPLE_DESC','PLOT','resolution','taxon'}
# ,but is missing elements {'valid_name','ABUNDANCE','BIOMASS','resolution','taxon'}.

# Trying to find equivalent columns for the downloaded subset
columnsNeeded <- c(
  'valid_name','STUDY_ID','DAY','MONTH','YEAR','LATITUDE','LONGITUDE','ABUNDANCE',
  'BIOMASS','SAMPLE_DESC','PLOT','resolution','taxon'
)

# Load data from BioTIMEr package to see what data is provided in those required columns
data('BTsubset_data')
data('BTsubset_meta')

# Column names for data
colnames(biotime) %>% sort()
colnames(BTsubset_data) %>% sort()

# Column names for metadata
colnames(biotimeMetadata) %>% sort()
colnames(BTsubset_meta) %>% sort()

# Rectify mismatches - try to find if equivalent of column exists in downloaded BioTIME data frame
biotimeObservations <- BTsubset_data %>%
  group_by(STUDY_ID) %>%
  count()

# Get a study from the provided subset with minimal observations so it's easier to investigate
biotimeObservations %>% arrange(n)

btSubset <-  BTsubset_data %>% filter(STUDY_ID == 53)
downloadedSubset <- biotime %>% filter(STUDY_ID == 53)

# Check they have the same number of observations
nrow(btSubset) == nrow(downloadedSubset)

# valid_name = genus species (sp for species if just a genus level)
# ABUNDANCE = count of species at location
# BIOMASS = total biomass of species at location
# resolution = species or genus level
# taxon = mammals, plants, etc.

# Rename downloaded subset columns to match bt subset columns
downloadedSubsetColumns <- downloadedSubset %>%

  # bt subset using sp rather than spec for Genus sp when just genus level
  mutate(SPECIES = ifelse(SPECIES == 'spec', 'sp', SPECIES)) %>%

  # Create valid_name column using ^ nomenclature
  mutate(valid_name = paste(GENUS, SPECIES)) %>%

  # If species is just sp, resolution is genus, otherwise species
  mutate(resolution = ifelse(SPECIES == 'sp', 'genus', 'species')) %>%

  # Rename abundance, biomass, and taxa to match
  rename(ABUNDANCE = sum.allrawdata.ABUNDANCE) %>%
  rename(BIOMASS = sum.allrawdata.BIOMASS) %>%
  rename(taxon = TAXA) %>%

  # Select just the columns needed for the gridding() function
  dplyr::select(columnsNeeded)

# Select just the needed columns
btSubsetColumns <- btSubset %>% dplyr::select(columnsNeeded)

df1 <- downloadedSubsetColumns %>% arrange(valid_name, YEAR)
df2 <- btSubsetColumns %>%  arrange(valid_name, YEAR)

# Compare results
# Makes sense DAY and MONTH don't match - downloaded subset is populated, bt subset is not
# Makes sense BIOMASS doesn't match - downloaded subset is 0, bt subset is NAs
# Makes sense SAMPLE_DESC (sample description) doesn't match, is more detailed for downloaded subset
# Make sense PLOT doesn't match, not populated for either
(df1 == df2) & !is.na(df1) & !is.na(df2)

# Why isn't valid_name matching?

# These two names are being disagreed on (both the Southern red-backed vole)
# Clethrionomys gapperi - Downloaded subset
# Myodes gapperi - bt subset

# https://en.wikipedia.org/wiki/Clethrionomys
# "Clethrionomys is a genus of small, slender voles.
# In recent years the genus name was changed to Myodes,
# however a 2019 paper found that Myodes was actually a junior synonym for Lemmus, thus making it unusable.
# As such, Clethrionomys is re-established as the proper genus name.[2]

# So, they are the same species but the downloaded subset has the most recently accepted name (great!)

# Why isn't ABUNDANCE matching? Same numbers, must be a format thing?
position <- 1

# Get the ABUNDANCE values at a row where they do not match
df1PositionAbundance <- df1 %>% slice(position) %>% pull(ABUNDANCE)
df2PositionAbundance <- df2 %>% slice(position) %>% pull(ABUNDANCE)

# Look at values
df1PositionAbundance
df2PositionAbundance

# Appear the same, but not equal
df1PositionAbundance == df2PositionAbundance

# Print full values - downloaded subset has much higher precision!!
print(df1PositionAbundance, digits = 15)
print(df2PositionAbundance, digits = 15)

# Round to see if ABUNDANCE now matches
df1 <- downloadedSubsetColumns %>%
  arrange(valid_name, YEAR) %>%
  mutate(ABUNDANCE = round(ABUNDANCE, digits = 3))

df2 <- btSubsetColumns %>%
  arrange(valid_name, YEAR)  %>%
  mutate(ABUNDANCE = round(ABUNDANCE, digits = 3))

# Yes, they do!
(df1 == df2) & !is.na(df1) & !is.na(df2)

# Difference have been investigated and explained, necessary columns selected

#--------------------------------------------------------------------------------------

# FIX THE COLUMNS TO MATCH THOSE NEEDED FOR THE GRIDDING() FUNCTION

biotimeRenamed <- biotime %>%

  # bt subset using sp rather than spec for Genus sp when just genus level
  mutate(SPECIES = ifelse(SPECIES == 'spec', 'sp', SPECIES)) %>%

  # Create valid_name column using ^ nomenclature
  mutate(valid_name = paste(GENUS, SPECIES)) %>%

  # If species is just sp, resolution is genus, otherwise species
  mutate(resolution = ifelse(SPECIES == 'sp', 'genus', 'species')) %>%

  # Rename abundance, biomass, and taxa to match
  rename(ABUNDANCE = sum.allrawdata.ABUNDANCE) %>%
  rename(BIOMASS = sum.allrawdata.BIOMASS) %>%
  rename(taxon = TAXA) %>%

  # Remove the metadata columns because gridding() joins with metadata by selecting only what is needed
  # That then make repeat columns end with .x or .y and throws an error saying columns don't exist
  dplyr::select(columnsNeeded)

# #--------------------------------------------------------------------------------------

# GRID THE DATA - ALL BIOTIME DATA (ADDS CELL ID TO DATASET)

# Manually setting the resolution, change resByData to TRUE to have it suggest best size
# Turn into data frame because it's easier
# Separate - STUDY_ID just repeats but need cell by itself to match up with hexagons shapefile
biotimeRenamedGridded <- gridding(meta = biotimeMetadata, btf = biotimeRenamed, res = griddingResolution, resByData = FALSE) %>%
  mutate(YEAR = as.integer(YEAR)) %>%
  as.data.frame() %>%
  separate(assemblageID, into = c('STUDY_ID_REPEAT', 'cell'), remove = FALSE)

# Export for easier loading
write.csv(biotimeRenamedGridded, paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution, '.csv', sep = ''))

# See differences in column names between vignette and downloaded data
setdiff(colnames(BTsubset_data), colnames(biotimeRenamed))

# Make sure the number of rows is the same
nrow(biotimeRenamed) - nrow(biotimeRenamedGridded)

# Make sure the studies match
setdiff(
  biotimeRenamed %>% distinct(STUDY_ID) %>% pull(STUDY_ID),
  biotimeRenamedGridded %>% distinct(STUDY_ID) %>% pull(STUDY_ID)
)

#--------------------------------------------------------------------------------------

# EXTRACTING HEXAGONS FROM GRIDDING() FUNCTION (GET A SHAPEFILE OF GRIDS WITH CELLS IDS)

# Run the gridded function again
gridded <- gridding(meta = biotimeMetadata, btf = biotimeRenamed, res = griddingResolution, resByData = FALSE) %>%
  mutate(YEAR = as.integer(YEAR)) %>%
  as.data.frame()

# Set the grid resolution for the shapefile
dgg <- dggridR::dgconstruct(res = griddingResolution)

# Construct the hexagons using the appropriate grid resolution
biotimeRenamedGriddedCells <- dggridR::dgcellstogrid(dgg, gridded$cell) %>%
  rename(cell = seqnum) %>%
  centerPolygonsDFWGS()

# Export for easier loading, uploading to Google Earth Engine (prefers data be in WGS)
st_write(biotimeRenamedGriddedCells, paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution, '.shp', sep = ''))

#--------------------------------------------------------------------------------------

# CREATE A NORTH AMERICA (NAM) FILTER DEFINED AS...
# CLOSEST COUNTRY IS NAM AND/OR INTERSECTS A MARINE OR TERRESTRIAL ECOREGION THAT INTERSECTS A NAM COUNTRY

# Collate a list of countries in NAM to use to filter BioTIME data based on nearest country
countriesNAMList <- c('Canada', 'Mexico', 'United States of America')

# Filter countries to just those in NAM
countriesNAM <- countriesWGS %>% filter(GEOUNIT %in% countriesNAMList)

# NAM countries are gray, non-NAM countries are red so it's easy to see if you missed any
ggplot() +
  geom_sf(data = countriesRobinson %>% st_wrap_dateline(), color = 'red', fill = 'red') +
  geom_sf(data = countriesNAM %>% st_wrap_dateline()) +
  theme_classic()

# Save NAM countries of global map
ggsave('Figures/Countries_NAM.png', width = 4000, height = 2000, units = 'px')

# Load ecoregion units for terrestrial realms
terrestrialEcoregions <- st_read('GIS/tnc_terr_ecoregions.shp') %>%
  st_transform(4326) %>%
  centerPolygonsDFWGS()

# Load ecoregion units for marine realms
marineEcoregions <- st_read('GIS/meow_ecos.shp') %>%
  st_transform(4326) %>%
  centerPolygonsDFWGS()

# Filter MEOWs to just countries that overlap NAM countries
terrestrialEcoregionsNAM <- st_filter(terrestrialEcoregions, countriesNAM, .predicate  = st_intersects)
marineEcoregionsNAM <- st_filter(marineEcoregions, countriesNAM, .predicate  = st_intersects)

# Make plots of the overlapping ecoregions compared to all ecoregions
terrestrialEcoregionsNAMPlot <- ggplot() +
  geom_sf(data = countriesWGS %>% st_wrap_dateline()) +
  geom_sf(data = terrestrialEcoregionsNAM %>% st_wrap_dateline(), fill = 'green', color = 'black', linewidth = 0.1) +
  ggtitle('Terrestrial Ecoregions Intersecting NAM Countries\n') +
  theme_bw()

marineEcoregionsNAMPlot <- ggplot() +
  geom_sf(data = countriesWGS %>% st_wrap_dateline()) +
  geom_sf(data = marineEcoregionsNAM %>% st_wrap_dateline(), fill = 'cyan', color = 'black', linewidth = 0.1) +
  ggtitle('Marine Ecoregions Intersecting NAM Countries\n') +
  theme_bw()

# Plot and save
(terrestrialEcoregionsNAMPlot / marineEcoregionsNAMPlot)

ggsave('Figures/Terrestrial_Marine_Ecoregions_World_NAM.png', width = 3000, height = 3500, units = 'px')

# Combine to get overall coverage and spatial extent of NAM study area
ecoregionsNAM <- st_union(rbind(terrestrialEcoregionsNAM %>% dplyr::select(geometry), marineEcoregionsNAM %>% dplyr::select(geometry)))

ggplot() +
  geom_sf(data = countriesWGS) +
  geom_sf(data = ecoregionsNAM, fill = 'gold', color = 'black', linewidth = 0.1) +
  ggtitle('Ecoregions Intersecting NAM Countries\n') +
  theme_bw()

ggsave('Figures/Ecoregions_World_NAM.png', width = 3000, height = 3500/2, units = 'px')

# Get the centroid of each cell to use to find the nearest country
cellsPoints <- biotimeRenamedGriddedCells %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(LONGITUDE = st_coordinates(centroid)[,1], LATITUDE = st_coordinates(centroid)[,2])

# Using the countries shapefile, find the closest country to each of the BioTIME cells
# Select just needed columns
cellsClosestCountry <- st_join(cellsPoints, countriesWGS, join = st_nearest_feature) %>%
  as.data.frame() %>%
  dplyr::select(COUNTRY, cell)

# Get a list of the cells whose closest country is in NAM
cellsClosestCountryNAM <- cellsClosestCountry %>%
  filter(COUNTRY %in% countriesNAMList) %>%
  pull(cell)

# Get the number of nearest cells per country
cellsClosestCountrySF <- countriesWGS %>%
  left_join(cellsClosestCountry %>% as.data.frame() %>% group_by(COUNTRY) %>% count(), by = c('COUNTRY'))

# Map of nearest study sites per country
ggplot() +
  geom_sf(data = cellsClosestCountrySF, aes(fill = log10(n))) +
  theme_classic() +
  scale_fill_viridis_c(name = 'Log10 n Nearest Grid Cells', option = 'magma')

ggsave('Figures/BioTIME_Cells_Nearest_Countries_Global.png', width = 4000, height = 2000, units = 'px')

# Filter to just cells within the NAM ecoregions, closest country is NAM
cellsEcoregionsNAM <- st_filter(x = biotimeRenamedGriddedCells, y = ecoregionsNAM, .predicate = st_intersects)
cellsClosestNAM <- biotimeRenamedGriddedCells %>% filter(cell %in% cellsClosestCountryNAM)

# Filter of cells from both filters (pass either/both filters) by pulling IDs from both filters
cellsEcoregionsNAMIDs <- cellsEcoregionsNAM %>% pull(cell)
cellsClosestNAMIDs <- cellsClosestNAM %>% pull(cell)

# Filter to cells to either of the two filter definition
cellsEcoregionsClosestNAM <- biotimeRenamedGriddedCells %>%
  filter(cell %in% cellsEcoregionsNAMIDs | cell %in% cellsClosestNAMIDs)

# Number of cells remaining
nCellsEcoregionsClosestNAM <- nrow(cellsEcoregionsClosestNAM)

# Map all cells and the cells that passed the criteria for NAM
ggplot() +
  geom_sf(data = countriesWGS) +
  geom_sf(data = biotimeRenamedGriddedCells, fill = 'red', color = 'red') +
  geom_sf(data = cellsEcoregionsClosestNAM, fill = 'black', color = 'black') +
  ggtitle(paste('Remaining Cells (', nCellsEcoregionsClosestNAM, ') After Filtering by Closest to NAM Country or NAM Ecoregion\n', sep = '')) +
  theme_classic()

ggsave('Figures/BioTIME_Cells_Ecoregions_NAM_Closest_NAM.png', width = 3000, height = 1750, units = 'px')

# Get and export a csv of cells that you can pipe into Google Earth Engine and run cell-based analyses/data grabbing
# Easier to have all cells and then just subset from these lists rather than use filtered shapefiles below
write.csv(cellsEcoregionsClosestNAM %>% as.data.frame() %>% distinct(cell), paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution,'_Ecoregions_NAM_Closest_NAM_Cells.csv', sep = ''), row.names = FALSE)

# Export cells
st_write(cellsEcoregionsClosestNAM, paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution,'_Ecoregions_NAM_Closest_NAM.shp', sep = ''), append = FALSE)

#--------------------------------------------------------------------------------------

# CREATE SERIES OF INDIVIDUAL FILTERS TO JUST DESIRED DATA

# Get NAM cells
filterNAM <- biotimeRenamedGridded %>%
  filter(cell %in% (cellsEcoregionsClosestNAM %>% pull(cell))) %>%
  distinct(cell) %>%
  pull(cell)

#--------------------------------------------------------------------------------------

# Filter to just samples more recent than earliest year
# See how much data retained
filterEarliestYear <- biotimeRenamedGridded %>% 
  filter(YEAR >= 1965) %>% 
  distinct(SAMPLE_DESC) %>% 
  pull(SAMPLE_DESC)

#--------------------------------------------------------------------------------------

# CUMULATIVE FILTERING

# ORDER OF FILTERS
# f0 = Start
# f1 = NAM - Cell-wise
# f2 = Only since 1965 - Sample-wise
# f3 = Years in time series - Time series-wise

f0 <- biotimeRenamedGridded
f1 <- f0 %>% filter(cell %in% filterNAM)
f2 <- f1 %>% filter(SAMPLE_DESC %in% filterEarliestYear)

# Get assemblageIDs with minYearsTimeSeries years or more
f3Filter <- f2 %>%
  group_by(assemblageID) %>%
  summarise(n = n_distinct(YEAR)) %>%
  filter(n >= minYearsTimeSeries) %>%
  pull(assemblageID)

f3 <- f2 %>% filter(assemblageID %in% f3Filter)

#--------------------------------------------------------------------------------------

# GET REMAINING FILTERED DATA AND CELLS, CREATE NEW ASSEMBLAGE ID YEAR START YEAR END, EXPORT

biotimeRenamedGriddedFiltered <- f3

# Get the names of the remaining assemblage IDs
biotimeRenamedGriddedFilteredIDs <- biotimeRenamedGriddedFiltered %>%
  distinct(cell) %>%
  pull(cell)

# Filter the grid to just those cells
biotimeRenamedGriddedFilteredCells <- biotimeRenamedGriddedCells %>%
  filter(cell %in% biotimeRenamedGriddedFilteredIDs)

# Percent of BioTIME data retained after filters
paste(
  'Perfect of BioTIME data retained after all filters: ',
  round((nrow(biotimeRenamedGriddedFiltered) / nrow(biotimeRenamedGridded)) * 100, 2), '%', sep = ''
)

# Based on the filtered cells, get assemblages and start year end year, might have changed during filtering
# Need this for years specific trend analysis for temperature and precipitation
biotimeRenamedGriddedFilteredAssemblageIDStartYearEndYear <- biotimeRenamedGriddedFiltered %>%
  group_by(assemblageID, STUDY_ID) %>%
  summarise(
    START_YEAR = min(YEAR),
    END_YEAR = max(YEAR)
  )

# Join with original metadata, except replace with new filtered start and end year, to study-wide info
biotimeRenamedGriddedFilteredAssemblageIDStartYearEndYearMetadata <- biotimeRenamedGriddedFilteredAssemblageIDStartYearEndYear %>%
  left_join(biotimeMetadata %>% dplyr::select(-START_YEAR, -END_YEAR), by = c('STUDY_ID'))

# Export
st_write(biotimeRenamedGriddedFilteredCells, paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution, '_Filtered.shp', sep = ''), append = FALSE)
write.csv(biotimeRenamedGriddedFiltered, paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution, '_Filtered.csv', sep = ''))
write.csv(biotimeRenamedGriddedFilteredAssemblageIDStartYearEndYear, paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution, '_Filtered_New_Assemblage_Start_End_Year.csv', sep = ''))
write.csv(biotimeRenamedGriddedFilteredAssemblageIDStartYearEndYearMetadata, paste('Raw Data/Biotime_Renamed_Gridded_', griddingResolution, '_Filtered_New_Assemblage_Start_End_Year_Metadata.csv', sep = ''))

#--------------------------------------------------------------------------------------