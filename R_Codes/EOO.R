### This code replicates the construction and estimation of the geographic ranges (polynomial methods only) 
# for fish of the Curimatidae family from the paper titled: 
# "Comparing methods for estimating geographic ranges in freshwater fishes"
# by Valencia-Rodríguez et al. Submitted to Freshwater Biology

#--------------------------------------------------------------#
####                1.1 Optimal alpha value                 ####
#--------------------------------------------------------------#

# Load packages
library(readr)
library(dplyr)
library(rangeBuilder)

# Set your working directory
setwd("~/yourpath")

### Read and prepare data ###
# Read the clean species occurrence database. 
# This database only requires three columns: (i) species name, (ii) longitude, and (iii) latitude
data_ocurr <- read_csv("curimatidae_occ_2023.csv")

# Note: As our analysis involves comparing methods, we selected only species with 20 or more occurrence records, 
# since this is the minimum number of records used to generate the SDMs
# If you require range estimation with fewer than three occurrence records, 
# consult the work of Andrade-Garcia et al. (2021) https://doi.org/10.1111/geb.13299, 
# or their R code for the analyses available at https://github.com/fabro/Poeciliidae_LDG

# Selecting species with 20 or more occurrence records
spp_counts <- data_ocurr %>% 
  group_by(species) %>% # Column with the name of the species
  summarise(records = n()) %>%
  filter(records >= 19)
# Filter the original data to include only species with at least 20 records
filter_data <- data_ocurr %>% 
  filter(species %in% spp_counts$species)
# Create the final list and assign names to each item
spp.list <- filter_data %>% 
  group_split(species) %>%
  setNames(spp_counts$species)

# Errors might occur if the occurrences are too close geographically or if two or more records share 
# the same value on the longitude or latitude. Therefore, we adjust the decimal points of the longitude 
# and latitude coordinates to avoid potential failures in the 'ahull' function
spp.list.2 <- list()
for (i in 1:length(spp.list)) {
  spp.list.2[[i]] <- spp.list[[i]] %>%
    dplyr::select("longitude", "latitude") %>%
    mutate(longitude=ifelse(duplicated(spp.list[[i]]$longitude), longitude + rnorm(1, mean = 0, sd = 0.0001), longitude)) %>%
    mutate(latitude=ifelse(duplicated(spp.list[[i]]$latitude), latitude+rnorm(1, mean=0, sd=0.0001), latitude)) %>%  
    hsi::clean_dup(longitude = "longitude", latitude = "latitude", threshold = 0.041665)
  names(spp.list.2)[i] <- names(spp.list)[i]
}

# Note: The alphas are computed using the rangeBuilder package (Davis Rabosky et al., 2016). 
# However, some species may trigger errors during the execution of the 'getDynamicAlphaHull' function, 
# due to the arrangement of the records in the geography.
# The subsequent function operates on the 'getDynamicAlphaHull' function of the rangeBuilder package
# Each time it encounters an error in alpha value calculation, this new function overcomes the error 
# by adding a unit to the alpha value where the fault was generated

# Function to calculate alpha given a number of records
calculateAlpha <- function(spp_data) {
  num_records <- nrow(spp_data)
  if (num_records < 3) {
    return(NA)  # Species must have at least three occurrences to allow creation of a polygon
  } else {
    maxAlpha <- 400  # Set the maximum limit of alpha
    for (initialAlpha in 1:maxAlpha) {
      tryCatch({
        alpha_hull <- getDynamicAlphaHull(spp_data, fraction = 1,
                                          partCount=1, initialAlpha=initialAlpha, 
                                          coordHeaders=c("longitude", "latitude"),
                                          verbose=TRUE, alphaIncrement=0.5, clipToCoast="no")
        return(alpha_hull$alpha)
      }, error = function(err) {
        # Do nothing and proceed with the next attempt
      })
    }
    return("400")  # Assign '400' if the maximum limit of attempts is reached
  }
}

# Create a data frame to store the results
results_dym <- data.frame(spp = character(), alpha = numeric(), stringsAsFactors = FALSE)
# Calculate alpha values for each species
species_names <- names(spp.list.2)

for (spp_name in species_names) {
# Retrieve the data for the current species
  spp_data <- spp.list.2[[spp_name]]
# Calculate the alpha value
  alpha_val <- calculateAlpha(spp_data)
# Assign maximum alpha value in the table if the species does not meet the conditions
  alpha_val <- ifelse(alpha_val == "alphaMCH", "400", sub("alpha", "", alpha_val))
# Append the results to the results_dym data frame
  results_dym <- rbind(results_dym, data.frame(spp = spp_name, alpha = alpha_val))
}

#--------------------------------------------------------------#
####                  1.2 Dynamic alpha                     ####
#--------------------------------------------------------------#

# Load packages
library(terra)
library(sf)
library(dplyr)
library(alphahull)

# Load function that convert the alphahull in a spatial object
# Download this function by Andrew Bevan from https://stat.ethz.ch/pipermail/r-sig-geo/2012-March/014409.html
source(file = "ah2sf.R")

### Read and prepare data ###

## Note: The object that contains the species occurrence records, "spp.list.2", is found in section 1.1

# Read the reference raster layer for fitting polygons to water bodies
# we used as background a raster generated by Domisch et al. (2015), available at https://www.earthenv.org/
mask<- rast(".")
mask[!is.na(mask)]<-1

# Read the shapefile of HydroBASINS level 8 layer (Lehner & Grill, 2013), 
# available in https://www.hydrosheds.org/products/hydrobasins 
basins <- st_read(dsn = ".") %>% # Load the basin layer
  st_transform(4326) %>% # Transform the coordinate reference system (CRS)
  st_make_valid() # Ensure validity of geometries

# Generate lists and data frame to store the polygons and their associated information
range_union_list <- list() # Dynamic alpha overlapped with the sub-basins (unrestricted polygon)
alphas <- list() # Unrestricted polygons in raster format, used to count pixel number
alpha_crop <- list() # Restricted polygons to water bodies

# Table to store numerical values
areas_alpha <- data.frame(ID = integer(), species = character(), 
                          area.alpha.dym = numeric(), area.alpha.dym.rest = numeric(),
                          stringsAsFactors = FALSE)  

for (i in 1:length(spp.list.2)) {
  print(i)
  species_name <- names(spp.list.2)[i]
# Obtaining the alpha value for the current species from results_dym data.frame
  alpha <- results_dym$alpha[results_dym$spp == species_name] %>% 
    as.numeric()
# Generate the alpha polygon and convert it to an sf object
  spp_alphahull <- ahull(spp.list.2[[i]], alpha = alpha) %>% 
    ah2sf()
# Intersect the alpha polygon with the sub-basins
  range_intersection <- st_intersection(basins, spp_alphahull)
  range_hydrobasins <- basins[(range_intersection),]
  range_union <- st_union(range_hydrobasins)
# Store the range_union polygon in the list
  range_union_list[[species_name]] <- range_union
# Rasterize the polygons
  alphas[[i]] <- rasterize(st_as_sf(vect(range_union)), mask, getCover = TRUE)
  names(alphas)[i] <- species_name
# Crop the raster to the water body
  alpha_crop[[i]] <- crop(alphas[[i]], mask, mask = TRUE, snap = "near")
  names(alpha_crop)[i] <- species_name
# Assign columns to the data.frame and count the number of pixels per species
  areas_alpha <- rbind(areas_alpha, data.frame(ID = i, 
                                               species = species_name, 
                                               area.alpha.dym = freq(alphas[[i]])[, 3], 
                                               area.alpha.dym.rest = freq(alpha_crop[[i]])[, 3]))
}

## Save the data frame containing numerical values ##
write.csv(areas_alpha, file = ".")

## Save the alpha shape polygons in the path work ##
# Combine all polygons into a single sf object
range_union_combined <- do.call(rbind, range_union_list)
# Convert the polygons to an sfc object
geom <- range_union_combined[, 1]
geometry <- st_sfc(geom)
# Create a data frame with the species
species_df <- data.frame(species = rownames(range_union_combined))
# Combine geometry and data into an sf object
range_union_sf <- st_sf(geometry = geometry, data = species_df)
# Define the coordinate reference system (CRS)
crs <- "+proj=longlat +datum=WGS84 +no_defs"
# Assign the CRS to the sf object
range_union_sf <- st_set_crs(range_union_sf, crs)
# Write the combined multipolygon file
st_write(range_union_sf, dsn = ".")

#--------------------------------------------------------------#
####                   1.3 Static alpha                     ####
#--------------------------------------------------------------#

# Load packages
library(terra)
library(sf)
library(dplyr)
library(alphahull)

# Load function that convert the alphahull in a spatial object
source(file = "ah2sf.R")

## Note: The data used to generate static alpha polygons, such as the background raster and sub-basins shapefile, 
# are the same files used in section 1.2 (Dynamic alpha polygons). In this case, the objects are named mask, 
# the raster layer representing water bodies, and basins, the level 8 sub-basins layer from HydroBASINS. 
# The object that contains the species occurrence records, "spp.list.2", is found in section 1.1.

# Generate lists and data frame to store the polygons and their associated information
static.union.list <- list() # Static alpha overlapped with the sub-basins (unrestricted polygon)
static.alpha <- list() # Unrestricted polygons in raster format, used to count pixel number
static.crop <- list() # Restricted polygons to water bodies

# Table to store numerical values
static.range <- data.frame(ID = integer(), species = character(), 
                          range.static = numeric(), range.static.rest = numeric(),
                          stringsAsFactors = FALSE)  

for (i in 1:length(spp.list.2)) {
  print(i)
  species_name <- names(spp.list.2)[i]
  
  # Generate the alpha polygon and convert it to an sf object
  spp_alphahull <- ahull(spp.list.2[[i]], alpha = 6) %>% 
    ah2sf()
  # Intersect the alpha polygon with the sub-basins
  range_intersection <- st_intersection(basins, spp_alphahull)
  range_hydrobasins <- basins[(range_intersection),]
  range_union <- st_union(range_hydrobasins)
  # Store the range_union polygon in the list
  static.union.list[[species_name]] <- range_union
  # Rasterize the polygons
  static.alpha[[i]] <- rasterize(st_as_sf(vect(range_union)), mask, getCover = TRUE)
  names(static.alpha)[i] <- species_name
  # Crop the raster to the water body
  static.crop[[i]] <- crop(static.alpha[[i]], mask, mask = TRUE, snap = "near")
  names(static.crop)[i] <- species_name
  # Assign columns to the data.frame and count the number of pixels per species
  static.range <- rbind(static.range, data.frame(ID = i, 
                                               species = species_name, 
                                               range.static = freq(static.alpha[[i]])[, 3], 
                                               range.static.rest = freq(static.crop[[i]])[, 3]))
}

## Save the data frame containing numerical values ##
write.csv(static.range, file = ".")

## Save the static alpha in the path work ##
# Combine all polygons into a single sf object
static.combined <- do.call(rbind, static.union.list)
# Convert the polygons to an sfc object
geom.static <- static.combined[, 1]
geometry.static <- st_sfc(geom.static)
# Create a data frame with the species
species.df.static <- data.frame(species = rownames(static.combined))
# Combine geometry and data into an sf object
static.sf <- st_sf(geometry = geometry.static, data = species.df.static)
# Define the coordinate reference system (CRS)
crs.static <- "+proj=longlat +datum=WGS84 +no_defs"
# Assign the CRS to the sf object
static.sf <- st_set_crs(static.sf, crs.static)
# Write the combined multipolygon file
st_write(static.sf, dsn = ".")

#--------------------------------------------------------------#
####                     1.4 Convex hull                    ####
#--------------------------------------------------------------#

# Load packages
library(terra)
library(sf)
library(dplyr)

## Note: The data used to generate convex-hull polygons, such as the background raster and sub-basins shapefile, 
# are the same files used in section 1.2 (Dynamic alpha polygons). In this case, the objects are named mask, 
# the raster layer representing water bodies, and basins, the level 8 sub-basins layer from HydroBASINS. 
# The object that contains the species occurrence records, "spp.list", is found in section 1.1.

# Generate lists and data frame to store the polygons and their associated information
convex_union_list <- list() # Convex hull polygons overlapped with the sub-basins (unrestricted polygon)
convex_polygon <- list() # Unrestricted convex hull polygons in raster format
convex_crop <- list() # Convex hull restricted to water bodies

# Table to store numerical values
areas_convex <- data.frame(ID = integer(), species = character(), 
                           area.chull = numeric(), area.chull.rest = numeric(),
                           stringsAsFactors = FALSE)

for (i in 1:length(spp.list)) {
  print(i)
# Transform tables into sf objects and add attributes to spatial object
  spp.list.convex <- spp.list[[i]] %>%
    select("longitude", "latitude", "species") %>%
    hsi::clean_dup(longitude="longitude",latitude="latitude",threshold=0.041665) %>%
    st_as_sf(coords=c('longitude', 'latitude'), crs="EPSG: 4326")
# Draw convex hull
  polygon <- st_convex_hull(st_union(spp.list.convex))
# Intercept convex hull with sub-basins
  hull_intersection <- st_intersection(basins, polygon)
  hull_hydrobasins <- basins[(hull_intersection),]
  hull_union <- st_union(hull_hydrobasins)
# Store the hull_union polygon in the list
  convex_union_list[[names(spp.list)[i]]] <- hull_union
# Rasterize the polygons
  convex_polygon[[i]] <- rasterize(st_as_sf(vect(hull_union)), mask, getCover = TRUE)
  names(convex_polygon)[i] <- names(spp.list)[i]
# Crop the raster to the water body
  convex_crop[[i]] <- crop(convex_polygon[[i]], mask, mask = TRUE, snap = "near")
  names(convex_crop)[i] <- names(spp.list)[i]
# Assign columns to the data.frame and count the number of pixels per spp
  areas_convex <- rbind(areas_convex, data.frame(ID = i, 
                                                 species = names(spp.list)[i], 
                                                 area.chull = freq(convex_polygon[[i]])[, 3], 
                                                 area.chull.rest = freq(convex_crop[[i]])[, 3]))
}

## Save the data frame containing numerical values ##
write.csv(areas_convex, file = ".")

## Save the convex hull in the path work ##
# Combine all convex hull polygons into a single sf object
convex_combined <- do.call(rbind, convex_union_list)
# Convert the polygons to an sfc object
geom_convex <- convex_combined[, 1]
geometry_convex <- st_sfc(geom_convex)
# Create a data frame with the species
species_df_convex <- data.frame(species = rownames(convex_combined))
# Combine geometry and data into an sf object
convex_sf <- st_sf(geometry = geometry_convex, data = species_df_convex)
# Define the coordinate reference system (CRS)
crs_convex <- "+proj=longlat +datum=WGS84 +no_defs"
# Assign the CRS to the sf object
convex_sf <- st_set_crs(convex_sf, crs_convex)
# Write the combined multipolygon file
st_write(convex_sf, dsn = ".")

#--------------------------------------------------------------#
####                     1.5 Expert maps                    ####
#--------------------------------------------------------------#
# Load packages
library(terra)
library(sf)
library(dplyr)

# Set your working directory
setwd("~/yourpath")

### Read and prepare data ###
# The expert maps were obtained from the IUCN website, available at https://www.iucnredlist.org/
# Note: All shapefiles (expert maps) must be in the same folder
exp.maps <- list.files(path = ".", pattern = ".shp$", full.names = TRUE) # Load expert maps
species_names <- sapply(strsplit(basename(exp.maps), "\\."), `[`, 1) # Extract species names from file names

# Generate lists and data frame to store the polygons and their associated information
rast_expert_map <- list() # Unrestricted expert maps polygons in raster format
expert_crop <- list() # Expert maps restricted to water bodies

# Table to store numerical values
expert_map_areas <- data.frame(ID = integer(), species = character(), 
                               area_exp_map = numeric(), area_exp_map_rest = numeric(), 
                               stringsAsFactors = FALSE)

for (i in 1:length(exp.maps)) {
  print(i)
# Get species name
  species_name <- species_names[i]
# Rasterize the polygons
  rast_expert_map[[species_name]] <- rasterize(st_as_sf(vect(exp.maps[[i]],
                                     crs="+proj=longlat +ellps=WGS84 +no_defs")), mask, getCover=TRUE)
# Cut raster to water body
expert_crop[[species_name]] <- crop(rast_expert_map[[species_name]], mask, mask = T, snap = "near")
# Assign columns to the data frame and count the number of pixels per spp 
  expert_map_areas <- rbind(expert_map_areas, data.frame(
    ID = i, species = species_name,
    area_exp_map = freq(rast_expert_map[[species_name]])[, 3],
    area_exp_map_rest = freq(expert_crop[[species_name]])[, 3]
  ))
}

## Save the data frame containing numerical values ##
write.csv(expert_map_areas, file = ".")


### Join the tables containing the numerical values for the restricted and unrestricted ranges of the water bodies
ranges_hulls <- spp_counts %>%
  left_join(results_dym, by = c("species" = "spp")) %>%
  left_join(areas_alpha, by = c("species" = "species")) %>%
  left_join(static.range, by = c("species" = "species")) %>% 
  left_join(areas_convex, by = c("species" = "species")) %>% 
  left_join(expert_map_areas, by = c("species" = "species")) %>% 
  select(-ID.x, -ID.y, -ID.x.x, -ID.y.y) %>% 
  rename(Species=species,'Number of occurrence records'=records, 'Alpha value'=alpha,
    Convex=area.chull, Convex_network=area.chull.rest, Static=range.static, 
    Static_network=range.static.rest, Dynamic=area.alpha.dym,
    Dynamic_network=area.alpha.dym.rest, Expert=area_exp_map, Expert_network=area_exp_map_rest)
  
write.csv(ranges_hulls, file = "ranges_hulls.csv")
