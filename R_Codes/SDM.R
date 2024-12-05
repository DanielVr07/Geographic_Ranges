### This code replicates the construction and estimation of the Species Distribution Models
# (SDMs) for fish of the Curimatidae family from the paper titled: 
# "Comparing methods for estimating geographic ranges in freshwater fishes"
# by Valencia-Rodríguez et al. Submitted to Freshwater Biology

#-------------------------------------------------------#
####              2.1 Occurrence records             ####
#-------------------------------------------------------#

# Load packages
library(readr)
library(dplyr)

# Set your working directory
setwd(~/yourpath)

## Note: To develop the SDMs, we must separate the occurrence records for each species. 
# Additionally, to develop the SDMs, we use the protocol proposed by Cobos et al. (2019) 
# https://doi.org/10.7717/peerj.6281 from the kuenm package. This protocol requires that 
# the data for each species be divided into three files: (i) sp_train, (ii) sp_test, and (iii) sp_joint.

### Read and prepare data ###
# Read the clean species occurrence database. 
# This database only requires three columns: (i) species name, (ii) longitude, and (iii) latitude
data_ocurr <- read_csv ("curimatidae_occ_2023.csv")

# Create three .csv files with the occurrence records
# sp_train: contains the training points for the models (80% of the total)
# sp_test: contains the testing points for the candidate models (20% of the total)
# sp_joint: contains all the complete records

# Create folders and files for each species
for (i in unique(data_ocurr$species)) {
  # Create a folder in the directory for each species
  dir.create(paste0("models/", i, "/"), recursive = TRUE)
  
  # Filter the occurrence data for the current species
  data_joint <- data_ocurr %>%
    filter(species == i) %>%
    select(species, longitude, latitude) %>%
    rename(specie = species)
  
  # Split the data into training and testing sets
  samples <- sample(2, nrow(data_joint), replace = TRUE, prob = c(0.8, 0.2))
  data_train <- data_joint[samples == 1, ]
  data_test <- data_joint[samples == 2, ]
  
  # Save the .csv files
  write.csv(data_train, paste0("models/", i, "/", "Sp_train.csv"), row.names = FALSE)
  write.csv(data_test, paste0("models/", i, "/", "Sp_test.csv"), row.names = FALSE)
  write.csv(data_joint, paste0("models/", i, "/", "Sp_joint.csv"), row.names = FALSE)
}

## Note: Please note that the data separation is random, so each data split may generate different 
# results in the models, as the training and evaluation data are not necessarily the same. 
# However, in the "models" folder, we provide the data used to build the models for our study.

#-------------------------------------------------------#
####              2.2 Variable selection             ####
#-------------------------------------------------------#

# Load packages
library(terra)
library(readr)
library(usdm)
library(dplyr)
library(stringr)
library(tibble)

# Set your working directory
setwd("~/yourpath")

# This section replicates the process of selecting variables (climatic and environmental). 
# Since the process is specific to each species, the code only exemplifies the procedure for one species.
# However, variable selection for each species is available in the predictors file.


### Read and prepare data ###
# Load and stack the climatic and environmental variables.

# Note: All raster files corresponding to a category (climatic or environmental) 
# should be placed in the same folder

climatic <- rast(list.files(".", pattern='.tif', full.names=TRUE))
names(climatic) <- list.files(".", pattern='.tif', full.names=FALSE)
environmental <- rast(list.files(".", pattern='.tif', full.names=TRUE))
names(environmental) <- list.files(".", pattern='.tif', full.names=FALSE)

# Load species occurrence records
# The file 'Sp_join' contains the occurrence records for each species, saved in the folder 'models,' 
# obtained in section in section 2.1.
occ <- read_csv("models/Curimata aspera/Sp_joint.csv")
occ_cur <- occ[, c(2, 3)] # Select the columns "longitude and latitude"

## Check for collinearity among climatic variables
# Extract the values of the climatic variables for each record
pred_clim <- extract(climatic, occ_cur) 
# Test for collinearity among the predictors
no_cor_clim <- vifstep(pred_clim[, 2:16], th=5)
# Extract variables with low collinearity
clim_vars <- no_cor_clim@results[["Variables"]]

## Check for collinearity among environmental variables
# Extract the values of the environmental variables for each record
pred_amb <- extract(environmental, occ_cur)
# Test for collinearity among the predictors
no_cor_amb <- vifstep(pred_amb[, 2:33], th = 5)
# Extract variables with low collinearity
amb_vars <- no_cor_amb@results[["Variables"]]

# Combine the variables lists
result <- occ %>%
  pull(specie) %>%                 # Extract the species name
  unique() %>%                     # Ensure the species name is unique
  {tibble(Species = .,             # Create a tibble with the species name
          Predictors = c(clim_vars, amb_vars) %>%  # Combine the climatic and environmental variables
            str_remove_all("\\.tif$") %>%          # Remove the ".tif" extension from the variable names
            paste(collapse = ","))}                # Combine the variable names into a single comma-separated string


write.csv(result, dsn = ".")

#-------------------------------------------------------#
####           2.3 Processing of predictors          ####
#-------------------------------------------------------#

# Load packages
library(readr)
library(terra)
library(sf)

# Set your working directory
setwd("~/yourpath")


## The selected variables for each species are cropped based on their accessible area. 
# Additionally, each variable is exported to the species' subfolder within the "models" folder.

### Read and prepare data ###

# Load shapefiles with the accessible areas of each species 
# Note: The name of the accessible area file for each species should match the species name. 
# Also, all shapes must be in the same folder.
spp.shape <- list.files(".", full.names = TRUE, pattern=".shp$") 

# Load selected variables for each species, available in the predictors file
spp.set <- read_csv("Predictors.csv") 
# List all files within each species' folder in the "models/" directory
spp.carp.list <- list.files("models/", full.names = TRUE)

for (i in 1:length(spp.carp.list)) {
  print(spp.carp.list[i])
  
  # Create the M_variables folder in the folder of each species
  dir.create(paste(spp.carp.list[i], "/M_variables", sep = ""))
  # Identify variables for each species
  species <- gsub("_", " ", substr(basename(spp.carp.list[i]), 1, 
                                   (nchar(basename(spp.carp.list[i])))))
  
  # Filter spp.set for the species
  subs.spp.actual <- spp.set[spp.set$Species == species, ]
  # Get the predictors for the species
  predictors.species <- unlist(strsplit(subs.spp.actual$Predictors, ","))
  # Filter the list of predictors according to those selected for the species
  selected.predictors <- predictors.species[predictors.species %in% predictors.species]
  
  # Load only the selected predictors and stack them
  raster_stack <- rast(paste0("Variables/", selected.predictors, ".tif"))
  names(raster_stack) <- selected.predictors
  
  # Load the accessible area "M" for the species and assign CRS
  shp.spp <- vect(spp.shape[i])
  crs(shp.spp) <- raster_stack[[1]]
  
  # Create a directory for the stack of variables
  stack_dir <- file.path(spp.carp.list[i], "M_variables", "Set_01")
  dir.create(stack_dir, showWarnings = FALSE)
  
  # Crop the stack to the extent of the accessible area
  subs.stack.crop <- crop(raster_stack, shp.spp, mask = TRUE, touches = TRUE)
  # Save variables
  writeRaster(subs.stack.crop, paste0(spp.carp.list[i], "/M_variables/Set_01/",
                names(subs.stack.crop), ".asc"), overwrite = TRUE, NAflag = -9999)
}

#-------------------------------------------------------#
####      2.4 Species Distribution Models (SDM)      ####
#-------------------------------------------------------#

# Load packages
library(kuenm)

# Set your working directory
setwd("~/yourpath")

# Note: The entire modeling protocol of the kuenm package is available at 
# https://github.com/marlonecobos/kuenm. 
# In this section, we only share the parameters used to model the Curimatidae fish.

### Read data ###
# Load the directory containing the species to be modeled, in this case, 
# Curimatidae fish located in the "models/" folder
maxent.models <- paste("models/", list.files("models/"), sep="")
 
for(i in 1:length(maxent.models)){
  setwd(maxent.models[i]) 
  # The next chunk of code is for preparing the arguments for using the function following the modularity principle. These variables can be changed according to each case.
  occ_joint <- "Sp_joint.csv"
  occ_tra <- "Sp_train.csv"
  M_var_dir <- "M_variables"
  batch_cal <- "Candidate_models"
  out_dir <- "Candidate_Models"
  reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
  f_clas <- c("l", "lq", "lqp")
  background <- 10000
  maxent_path <- "C:/maxent"
  wait <- FALSE
  run <- TRUE
  
  # The following is the code for using the function.
  kuenm_cal(occ.joint=occ_joint, occ.tra=occ_tra, M.var.dir=M_var_dir, batch=batch_cal, out.dir=out_dir, 
            reg.mult=reg_mult, f.clas=f_clas, maxent.path=maxent_path, wait=wait, run=run)
  
  ## Evaluation and selection of best models
  occ_test <- "Sp_test.csv"
  out_eval <- "Calibration_results"
  threshold <- 5
  rand_percent <- 50
  iterations <- 500
  kept <- TRUE
  selection <- "OR_AICc"
  paral_proc <- FALSE 
  
  # This code also allows evaluating candidate models that were created previously, 
  # selecting those with best performance based on the three criteria.
  cal_eval <- kuenm_ceval(path=out_dir, occ.joint=occ_joint, occ.tra=occ_tra, occ.test=occ_test, 
                          batch=batch_cal, out.eval=out_eval, threshold=threshold, rand.percent=rand_percent,
                          iterations=iterations, kept=kept, selection=selection, parallel.proc=paral_proc)
  
  ### Final model creation
  # For preparing the arguments for this function use the following chunk of code.
  batch_fin <- "Final_models"
  mod_dir <- "Final_Models"
  rep_n <- 10
  rep_type <- "Bootstrap"
  jackknife <- FALSE
  out_format <- "logistic"
  project <- FALSE
  G_var_dir <- "G_variables"
  ext_type <- "all"
  write_mess <- FALSE
  write_clamp <- FALSE
  wait1 <- FALSE
  run1 <- TRUE
  args <- NULL
  
  # The kuenm_mod function has the following syntax:
  kuenm_mod(occ.joint=occ_joint, M.var.dir=M_var_dir, out.eval=out_eval, batch=batch_fin, rep.n=rep_n, 
            rep.type=rep_type, jackknife=jackknife, out.dir=mod_dir, out.format=out_format, project=project, 
            G.var.dir=G_var_dir, ext.type=ext_type, write.mess=write_mess, write.clamp=write_clamp, 
            maxent.path=maxent_path, args=args, wait=wait1, run=run1)
}

#-------------------------------------------------------#
####          2.5 Extraction of SDM Results          ####
#-------------------------------------------------------#

## This section extracts the results of the SDM evaluations for each species and stores them in a data frame.
# The data are extracted from the tables containing the model evaluations performed with the kuenm package. 

# Set your working directory
setwd("~/yourpath")

### Read and prepare data ###

# Note: The models generated for each species were manually reviewed,
# ensuring that the best statistical model selected was consistent with the biological 
# knowledge obtained for the species.

# Table with the species names and the selected model, available in the best.model file
tab.mods <- read.csv("best.model.csv")

# Directory where the models are located
models.dir <- list.files("models/", full.names = TRUE)
spp.mod.list <- list.files("models/", full.names = FALSE) # Recognize species names  

#  Construct data.frame to store results of 10 percentile training presence logistic threshold
table.var.bin <- data.frame(spp=list.files("models/", full.names=FALSE),
      '10ptp.median'=NA, '10ptp.mean'=NA, '10ptp.max'=NA, '10ptp.rep.max'=NA,
      '10ptp.min'=NA, '10ptp.rep.min'=NA, Mean_AUC_ratio=NA, pval_pROC=NA,
      Omiss.rat=NA, AICc=NA, delta_AICc=NA, W_AICc=NA, num_par=NA)

## Extract values for each species and some metrics from the model evaluation

for (i in 1:length(tab.mods[,1])){
  # Open folder with final models
  carp.fm <- list.files(paste0(models.dir[i], "/Final_Models"), full.names=FALSE)
  
  # Apply an if condition to handle the discrepancy between the folder names and the names in the dataset
  if (tab.mods[i,1] == substr(spp.mod.list,1,nchar(spp.mod.list))[i]) { 
    ## Calculate the median
    # Open Maxent results
    max.res <- read.csv(paste0(paste0(models.dir[i], "/Final_Models/", tab.mods[i, 2], "/maxentResults.csv")))
    
    # Calculate statistics for the 10 percentile training presence logistic threshold
    ptp.lh<-max.res$X10.percentile.training.presence.Logistic.threshold
    table.var.bin[i,2]<-median(ptp.lh[1:10]) 
    table.var.bin[i,3]<-mean(ptp.lh[1:10]) 
    table.var.bin[i,4]<-max(ptp.lh[1:10]) 
    table.var.bin[i,5]<-max.res[which(max.res$X10.percentile.training.presence.Logistic.threshold==max(ptp.lh[1:10])),1]
    table.var.bin[i,6]<-min(ptp.lh[1:10])
    table.var.bin[i,7]<-max.res[which(max.res$X10.percentile.training.presence.Logistic.threshold==min(ptp.lh[1:10])),1]
    
    # Calculate calibration results
    cal.res <- 
    read.csv(paste0(paste0(models.dir[i], "/Calibration_results/calibration_results.csv")))
    table.var.bin[i,8]<-cal.res[cal.res[,1] %in% tab.mods[i,2],][,2]
    table.var.bin[i,9]<-cal.res[cal.res[,1] %in% tab.mods[i,2],][,3]
    table.var.bin[i,10]<-cal.res[cal.res[,1] %in% tab.mods[i,2],][,4]
    table.var.bin[i,11]<-cal.res[cal.res[,1] %in% tab.mods[i,2],][,5]
    table.var.bin[i,12]<-cal.res[cal.res[,1] %in% tab.mods[i,2],][,6]
    table.var.bin[i,13]<-cal.res[cal.res[,1] %in% tab.mods[i,2],][,7]
    table.var.bin[i,14]<-cal.res[cal.res[,1] %in% tab.mods[i,2],][,8]
  }
}

# Save table with results
write.csv(tabla.var.bin, "Eval_SDM.csv", row.names = FALSE)

#-------------------------------------------------------#
####        2.6 Reclassification of the models       ####
#-------------------------------------------------------#

# Load packages
library(dplyr)
library(readr)
library(terra)

# Set your working directory
setwd("~/yourpath")

### Read and prepare data ###

# Generate a list of directories with the models to be reclassified
Mods.dir <- list.files("models/", full.names = TRUE)

# Generate paths for the folders of the final models to be reclassified
pathway <- c()
for (i in 1:length(Mods.dir)) {
  complement <- list.files(paste0("", Mods.dir[i], "/Final_Models"))[1]
  pathway[i] <- paste0("", Mods.dir[i], "/Final_Models/", complement)
}

# Generate paths for the raster to be reclassified (in this case median.asc)
median <- c()
for (i in 1:length(pathway)) {
  complement.2 <- list.files(pathway[i])
  median[i] <- paste0(pathway[i], "/", complement.2[which(grepl('median', complement.2))])
}

# Load the file with the species names and the cut-off value, in this study the 10 percentile training presence
## Note: The "Eval_SDM" file was generated in section 2.5
threshold_data <- read_csv("~/yourpath/Eval_SDM.csv") %>%
  mutate(spp = gsub(" ", "_", spp)) %>%
  mutate(spp = paste0(spp, "_median")) %>% 
  select(spp, X10ptp.median)

# Obtain a unique list of species names
species <- unique(threshold_data$spp) 
# Create a folder to store the binarized models
dir.create(paste0("Binary/")) 
# Create a table for the numerical values of the range size
sdm.range.size <- data.frame(Species=character(), Range.size=numeric(), stringsAsFactors=FALSE) 

for (i in 1:length(species)) {
  print(species[i])
  
  # Load the raster
  raster_original <- rast(median[i])
  umbral <- threshold_data$X10ptp.median[threshold_data$spp == species[i]]
  
  # Reclassify raster using the cut-off value
  raster_original[raster_original >= umbral] <- 1
  raster_original[raster_original < umbral] <- 0
  
  # Calculate the species range size
  sdm.range.size[i, "Species"] <- substr(species[i], 1, nchar(species[i])-7)
  sdm.range.size[i, "Range.size"] <- freq(raster_original)[2, 3]
  
  # Save binarized raster
  writeRaster(raster_original, paste0("Binary/", 
              paste0(substr(species[i], 1, nchar(species[i])-7), ".tif")))
}

# Replace underscores with spaces in the "Species" column
sdm.range.size$Species <- gsub("_", " ", sdm.range.size$Species)

# Save the table with range size values of the SDMs 
write.csv(sdm.range.size, "range_size_SDM.csv", row.names = FALSE)


### Load the table with the range sizes estimated using polygonal methods and 
# merge it with the range sizes estimated using SDMs

# Load the table with the range sizes from the polygonal methods
# Note: This table is the final result from section 1.4
range.polygons <- read_csv("~/yourpath/ranges_hulls.csv")

# Merge the tables containing the range size values from the polygonal methods and SDMs
ranges_size <- range.polygons %>% 
  left_join(sdm.range.size, by = c("Species" = "Species")) %>% 
  rename(SDM = Range.size)

# Save the table with the range size values estimated for the five methods
write.csv(ranges_size, "ranges_size.csv", row.names = FALSE)
