# To redirect the output to a textfile, do:
# Rscript run_daisy.R > output.txt

# order
# modelpft assigns the PFTs to their position in the kk sized matrix for CLASS and CTEM:
# modelpft(1:3)  = 1,     1,     0,      ! CLASS PFT 1 NDL (NdlTr)
#               EVG    DCD
# modelpft(4:6) = 1,     1,     1,       ! CLASS PFT 2 BDL (BdlTr)
#               EVG  DCD-CLD DCD-DRY
# modelpft(7:9)= 1,     1,     0,        ! CLASS PFT 3 CROP (Crop)
#              C3      C4
# modelpft(10:12)= 1,     1,   0,   ! CLASS PFT 4 GRASS (Grass)
#             C3      C4

# vmax(1:3)   = 42.0e-06, 47.0e-06, 0.00e-06,
# vmax(4:6)   = 35.0e-06, 57.0e-06, 40.0e-06,
# vmax(7:9)   = 55.0e-06, 40.0e-06, 0.00e-06,
# vmax(10:12) = 55.0e-06, 15.0e-06, 0.00e-06,


.libPaths(c("/home/cseiler/daisy/renv", .libPaths()))
library(daisy)
rm(list = ls())

setwd('/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0')

# Copy default parameter run_parameters.nml file to current work directory
# This file will be overwritten with new values before CLASSIC is run

# Number of nodes to run meta-jobs on.
parallel <- 20 # 20
# Time for each meta-job.
metajobTime <- "48:00:00"
# The farm name.
farmName <- "CLASSIC_meta"
# Directory path for the data-assimilation-CLASSICv2.0 folder
dataAssimPath <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0"
# Archive folder for meta runs.
archiveFolder <- paste0(farmName, "_ARCHIVE")
# Make the archive folder if it does not exist already.
if (!file.exists(paste0(farmName, "_ARCHIVE"))) {system(paste0("mkdir ", farmName, "_ARCHIVE"))}
Sys.sleep(1)

# Check if there is an object.rds file in the farm directory from a previous run. If there is,
# copy it to the data-assimilation-CLASSICv2.0 folder before deleting the farm directory. Since the file
# will always be called object.rds, let ga_daisy know that an object file exists in the data-assimilation-CLASSICv2.0
# folder.
if (file.exists(paste0(farmName, "/object.rds"))) {

  # Remove a previous object.rds file from the dataAssimPath and replace it with
  # this new one.
  if (file.exists(paste0(dataAssimPath, "/object.rds"))) {
    system(paste0("rm ", dataAssimPath, "/object.rds"))
    Sys.sleep(2)
  }
  
  system(paste0("mv ", farmName, "/object.rds ", dataAssimPath))
  previousObject <- paste(dataAssimPath, "object.rds", sep = "/")
  Sys.sleep(3)
  
} else {
  previousObject <- NULL
}

# Remove folder if it already exists from a previous run.
if (file.exists(farmName)) {system(paste("rm -rf", farmName))}
Sys.sleep(1)

system('cp /home/cseiler/CLASSICv2.0/classic/configurationFiles/run_parameters.nml run_parameters.nml')

# Select parameters and obtain their parameter values
parameterFile <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/run_parameters.nml'

parameterNames <- list(
  "nresorp_coeff",
  "c2n_lmax",
  "c2n_lmin",
  "c2n_smax",
  "c2n_smin",
  "c2n_rmax",
  "c2n_rmin",
  "alpha_vcmax",
  "beta_vcmax",
  "r_bnf_s",
  "b_bnf",
  "max_frac_bnf",
  "coeff_a",
  "bnfdpth",
  "k_redcoeff_vcmax",
  "r_bnf_f",
  "c_cost_bnf_s",
  "nvol_coeff",
  "mu_nit",
  "topt_nit",
  "nitdenitdpth",
  "mu_no_nit",
  "mu_n2o_nit",
  "mu_denit",
  "topt_denit",
  "mu_no_denit",
  "mu_n2o_denit",
  "denitswthrsh",
  "ntolerance",
  "nleach_coeff",
  "coeff_p_nh4",
  "coeff_p_no3",
  "b_nh4",
  "b_no3",
  "kphalf",
  "const_c2n_humus",
  "ccost_coeff")

parameterValues <- list()

for (i in parameterNames) {
  parameters <- getParameterValues(parameterFile = parameterFile, parameterName = i)
  parameterValues[[i]] <- parameters
}

# CLASSIC has 12 PFTs, however, the default configuration only uses 9 PFTs.
# The PFTs that are not being used have parameter values that are equal to zero.
# Those values need to be excluded, otherwise the optimization will not work

# For single parameter
# non_zero_indices <- which(parameterValues != 0, arr.ind = TRUE)
# parameterValues <- parameterValues[non_zero_indices]

# Apply which() function to each element of the list
# non_zero_indices <- lapply(parameterValues, function(x) which(x != 0))

# Extract non-zero values using the obtained indices
# non_zero_values <- lapply(seq_along(parameterValues), function(i) parameterValues[[i]][non_zero_indices[[i]]])

# Print the result
# print(non_zero_values)

# non_zero_indices <- lapply(parameterValues, function(x) which(x != 0, arr.ind = TRUE))
# non_zero_values <- lapply(seq_along(parameterValues), function(i) parameterValues[[i]][non_zero_indices[[i]]])

# parameterValues <- non_zero_values

# uncertainty ranges
min.csv <- paste(parameterNames, "min.csv", sep = ".")
max.csv <- paste(parameterNames, "max.csv", sep = ".")

# Read all CSV files into a list of data frames
min.csv.list <- lapply(min.csv, read.csv)
max.csv.list <- lapply(max.csv, read.csv)

# Check if number of parameter values per parameter is consistent between min and max values:
length.min <- length(unlist(min.csv.list))
length.max <- length(unlist(max.csv.list))

if (length.max != length.min) {
  stop("Error: The lengths of the minimum and maximum parameter values are not equal. Script stopped.")
}

upperBound <- max.csv.list
lowerBound <- min.csv.list

# Convert 0 to NA

 upperBound <- lapply(upperBound, function(x) {
  if(length(x) == 12){
   x[x == 0] <- NA
  }
  return(x)
})

lowerBound <- lapply(lowerBound, function(x) {
  if(length(x) == 12){
    x[x == 0] <- NA
  }
  return(x)
})

normalization <- list()

for (p in 1:length(parameterValues)) {
  df <- parameterValues[p] # default parameter values
  ub <- upperBound[p] # upper bounds per parameter
  lb <- lowerBound[p] # lower bounds per parameter
  
  # Check whether l <= df <= u
  condition_satisfied <- mapply(function(lb, df, ub) lb <= df & df <= ub, lb, df, ub)
  
  if (all(condition_satisfied, na.rm = TRUE)) {
    print(parameterNames[[p]])
    print("low < def < up satisfied")
  } else {
    print(parameterNames[[p]])
    print("low < def < up not satisfied")
  }
  
  result <- mapply(intFun.normalize, df, ub, lb)
  result <- unlist(result)
  result <- result[!is.na(result)]
  normalization[[p]] <- result
}

parameterValueLength <- sapply(normalization, length)
normParameterValues <- unlist(normalization) 

upperBound <- unlist(upperBound)
lowerBound <- unlist(lowerBound) 

# Drop all NaNs
normParameterValues <- normParameterValues[!is.na(normParameterValues)]
upperBound <- upperBound[!is.na(upperBound)]
lowerBound <- lowerBound[!is.na(lowerBound)]

upper <- upperBound-upperBound + 1
lower <- lowerBound-lowerBound

# Omit column names
names(normParameterValues)  <- NULL
names(upperBound) <- NULL
names(lowerBound) <- NULL
names(upper) <- NULL
names(lower) <- NULL

write.table(x = normParameterValues, file = "normParameterValues")

run_classic_file <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/run_classic.sh"
# run_classic_file <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/test.sh"

dir.mod <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/daisyRun/netcdf_files"
dir.ref <- "/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63"

nc.mod01 <- file.path(dir.mod, "gpp_monthly.nc")
nc.mod02 <- file.path(dir.mod, "lai_monthly.nc")
nc.mod03 <- file.path(dir.mod, "hfls_monthly.nc")
nc.mod04 <- file.path(dir.mod, "hfss_monthly.nc")
nc.mod05 <- file.path(dir.mod, "ts_monthly.nc")
nc.mod06 <- file.path(dir.mod, "albs_monthly.nc")

nc.mod07 <- file.path(dir.mod, "burntFractionAll_monthly.nc")
nc.mod08 <- file.path(dir.mod, "cVeg_annually.nc")
nc.mod09 <- file.path(dir.mod, "cSoil_annually.nc")
nc.mod10 <- file.path(dir.mod, "bnf_tot_annually.nc")
nc.mod11 <- file.path(dir.mod, "c2n_wp_annually.nc")
nc.mod12 <- file.path(dir.mod, "c2n_humus_annually.nc")

# GPP
nc.ref01a <- file.path(dir.ref, "gpp_GOSIF_128x64.nc")
nc.ref01b <- file.path(dir.ref, "GPP.ensembleMedian.CRUNCEPv6.monthly_128x64.nc")

# LAI
nc.ref02a <- file.path(dir.ref, "lai_AVHRR_128x64.nc")
nc.ref02b <- file.path(dir.ref, "lai_copernicus_128x64.nc")
nc.ref02c <- file.path(dir.ref, "lai_MODIS_128x64_LiboWang.nc")

# HFLS
nc.ref03a <- file.path(dir.ref, "CLASS_v1.1_hfls.nc")
nc.ref03b <- file.path(dir.ref, "LE.RS_METEO.EBC-ALL.MLM-ALL.METEO-ALL.720_360.monthly_128x64.nc")

# HFSS
nc.ref04a <- file.path(dir.ref, "CLASS_v1.1_hfss.nc")
nc.ref04b <- file.path(dir.ref, "H.RS_METEO.EBC-ALL.MLM-ALL.METEO-ALL.720_360.monthly_128x64.nc")

# LST
nc.ref05a <- file.path(dir.ref, "MOYD11C3_128x64.nc")

# ALBS
nc.ref06a <- file.path(dir.ref, "albs_CERES_128x64.nc")
nc.ref06b <- file.path(dir.ref, "albs_GEWEXSRB_128x64.nc")
nc.ref06c <- file.path(dir.ref, "albedo_MODIS_128x64.nc")

# burnt (fraction)
nc.ref07a <- file.path(dir.ref, "burntFractionAll_GFED4S_128x64.nc")
nc.ref07b <- file.path(dir.ref, "fractionBurnt_ESACCI-L4_FIRE-BA-MODIS-fv5.1.nc")

# cVeg (kgC m$^{-2}$)
nc.ref08 <- file.path(dir.ref, "Huang2021_CVEG_128x64.nc")

# csoil (kgC m$^{-2}$)
nc.ref09 <- file.path(dir.ref, "soc_OCSTHA_M_depth_0-100cm_128x64.nc")

# BNF 
# nc.ref10 <- file.path(dir.ref, "bnfmap_ilamb_128x64.nc") # (gN m$^{-2}$ y$^{-1}$)
nc.ref10 <- file.path(dir.ref, "BNF_USGS_T63.nc") # (kgN ha$^{-1}$ y$^{-1}$)

# C:N veg (kgC kgN$^{-1}$)
# nc.ref11 <- file.path(dir.ref, "CNveg_128x64.nc")
nc.ref11 <- file.path(dir.ref, "CNveg_T63.nc")

# C:N soil (kgC kgN$^{-1}$)
nc.ref12 <- file.path(dir.ref, "CNsoil_128x64.nc")


mod.list <- list(
  nc.mod01, nc.mod01,
  nc.mod02, nc.mod02, nc.mod02,
  nc.mod03, nc.mod03,
  nc.mod04, nc.mod04,
  nc.mod05, 
  nc.mod06, nc.mod06, nc.mod06,
  nc.mod07, nc.mod07,
  nc.mod08,
  nc.mod09,
  nc.mod10,
  nc.mod11,
  nc.mod12)

ref.list <- list(
  nc.ref01a, nc.ref01b,
  nc.ref02a, nc.ref02b, nc.ref02c,
  nc.ref03a, nc.ref03b,
  nc.ref04a, nc.ref04b,
  nc.ref05a,
  nc.ref06a, nc.ref06b, nc.ref06c,
  nc.ref07a, nc.ref07b,
  nc.ref08,
  nc.ref09,
  nc.ref10,
  nc.ref11,
  nc.ref12)

# Unit conversion factors for reference data

ref.unit.conv.list <- list(
  1/86400000, 1/86400000,
  1,1,1,
  1, (10^6/86400),
  1, (10^6/86400),
  1,
  1,1,1,
  1,1,
  1,
  1,
  1/(10000*365.25*86400),
  1,
  1)

ref.id.list <- list("GPP-GOSIF", "GPP-FLUXCOM",
                    "LAI-AVHRR", "LAI-Copernicus", "LAI-MODIS",
                    "HFLS-CLASSr", "HFLS-FLUXCOM",
                    "HFSS-CLASSr", "HFSS-FLUXCOM",
                    "TS-MODIS",
                    "ALBS-CERES", "ALBS-GEWEXSRB", "ALBS-MODIS",
                    "BURNT-GFED4S", "BURNT-ESACCI",
                    "CVEG-H2021",
                    "CSOIL-SG250m",
                    "BNF-USGS",
                    "CNVeg-KG2023",
                    "CNSoil-KG2023")

# modelOutputFolder is the folder for the model while the model is running the individual.
modelOutputFolder <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/daisyRun"

# By default, keepOutput is FALSE and simFolder is NULL.
# keepOutput is a logical argument determining whether output is kept.
keepOutput <- FALSE
# simFolder is the folder where finished model output is moved to.
simFolder <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations"

library("GA")

source("ga_daisy.R")

setwd("/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0")

# Different options for selection, crossover, and mutation:

selection <- c(
  "gareal_lrSelection",
  "gareal_nlrSelection", 
  "gareal_rwSelection",
  "gareal_tourSelection", 
  "gareal_lsSelection", 
  "gareal_sigmaSelection")

crossover <- c(
  "gareal_spCrossover",
  "gareal_waCrossover",
  "gareal_laCrossover",
  "gareal_blxCrossover",
  "gareal_laplaceCrossover")

mutation <- c(
  "gareal_raMutation",
  "gareal_nraMutation",
  "gareal_rsMutation",
  "gareal_powMutation")


# options S1-6, C1-5, M1-4

# Find the importance of selection, crossover, and mutation
S1C1M1 <- c(selection[1], crossover[1], mutation[1])
S2C1M1 <- c(selection[2], crossover[1], mutation[1])
S3C1M1 <- c(selection[3], crossover[1], mutation[1])
S4C1M1 <- c(selection[4], crossover[1], mutation[1])
S5C1M1 <- c(selection[5], crossover[1], mutation[1])
S6C1M1 <- c(selection[6], crossover[1], mutation[1])
S1C2M1 <- c(selection[1], crossover[2], mutation[1])
S1C3M1 <- c(selection[1], crossover[3], mutation[1])
S1C4M1 <- c(selection[1], crossover[4], mutation[1])
S1C5M1 <- c(selection[1], crossover[5], mutation[1])
S1C1M2 <- c(selection[1], crossover[1], mutation[2])
S1C1M3 <- c(selection[1], crossover[1], mutation[3])
S1C1M4 <- c(selection[1], crossover[1], mutation[4])

# Choosing M = 2 and setting C = 1, find the best method for S
S2C1M2 <- c(selection[2], crossover[1], mutation[2])
S3C1M2 <- c(selection[3], crossover[1], mutation[2])
S4C1M2 <- c(selection[4], crossover[1], mutation[2])
S5C1M2 <- c(selection[5], crossover[1], mutation[2])
S6C1M2 <- c(selection[6], crossover[1], mutation[2])

# Choosing M = 2 and S = 5, find the best method for C
S5C2M2 <- c(selection[5], crossover[2], mutation[2])
S5C3M2 <- c(selection[5], crossover[3], mutation[2])
S5C4M2 <- c(selection[5], crossover[4], mutation[2])
S5C5M2 <- c(selection[5], crossover[5], mutation[2])

# selCroMut <- S5C1M2
selCroMut <- S1C1M1
result <- ga_daisy(
    type = "real-valued",
    fitness = cost.fun,
    lowerBound = lowerBound,
    upperBound = upperBound,
    parameterValueLength = parameterValueLength,
    parameterNames = parameterNames, 
    parameterFile = parameterFile,
    mod.list = mod.list,
    ref.list = ref.list,
    ref.id.list = ref.id.list,
    ref.unit.conv.list = ref.unit.conv.list,
    run_classic_file = run_classic_file,
    modelOutputFolder = modelOutputFolder,
    lower = lower,
    upper = upper,
    
  #  selection = selCroMut[1],
  #  crossover = selCroMut[2],
  #  mutation = selCroMut[3],
    
    popSize = 100, # 100
    elitism = 4, # 4,
    maxiter = 25, # 25
    run = 25, # 25
    maxFitness = 1,
    parallel = parallel,
    jobTime = metajobTime,
    farmName = farmName,
    dataAssimPath = dataAssimPath,
    previousObject = previousObject,
    suggestions = normParameterValues,
    keepBest = TRUE,
    keepOutput = keepOutput,
    archiveFolder = archiveFolder,
    finalOutputFolder = simFolder,
    seed = 1)
saveRDS(result, "result.rds")

