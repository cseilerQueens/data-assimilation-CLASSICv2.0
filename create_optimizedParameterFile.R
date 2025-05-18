.libPaths("/home/cseiler/daisy/renv")
library(daisy)
library(GA)

rm(list = ls())

setwd('/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0')

# (1) Get the parameterValueLength

# Copy default parameter run_parameters.nml file to current work directory
# This file will be overwritten with new values before CLASSIC is run

system('cp /home/cseiler/CLASSICv2.0/classic/configurationFiles/default_run_parameters.nml run_parameters.nml')

# Select parameters and obtain their parameter values
parameterFile <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/run_parameters.nml'

parameterNames <- list("vmax", "kn", "thrprcnt", "alpha_phtsyn", "lfespany", "grescoef", 
                       "avertmas", "mxrtdpth", "sn", "XLEAF", "maxage", "ZOLNG", 
                       "coldthrs", "alpha", "gamma_w", "beta2", "kappa", "minslai", 
                       "coldlmt", "roothrsh", "abar", "vpd0", "albnir", "TCSAND", 
                       "TCCLAY", "ZOLNS", "minlvfr", "omega")

parameterValues <- list()

for (i in parameterNames) {
  parameters <- getParameterValues(parameterFile = parameterFile, parameterName = i)
  parameterValues[[i]] <- parameters
}

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

# This is the parameter value length where zeros that are place holders for PFTs are dropped
parameterValueLength <- sapply(normalization, length)
upperBound <- unlist(upperBound)
lowerBound <- unlist(lowerBound) 

# Drop all NaNs
# normParameterValues <- normParameterValues[!is.na(normParameterValues)]
upperBound <- upperBound[!is.na(upperBound)]
lowerBound <- lowerBound[!is.na(lowerBound)]

# Omit column names
names(upperBound) <- NULL
names(lowerBound) <- NULL

# (2) Convert the optimized normalized parameter values to actual parameter values and write them to a parameter file

object <- readRDS("analysis/ALBS-GPP-HFLS-HFSS-LAI-TS-PAR28-GC160-P100-I1701/objectFinal.rds")
normParameterValues <- object@solution
# object <- readRDS("../data-assimilation-CLASSICv2.0/object.rds")
# normParameterValues <- object@bestSol[[14]]

n <- length(parameterNames)

# Convert vector to list, where each element represents one parameter
breakpoints <- rep(seq_along(parameterValueLength), times = parameterValueLength)
normParameterValuesList <- split(x = normParameterValues, f = breakpoints)

upperBoundList <- split(x = upperBound, f = breakpoints)
lowerBoundList <- split(x = lowerBound, f = breakpoints)

# The loop edits the parameter file for each parameter

daisyOutput <- numeric()

for (i in 1:n) {
  
  # Get the parameter names
  parameterName <- parameterNames[[i]]
  
  # un-normalize parameter values
  parameterValues <- intFun.unnormalize(normParameterValuesList[[i]], upperBoundList[[i]], lowerBoundList[[i]])
  
  daisyOutput <- c(daisyOutput, parameterValues)
  
  # Get the default parameter values
  defaultParameterValues <- getParameterValues(parameterFile, parameterName)
  
  non_zero_indices <- which(defaultParameterValues != 0, arr.ind = TRUE)
  newValues <- defaultParameterValues
  
  # Replace prior with new parameter values, skipping all instances where
  # the parameter values equal zero in the parameter file
  
  newValues[non_zero_indices] <- defaultParameterValues[non_zero_indices] -
    defaultParameterValues[non_zero_indices] + parameterValues[!is.na(parameterValues)]
  
  # Overwrite the parameter file with the new parameter values
  lines <- editParameterFile(parameterFile, parameterName, parameterValues = newValues)
  writeLines(lines, parameterFile)
  
}
