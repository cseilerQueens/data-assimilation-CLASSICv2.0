rm(list = ls())
setwd("/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0")

# The PFT values in the file below are sorted by row rather than column
# My scripts expect parameter values sorted by column rather than row
# I therefore need to correct the order of PFT-specific parameter values
# This conversion is as follows:

# original position -> new position
# 1 -> 1
# 2 -> 5
# 3 -> 9
# 4 -> 2
# 5 -> 6
# 6 -> 10
# 7 -> 3
# 8 -> 7
# 9 -> 11
# 10 -> 4
# 11 -> 8
# 12 -> 12

# Same thing (original position -> new position)
#	1	->	1
#	4	->	2
#	7	->	3
#	10	->	4
#	2	->	5
#	5	->	6
#	8	->	7
#	11	->	8
#	3	->	9
#	6	->	10
#	9	->	11
#	12	->	12

# Define your desired order
new_order <- c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)

# Get uncertainty ranges from Raj Deepak S.N.
data <- read.csv("parameterUncertaintyRanges.csv", sep = ";")

# Subset the parameters you need

parameterNames <- list("vmax", "kn", "thrprcnt", "alpha_phtsyn", "lfespany", "grescoef", 
                       "avertmas", "mxrtdpth", "sn", "XLEAF", "maxage", "ZOLNG", 
                       "coldthrs", "alpha", "gamma_w", "beta2", "kappa", "minslai", 
                       "coldlmt", "roothrsh", "abar", "vpd0", "albnir", "TCSAND", 
                       "TCCLAY", "ZOLNS", "minlvfr", "omega")

# PFTs.max <- c("NdlEvgTr.max","BdlEvgTr.max","CropC3.max","GrassC3.max","NdlDcdTr.max","BdlDCoTr.max","CropC4.max","GrassC4.max","BdlDDrTr.max")
# PFTs.min <- c("NdlEvgTr.min","BdlEvgTr.min","CropC3.min","GrassC3.min","NdlDcdTr.min","BdlDCoTr.min","CropC4.min","GrassC4.min","BdlDDrTr.min")

data <- data[data$Parameter %in% parameterNames, ]

# Get the min and max values and write them to csv files, two for each parameter

for (i in parameterNames){
  
  # maximum
  parval <- data[data$Parameter == i, ]$Maximum
  parval <- as.numeric(strsplit(parval, ",")[[1]])
  
  if (length(parval) == 12){
    parval <- parval[new_order] # correct order of values
 #  parval <- parval[c(1,2,3,4,5,6,7,8,10)] # omit the three PFTs that are not used
  }
    parval <- data.frame(t(parval))

  fname <- paste(i, "max.csv", sep = ".")
  write.csv(x = parval, file = fname, row.names = FALSE)
  length.max <- length(parval)
  rm(parval)

  # minimum
  parval <- data[data$Parameter == i, ]$Minimum
  parval <- as.numeric(strsplit(parval, ",")[[1]])
  
  if (length(parval) == 12){
    parval <- parval[new_order] # correct order of values
#   parval <- parval[c(1,2,3,4,5,6,7,8,10)] # omit the three PFTs that are not used
  }
    parval <- data.frame(t(parval))

  fname <- paste(i, "min.csv", sep = ".")
  write.csv(x = parval, file = fname, row.names = FALSE)
  length.min <- length(parval)
  rm(parval)
  
  if (length.max != length.min) {
    stop("Error: The lengths of the minimum and maximum parameter values are not equal. Script stopped.")
  }
  
}

