rm(list = ls())
setwd("/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0")

# The PFT values in the file below are sorted by row rather than column

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



PFTs <- c("NdlEvgTr", "NdlDcdTr", NA, "BdlEvgTr", "BdlDCoTr", "BdlDDrTr", "CropC3", "CropC4", NA, "GrassC3", "GrassC4", NA)

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

# drop some columns and sort alphabetically
data <- subset(data, select = -c(ID,  X1stValue))
data <- data[order(data$Parameter), ]

Parameter <- paste(data$Parameter, data$units, sep = " ")
Parameter <- rep(Parameter, each = 3)

categories <- c("def", "min", "max")
categories <- rep(categories, nrow(data))

Parameter <- paste(Parameter, categories, sep = " ")

Default <- data$default
Minimum <- data$Minimum
Maximum <- data$Maximum

data <- data.frame(Default, Minimum, Maximum)
data <- t(data)

Values <- matrix(data, ncol = 1)
data <- data.frame(Parameter, Values)


# Print LeTeX table
data <- xtable::xtable(data)
xtable::print.xtable(data)
print(data, file="parameterUncertaintyRanges.tex")




