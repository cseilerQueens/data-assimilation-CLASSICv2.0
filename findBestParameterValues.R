.libPaths("/home/cseiler/daisy/renv")
library(daisy)

rm(list = ls())
setwd("/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0")

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

# default values

# vmax(1:3)   = 42.0e-06, 47.0e-06, 0.00e-06,
# vmax(4:6)   = 35.0e-06, 57.0e-06, 40.0e-06,
# vmax(7:9)   = 55.0e-06, 40.0e-06, 0.00e-06,
# vmax(10:12) = 55.0e-06, 15.0e-06, 0.00e-06,

# kn(1:3)   = 0.50, 0.50, 0.00,
# kn(4:6)   = 0.50, 0.50, 0.50,
# kn(7:9)   = 0.40, 0.48, 0.00,
# kn(10:12) = 0.46, 0.44, 0.00,


# vmax: 4.2e-05 3.5e-05 5.5e-05 5.5e-05 4.7e-05 5.7e-05 4.0e-05 1.5e-05 4.0e-05
# kn: 0.50 0.50 0.40 0.46 0.50 0.50 0.48 0.44 0.50
# "NDL-EVG", "BDL-EVG", "CROP-C3", "GRASS-C3", "NDL-DCD", "BDL-DCD-CLD", "CROP-C4", "GRASS-C4", "BDL-DCD-DRY"

# default values
vmax <- c(4.2e-05, 3.5e-05, 5.5e-05, 5.5e-05, 4.7e-05, 5.7e-05, 4.0e-05, 1.5e-05, 4.0e-05)
kn <- c(0.50, 0.50, 0.40, 0.46, 0.50, 0.50, 0.48, 0.44, 0.50)
default.list <- list(vmax, kn)

# uncertainty ranges
vmax.min <- read.csv("vmax.min.csv", header = TRUE)  * 10^(-6)
vmax.max <- read.csv("vmax.max.csv", header = TRUE)  * 10^(-6)
kn.min <- read.csv("kn.min.csv", header = TRUE)
kn.max <- read.csv("kn.max.csv", header = TRUE)
upperBound.list <- list(vmax.max, kn.max)
lowerBound.list <- list(vmax.min, kn.min)

parameterName.list <- list("vmax", "kn")

ylab.vmax <- latex2exp::TeX("Maximum Carboxylation Rate ($\\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$)")
ylab.kn <- "Canopy extinction coefficient (-)"
ylab.list <- list(ylab.vmax, ylab.kn)

factor.list <- list(10^6, 1) # unit conversion factor

legLoc.list <- list("bottomleft", "topleft")

PFTs <- c("NDL-EVG", "BDL-EVG", "CROP-C3", "GRASS-C3", "NDL-DCD", "BDL-DCD-CLD", "CROP-C4", "GRASS-C4", "BDL-DCD-DRY")

# Get results and omit duplicates that originate because runs have been restarted
# data <- read.table("daisyOutput_LAI-AVHRR")
data <- read.table("daisyOutput_GPP-GOSIF")
# data <- read.table("daisyOutput_GPP-CLASSIC")
# data <- read.table("daisyOutput_LAI-CLASSIC")
duplicated(data)
data <- data[!duplicated(data), ]
duplicated(data)

# The first 6 columns present the scores:
# "S", "S_bias", "S_rmse", "S_phase", "S_iav", "S_dist", 
# The next columns present the PFT-specific parameter values, with 9 columns per parameter
# "NDL-EVG", "BDL-EVG", "CROP-C3", "GRASS-C3", "NDL-DCD", "BDL-DCD-CLD", "CROP-C4", "GRASS-C4", "BDL-DCD-DRY"
# The last column gives the date when those values were written to the file

# Find run with highest mean score so far
id <- which.max(data[[1]])
best <- data[id,]
best.vmax <- best[7:15]
best.kn <- best[16:24]
best.list <- list(best.vmax, best.kn)
# best.list <- list(best.vmax)

colnames(best.vmax) <- PFTs
colnames(best.kn) <- PFTs

write.table(best.vmax, "best_vmax.txt")
write.table(best.kn, "best_kn.txt")

all.vmax <- data[7:15]
all.kn <- data[16:24]
all.list <- list(all.vmax, all.kn)
# all.list <- list(all.vmax)

# Create boxplots showing sampled parameter values, default parameter values, and optimized parameter value
for (i in 1:length(best.list)) {
  
  parameterName <- parameterName.list[[i]]
  best <- best.list[[i]]
  default <- default.list[[i]]
  all <- all.list[[i]]
  my.ylab <- ylab.list[[i]]
  factor <- factor.list[[i]]
  upperBound <- unlist(upperBound.list[[i]])
  lowerBound <- unlist(lowerBound.list[[i]])
  legLoc <- legLoc.list[[i]]
  
  bda <- rbind(best, default, upperBound, lowerBound, all) # bda = best, default, all
  colnames(bda) <- PFTs
  
  # Make a more intuitive order
  bda <- data.frame(
    bda["NDL-EVG"],
    bda["NDL-DCD"],
    bda["BDL-EVG"],
    bda["BDL-DCD-CLD"],
    bda["BDL-DCD-DRY"],
    bda["GRASS-C3"],
    bda["GRASS-C4"],
    bda["CROP-C3"],
    bda["CROP-C4"])

best <- as.vector(unlist(best))

# parameter values
pv <- bda[5:nrow(bda),]

q975 <- apply(X = pv, MARGIN = 2, FUN = function(x) {quantile(x,probs = c(.975),na.rm=TRUE)})
q75 <- apply(X = pv, MARGIN = 2, FUN = function(x) {quantile(x,probs = c(.75),na.rm=TRUE)})
q25 <- apply(X = pv, MARGIN = 2, FUN = function(x) {quantile(x,probs = c(.25),na.rm=TRUE)})
q025 <- apply(X = pv, MARGIN = 2, FUN = function(x) {quantile(x,probs = c(.025),na.rm=TRUE)})

quantiles <- rbind(q975, q75, q25, q025)

best <- unlist(bda[1,])
default <- unlist(bda[2,])
upperBound <- unlist(bda[3,])
lowerBound <- unlist(bda[4,])

my.ylim <- c(min(lowerBound), max(upperBound)) * factor

png(paste("parameterValues_", parameterName, ".png", sep = ""), width = 6, height = 6, units = "in", res = 300)
par(mar = c(8,5,1,1))

boxplot(quantiles * factor, lty = 1, las = 3, ylab = my.ylab, ylim = my.ylim)
points(default * factor, pch = 17, col = "blue", cex = 2)
points(best * factor, pch = 16, col = "red", cex = 2)
points(upperBound * factor, pch = 1, cex = 1)
points(lowerBound * factor, pch = 1, cex = 1)
legend(legLoc, pch = c(17, 16, 1), col = c("blue", "red", "black"), c("Default", "Optimized", "Range"), bty = "n")

dev.off()

}

# scores
scores <- data[1:6]

max <- apply(X = scores, MARGIN = 2, FUN = function(x) {max(x)})
q75 <- apply(X = scores, MARGIN = 2, FUN = function(x) {quantile(x,probs = c(.75),na.rm=TRUE)})
q25 <- apply(X = scores, MARGIN = 2, FUN = function(x) {quantile(x,probs = c(.25),na.rm=TRUE)})
min <- apply(X = scores, MARGIN = 2, FUN = function(x) {min(x)})

quantiles <- rbind(max,q75,q25,min)

colnames(quantiles) <- c("Mean Score", "Bias Score", "RMSE Score", "Seasonality Score", "Interannual Variability Score", "Spatial Correlation Score")

# add scores from best run
best <- data[id,]
best <- best[1:6]
best <- as.vector(unlist(best))

default <- unlist(read.table("DefaultParameterValuesScore"))

png("scores.png", width = 5, height = 5, units = "in", res = 300)
par(mar = c(13,5,1,1))
boxplot(quantiles, lty = 1, las = 2, ylab = "Score (-)")
# points(default, pch = 17, col = "blue", cex = 2)
points(best, pch = 16, col = "red")
legend("topleft", pch = c(17, 16), col = c("blue", "red"), c("Default", "Optimized"), bty = "n")

dev.off()

png("scores_difference.png", width = 5, height = 5, units = "in", res = 300)
par(mar = c(13,5,1,1))
delta <- best - default

my.ylab <- latex2exp::TeX("$\\Delta$ Score (opt. minus default) (-)")

my.xlab <- c("Mean Score", "Bias Score", "RMSE Score", "Seasonality Score", "Interannual Variability Score", "Spatial Correlation Score")
barplot(delta, names.arg = my.xlab, las = 2, ylab = my.ylab, ylim = c(0,0.05))
box()
dev.off()


library(GA)
object <- readRDS("objectFinal.rds")
df <- data.frame(object@summary)

png("score_generation.png", width = 6, height = 5, units = "in", res = 300)
par(mar = c(5,5,1,1))
plot(df$median, type = "l", xlim = c(1,50), ylim = c(0.65, 0.672),
xlab = "Generation", ylab = "Score (-)")
#points(df$median, pch = 16)
lines(df$q1, type = "l", col = "grey")
lines(df$q3, type = "l", col = "grey")

#lines(df$min, type = "l", lty = 2, col = "grey")
#lines(df$max, type = "l", lty = 2, col = "grey")

# abline(h = 0.7191567, col = "blue")
abline(h = 0.664981976460085, col = "blue") # Default 200 GC
# legend("bottomright", col = c("blue", "grey", NA), lty = 1, lwd = 2, 
# c("Score when using default vmax and kn values", "Interquartile range", "Population Size = 20"), bty = "n")

dev.off()

# Obtain solutions

parameterNames <- list("vmax", "kn", "thrprcnt", "alpha_phtsyn", "lfespany", "grescoef", 
                       "avertmas", "mxrtdpth", "sn", "XLEAF", "maxage", "ZOLNG", 
                       "coldthrs", "alpha", "gamma_w", "beta2", "kappa", "minslai", 
                       "coldlmt", "roothrsh", "abar", "vpd0", "albnir", "TCSAND", 
                       "TCCLAY", "ZOLNS", "minlvfr", "omega")


solution <- object@solution

parameterValues <- intFun.unnormalize(normParameterValuesList[[i]], upperBoundList[[i]], lowerBoundList[[i]])


