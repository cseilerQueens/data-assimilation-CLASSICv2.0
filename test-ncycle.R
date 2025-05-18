
.libPaths(c("/home/cseiler/daisy/renv", .libPaths()))
library(daisy)
library(terra)
rm(list = ls())

setwd('/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0')

soilpH <- terra::rast("NETCDF:/home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/initfiles/global/T63/rsFilev04.nc:soilpH")
tbar <- 273-50

soipH <- -1.175494e+38
KH <- 4.59 * tbar*exp(4092.* ( (1./tbar) - (1./298.15) ) )
KNH4 <- 5.67*1.E-10*exp(-6286.*((1. / tbar) - (1. / 298.15)))
D = 1.0 + KH + 10.**(- soilpH) * (KH / KNH4)
plot(D)



tbar.start <- -20 + 273.15
tbar.end <- 20 + 273.15
tbar <- seq(tbar.start, tbar.end, 0.1)

soilpH <- -1.175494e+38
KH <- 4.59 * tbar*exp(4092.* ( (1./tbar) - (1./298.15) ) )
KNH4 <- 5.67*1.E-10*exp(-6286.*((1. / tbar) - (1. / 298.15)))
D = 1.0 + KH + 10.**(- soilpH) * (KH / KNH4)
plot(D)



KH <- 4.59 * tbar*exp(4092.* ( (1./tbar) - (1./298.15) ) )

plot(tbar, KH, type = "l")
abline(v=273.15)

# tbar <- 271.738666793401
# KH <- 4734.68883219132
# KNH4 <- 7.305314159898535E-011

tbar <- 300

KH <- 4.59 * tbar*exp(4092.* ( (1./tbar) - (1./298.15) ) )
KNH4 <- 5.67*1.E-10*exp(-6286.*((1. / tbar) - (1. / 298.15)))
# KH <- 0
soilpH <- seq(0.1,8.8, 0.1)
D = 1.0 + KH + 10.**(- soilpH) * (KH / KNH4)

plot(soilpH, D, type = "l")
