
.libPaths("/home/cseiler/daisy/renv")
library(daisy)
setwd("/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0")

rm(list = ls())
# Obtain model and reference data
dir.mod <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/DefaultDaisyRun200GC/netcdf_files"
dir.ref <- "/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63"

nc.mod01 <- file.path(dir.mod, "gpp_monthly.nc")
nc.mod02 <- file.path(dir.mod, "lai_monthly.nc")
nc.mod03 <- file.path(dir.mod, "hfls_monthly.nc")
nc.mod04 <- file.path(dir.mod, "hfss_monthly.nc")
nc.mod05 <- file.path(dir.mod, "ts_monthly.nc")
nc.mod06 <- file.path(dir.mod, "albs_monthly.nc")

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

mod.list <- list(
  nc.mod01, nc.mod01,
  nc.mod02, nc.mod02, nc.mod02,
  nc.mod03, nc.mod03,
  nc.mod04, nc.mod04,
  nc.mod05,
  nc.mod06, nc.mod06, nc.mod06
)

ref.list <- list(
  nc.ref01a, nc.ref01b,
  nc.ref02a, nc.ref02b, nc.ref02c,
  nc.ref03a, nc.ref03b,
  nc.ref04a, nc.ref04b,
  nc.ref05a,
  nc.ref06a, nc.ref06b, nc.ref06c
)

# Unit conversion factors for reference data

ref.unit.conv.list <- list(
  1/86400000, 1/86400000,
  1,1,1,
  1, (10^6/86400),
  1, (10^6/86400),
  1,
  1,1,1
)


ref.id.list <- list("GPP-GOSIF", "GPP-FLUXCOM",
                    "LAI-AVHRR", "LAI-Copernicus", "LAI-MODIS",
                    "HFLS-CLASSr", "HFLS-FLUXCOM",
                    "HFSS-CLASSr", "HFSS-FLUXCOM",
                    "TS-MODIS",
                    "ALBS-CERES", "ALBS-GEWEXSRB", "ALBS-MODIS")

    # Check whether the length of the model and reference lists are equal
    condition <- length(mod.list) != length(ref.list)
    if (condition) {
      stop("Error: The length of the model and reference lists are not equal")
    }
    
    n <- length(ref.list)
    scores <- numeric(n)
    
    daisyOutput <- numeric()
    
    # Loop through all reference datasets
    for (i in 1:n) {
      
      nc.mod <- mod.list[[i]]
      nc.ref <- ref.list[[i]]
      ref.id <- ref.id.list[[i]]
      ref.unit.conv <- ref.unit.conv.list[[i]]
      
      myLevel <- 1
      
      mod <- raster::brick(nc.mod, level = myLevel)
      suppressWarnings(mod <- raster::rotate(mod))
      
      ref <- raster::brick(nc.ref)
      suppressWarnings(ref <- raster::rotate(ref))
      
      # model data dates
      dates.mod <- intFun.getZ(nc.mod)
      mod <- raster::setZ(mod, dates.mod)
      names(mod) <- dates.mod
      dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
      start.date.mod <- min(dates.mod)
      end.date.mod <- max(dates.mod)
      
      # reference data dates
      dates.ref <- intFun.getZ(nc.ref)
      ref <- raster::setZ(ref, dates.ref)
      names(ref) <- dates.ref
      dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
      start.date.ref <- min(dates.ref)
      end.date.ref <- max(dates.ref)
      
      # find common time period
      start.date <- max(start.date.mod, start.date.ref)
      end.date <- min(end.date.mod, end.date.ref)
      
      # subset common time period
      mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date &
                          format(as.Date(raster::getZ(mod)), "%Y-%m") <= end.date)]]
      ref <- ref[[which(format(as.Date(raster::getZ(ref)), "%Y-%m") >= start.date &
                          format(as.Date(raster::getZ(ref)), "%Y-%m") <= end.date)]]
      
      # get layer names
      mod.names <- base::names(mod)
      ref.names <- base::names(ref)
      
      # unit conversion if appropriate
      ref <- ref * ref.unit.conv
      
      #---------------------------------------------------------------------------
      
      # II Statistical analysis
      
      #---------------------------------------------------------------------------
      
      # (1) Bias
      
      #---------------------------------------------------------------------------
      # create a mask to excludes all grid cells that the model and reference
      # data do not have in common.  This mask varies in time.
      mask <- (mod * ref)
      mask <- mask - mask + 1
      mod <- mod * mask
      names(mod) <- mod.names  # this adds the corresponding dates
      ref <- ref * mask
      names(ref) <- ref.names  # this adds the corresponding dates
      # now mod and ref are based on the same grid cells
      
      mod.mean <- raster::mean(mod, na.rm = TRUE)  # time mean
      ref.mean <- raster::mean(ref, na.rm = TRUE)  # time mean
      bias <- mod.mean - ref.mean  # time mean
      
      mod.sd <- raster::calc(mod, fun = sd, na.rm = TRUE)  # standard deviation of model data
      ref.sd <- raster::calc(ref, fun = sd, na.rm = TRUE)  # standard deviation of reference data
      epsilon_bias <- abs(bias)/ref.sd
      epsilon_bias[epsilon_bias == Inf] <- NA  # relative error
      bias.score <- exp(-1 * epsilon_bias)  # bias score as a function of space
      S_bias <- mean(raster::getValues(bias.score), na.rm = TRUE)  # scalar score
      
      #---------------------------------------------------------------------------
      
      # (2) root mean square error (rmse)
      
      #---------------------------------------------------------------------------
      ESM.mode <- FALSE
      # Set ESM.mode to FALSE when forcing model with quasi-observed data
      # (e.g.  reanalysis)
      if (ESM.mode == FALSE) {
        rmse <- intFun.rmse(mod, ref)  # rmse
        mod.anom <- mod - mod.mean  # anomaly
        ref.anom <- ref - ref.mean  # anomaly
      }
      
      # Set ESM.mode to TRUE for online experiments or when forcing model
      # with earth system model data
      if (ESM.mode == TRUE) {
        index <- format(as.Date(names(ref), format = "X%Y.%m.%d"), format = "%m")
        index <- as.numeric(index)
        mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
        ref.clim.mly <- raster::stackApply(ref, index, fun = raster::mean)
        
        mod.cycle <- mod.clim.mly  # does not necessarily start in Jan or end in Dec
        ref.cycle <- ref.clim.mly  # does not necessarily start in Jan or end in Dec
        
        # Ensure that the order is correct
        JanToDec <- c("index_1", "index_2", "index_3", "index_4", "index_5",
                      "index_6", "index_7", "index_8", "index_9", "index_10", "index_11",
                      "index_12")
        mod.clim.mly <- mod.clim.mly[[JanToDec]]
        ref.clim.mly <- ref.clim.mly[[JanToDec]]
        
        rmse <- intFun.rmse(mod.clim.mly, ref.clim.mly)  # rmse
        
        mod.anom <- mod.clim.mly - mod.mean  # anomaly
        ref.anom <- ref.clim.mly - ref.mean  # anomaly
        
      }
      
      crmse <- intFun.crmse(mod.anom, ref.anom)  # centralized rmse
      
      #-------------------------------------------------------------
      epsilon_rmse <- crmse/ref.sd
      epsilon_rmse[epsilon_rmse == Inf] <- NA  # relative error
      rmse.score <- exp(-1 * epsilon_rmse)  # rmse score as a function of space
      S_rmse <- mean(raster::getValues(rmse.score), na.rm = TRUE)
      
      #---------------------------------------------------------------------------
      
      # (3) phase shift
      
      #---------------------------------------------------------------------------
      
      index <- format(as.Date(names(ref), format = "X%Y.%m.%d"), format = "%m")
      index <- as.numeric(index)
      mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
      ref.clim.mly <- raster::stackApply(ref, index, fun = raster::mean)
      
      mod.cycle <- mod.clim.mly  # does not necessarily start in Jan or end in Dec
      ref.cycle <- ref.clim.mly  # does not necessarily start in Jan or end in Dec
      
      # Ensure that the order is correct
      JanToDec <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6",
                    "index_7", "index_8", "index_9", "index_10", "index_11", "index_12")
      mod.clim.mly <- mod.clim.mly[[JanToDec]]
      ref.clim.mly <- ref.clim.mly[[JanToDec]]
      
      bias.clim.mly <- mod.clim.mly - ref.clim.mly
      
      # find month of seasonal peak
      
      # In most cases, we are interested in the timing of the seasonal
      # maximum value
      
      # In some cases, however, the seasonal peak is a minimum, e.g. NEE =
      # RECO - GPP
      
      phaseMinMax <- "phaseMax"
      #
      if (phaseMinMax == "phaseMax") {
        mod.max.month <- raster::which.max(mod.clim.mly)
        ref.max.month <- raster::which.max(ref.clim.mly)
      }
      
      if (phaseMinMax == "phaseMin") {
        mod.max.month <- raster::which.min(mod.clim.mly)
        ref.max.month <- raster::which.min(ref.clim.mly)
      }
      
      # get shortest time distance between these months
      abs.diff <- abs(mod.max.month - ref.max.month)  # absolute difference from 0 to 12 months
      phase <- raster::calc(abs.diff, intFun.theta)  # shortest distance from 0 to 6 months (theta)
      phase.score <- 0.5 * (1 + cos(2 * pi * phase/12))  # score from 0 (6 months) to 1 (0 months)
      S_phase <- mean(raster::getValues(phase.score), na.rm = TRUE)  # scalar score
      
      #---------------------------------------------------------------------------
      
      # (4) interannual variability
      
      #---------------------------------------------------------------------------
      
      years <- floor(raster::nlayers(mod)/12)  # total number of years
      months <- years * 12  # number of months considering complete years only
      mod.fullyear <- raster::subset(mod, 1:months)
      ref.fullyear <- raster::subset(ref, 1:months)
      c.mod <- raster::calc(mod.cycle, fun = function(x) {
        rep(x, years)
      })  # climatological cycle for all months (mod)
      c.ref <- raster::calc(ref.cycle, fun = function(x) {
        rep(x, years)
      })  # climatological cycle for all months (ref)
      mod.iav <- sqrt(raster::mean((mod.fullyear - c.mod)^2, na.rm = TRUE))  # interannual variability  (mod)
      ref.iav <- sqrt(raster::mean((ref.fullyear - c.ref)^2, na.rm = TRUE))  # interannual variability  (ref)
      # set values close to zero to NA
      # ref.iav.na <- ref.iav
      # ref.iav.na[ref.iav.na < 10^(-5)] <- NA
      # ref.iav.na[ref.iav.na < 10^(-10)] <- NA
      
      # epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav.na
      epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav
      
      epsilon_iav[epsilon_iav == Inf] <- NA  # I changed Eq. 26 so that epsilon_iav >= 0
      iav.score <- exp(-1 * epsilon_iav)  # iav score as a function of space
      S_iav <- mean(raster::getValues(iav.score), na.rm = TRUE)  # scalar score (not weighted)
      
      #---------------------------------------------------------------------------
      
      # (5) dist
      
      #---------------------------------------------------------------------------
      
      mod.sigma.scalar <- stats::sd(raster::getValues(mod.mean), na.rm = TRUE)  # standard deviation of period mean data
      ref.sigma.scalar <- stats::sd(raster::getValues(ref.mean), na.rm = TRUE)  # standard deviation of period mean data
      sigma <- mod.sigma.scalar/ref.sigma.scalar
      y <- raster::getValues(mod.mean)
      x <- raster::getValues(ref.mean)
      reg <- stats::lm(y ~ x)
      R <- sqrt(summary(reg)$r.squared)
      S_dist <- 2 * (1 + R)/(sigma + 1/sigma)^2
      
      #---------------------------------------------------------------------------
      
      # overall scores
      
      #---------------------------------------------------------------------------
      
      S <- mean(c(S_bias, S_rmse, S_phase, S_iav, S_dist), na.rm = TRUE)
      
      timeStamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      
      # Write current scores and parameter values to text file for analysis
      # purpose
      fileName <- paste("daisyOutput", ref.id, sep = "_")
      write.table(x = matrix(c(S, S_bias, S_rmse, S_phase, S_iav, S_dist, daisyOutput,
                               timeStamp), nrow = 1), file = fileName, append = TRUE, row.names = FALSE,
                  col.names = FALSE)
      
      scores[i] <- S
    }
    
    # Calculate the mean score across all evaluations
    S <- mean(scores, na.rm = TRUE)

    # Write current scores and parameter values to text file for analysis purpose
    write.table(x = matrix(c(S, S_bias, S_rmse, S_phase, S_iav, S_dist), nrow = 1), file = "DefaultParameterValuesScore200GC", row.names = FALSE, col.names = FALSE)
