.libPaths(c("/home/cseiler/daisy/renv", .libPaths()))
library(daisy)
rm(list = ls())

setwd('/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0')

library(raster)

nc.mod <- "CLASSIC_meta/RUN1/daisyRun/netcdf_files/gpp_monthly.nc"
nc.ref <- "/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63/gpp_GOSIF_128x64.nc"

period <- c("1980-01", "2020-12")

# Ensure that the model output has the correct spatial extent
# This is vital for CLASSIC-v2.0

nc <- ncdf4::nc_open(nc.mod)

# Extract longitude and latitude values
lon <- ncdf4::ncvar_get(nc, "longitude")
lat <- ncdf4::ncvar_get(nc, "latitude")

ncdf4::nc_close(nc)

xres <- mean(diff(lon))
yres <- mean(diff(lat))

xmin <- min(lon) - xres / 2
xmax <- max(lon) + xres / 2
ymin <- min(lat) - yres / 2
ymax <- max(lat) + yres / 2

mod.extent <- extent(xmin, xmax, ymin, ymax)




  
    #---------------------------------------------------------------------------  
    
    # Part (I) Compute scores for all variables with reference data that vary in time
    
    #---------------------------------------------------------------------------
    
mod <- raster::brick(nc.mod)
ref <- raster::brick(nc.ref)

raster::extent(mod) <- raster::extent(mod.extent)

raster::extent(mod) <- raster::extent(mod.extent)

mod.mean <- mean(mod)
ref.mean <- mean(ref)

plot(ref.mean)
plot(mod.mean, add = TRUE, col = "black")

z <- mod.mean + ref.mean

      
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
      
      # Part (II) Compute scores for all variables with reference data that do not vary in time
      
      #---------------------------------------------------------------------------
      
    } else {
      
      # model data are averaged over a period since the reference data set
      # has no time dimension
      start.date <- period[1]
      end.date <- period[2]
      
      # model data dates
      dates.mod <- intFun.getZ(nc.mod)
      mod <- raster::setZ(mod, dates.mod, name = "date")
      names(mod) <- dates.mod
      
      mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date &
                          format(as.Date(raster::getZ(mod)), "%Y-%m") <= end.date)]]
      
      mod <- raster::mean(mod, na.rm = TRUE)
      
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
      ref <- ref * mask
      
      #---------------------------------------------------------------------------
      
      mod.mean <- mod  # data is already time mean
      ref.mean <- ref  # data is already time mean
      bias <- mod.mean - ref.mean  # time mean
      epsilon_bias <- abs(bias)/abs(ref)
      epsilon_bias[epsilon_bias == Inf] <- NA  # relative error
      bias.score <- exp(-1 * epsilon_bias)  # bias score as a function of space
      S_bias <- mean(raster::getValues(bias.score), na.rm = TRUE)  # scalar score
      
      #---------------------------------------------------------------------------
      
      # (5) dist
      
      #---------------------------------------------------------------------------
      mod.sigma.scalar <- sd(raster::getValues(mod.mean), na.rm = TRUE)  # standard deviation of period mean data
      ref.sigma.scalar <- sd(raster::getValues(ref.mean), na.rm = TRUE)  # standard deviation of period mean data
      sigma <- mod.sigma.scalar/ref.sigma.scalar
      y <- raster::getValues(mod.mean)
      x <- raster::getValues(ref.mean)
      reg <- stats::lm(y ~ x)
      R <- sqrt(summary(reg)$r.squared)
      S_dist <- 2 * (1 + R)/(sigma + 1/sigma)^2  # weighting does not apply
      
      #---------------------------------------------------------------------------
      
      # Scores
      
      #---------------------------------------------------------------------------
      
      S_rmse <- NA
      S_phase <- NA
      S_iav <- NA
      
    }
    
    
    
    #---------------------------------------------------------------------------
    
    # overall scores
    
    #---------------------------------------------------------------------------
    
    S <- mean(c(S_bias, S_rmse, S_phase, S_iav, S_dist), na.rm = TRUE)
    
    timeStamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    
    # Since the cost.fun file gets executed from a different directory in parallel
    # than sequentially, new score files get created every run. To avoid this,
    # the directory is changed when creating these score files and then reverted
    # to what it was before.
    origDir <- getwd()
    setwd(dataAssimPath)
    Sys.sleep(1)
    # Write current scores and parameter values to text file for analysis
    # purpose
    fileName <- paste("daisyOutput", ref.id, sep = "_")
    write.table(x = matrix(c(S, S_bias, S_rmse, S_phase, S_iav, S_dist, daisyOutput,
                             timeStamp), nrow = 1), file = fileName, append = TRUE, row.names = FALSE,
                col.names = FALSE)
    # Revert to original directory.
    setwd(origDir)
    Sys.sleep(1)
    
    scores[i] <- S
  }
  
  # Move outputs to a new folder so that they are not overwritten.
  if (keepOutput) {
    NewModelOutputFolder <- paste(modelOutputFolder, paste0("RUN", identifier), sep = "_")
    NewModelOutputFolder <- paste(NewModelOutputFolder, timeStamp, sep = "_")
    moveModelOutputFolder <- paste("mv", modelOutputFolder, NewModelOutputFolder,
                                   sep = " ")
    system(moveModelOutputFolder)
    Sys.sleep(1)
    
    # Move the finished run into the final output folder.
    
    moveModelOutputFolder <- paste("mv", NewModelOutputFolder, finalOutputFolder)
    system(moveModelOutputFolder)
    Sys.sleep(1)
    # Make a file containing the name of the output folder for reference.
    runFolder <- unlist(strsplit(NewModelOutputFolder, "/"))
    runFolder <- runFolder[length(runFolder)]
    runFolder <- paste(finalOutputFolder, runFolder, sep = "/")
    write(runFolder, "simInfo.txt")
    Sys.sleep(1)
    
    # Delete output folder.
  } else {
    system(paste("rm -rf", modelOutputFolder))
    Sys.sleep(1)
  }
  
  # Calculate the mean score across all evaluations.
  S <- mean(scores, na.rm = TRUE)
  write(S, "simInfo.txt", append = TRUE)
  Sys.sleep(2)
  
  return(S)
}