# This script can be used to create parameter uncertainty ranges that are +/- a given percentage of the default value

.libPaths("/home/cseiler/daisy/renv")
library(daisy)
library(GA)

rm(list = ls())

setwd('/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0')

# (1) Get the parameterValueLength

# Copy default parameter run_parameters.nml file to current work directory
# This file will be overwritten with new values before CLASSIC is run

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

relUncRan <- 0.2 # relative uncertainty range
  
for (i in parameterNames) {
  parameters <- getParameterValues(parameterFile = parameterFile, parameterName = i)
  parval <- c(parameters)
  names(parval) <- paste0("X", 1:length(parval))
  parval <- as.data.frame(t(parval))
  
  absUncRan <- abs(parval) * relUncRan  # absolute uncertainty range
  
  parval.min <- parval - absUncRan
  parval.max <- parval + absUncRan

  fname <- paste(i, "min.csv", sep = ".")
  write.csv(x = parval.min, file = fname, row.names = FALSE)
  
  fname <- paste(i, "max.csv", sep = ".")
  write.csv(x = parval.max, file = fname, row.names = FALSE)
  }

