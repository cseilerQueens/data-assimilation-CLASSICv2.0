
module load StdEnv/2023 gcc/12.3 udunits/2.2.28 hdf/4.2.16 gdal/3.7.2 r/4.3.1 nco/5.1.7

rm GC_multipleGridCells.nc
rm masked_data.nc
rm *ncks.tmp
rm tmp*.nc

# cp /home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/initfiles/global/T63/rsFile_CLASSICv2.0_Ncycle_T63.nc rsFile.nc
cp /home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/spinup_CRUJRAv2.4.5-part04-NCycle-coupled/netcdf_files/rsFile_modified.nc rsFile.nc
Rscript selectGridcells.R

# module load openmpi/4.0.3
# module load cdo/1.9.8
# module load nco/5.0.6
# module load netcdf-fortran-mpi/4.5.2
# module load hdf5-mpi/1.10.6


# module load nco/5.0.6 StdEnv/2020 intel/2020.1.217
# module load cdo/2.0.5 StdEnv/2020 intel/2020.1.217 openmpi/4.0.3
ncks -4 tmp_GC.nc tmp01.nc

# cp /home/cseiler/projects/def-cseiler-ab/cseiler/classic_inputs/initfiles/global/T63/rsFile.nc .
# cp /home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/transient_CRUJRAv2.4.5_2000/netcdf_files/rsFile_modified.nc rsFile.nc
ncks -A -v FLND tmp01.nc rsFile.nc
ncap2 -s 'where(FLND < 1) FLND=0' rsFile.nc tmp02.nc
# convert zero to NA, which shows up as "_"
ncatted -a _FillValue,FLND,o,f,0.0 tmp02.nc tmp03.nc
# ncdump -v FLND tmp02.nc

mv tmp03.nc rsFile_daisy.nc
rm tmp*.nc
