# EnsembleForcing
Generates spatially and temporally correlated random patterns to create ensemble forcing. Tested with the UEB (Utah Energy Balance) snowmelt model and the Noah-MP land surface models.
The code creates multiplicative forcing error for precipitation, longwave radition, humidity, and wind speed; and addive forcing error for temperature and solar radiation.
Accounts for spatial and temporal correlation, and correlation among precipitation, solar and long wave raditions accounted for. 

# Dependencides
Eigen https://eigen.tuxfamily.org/index.php?title=Main_Page 

netCDF https://www.unidata.ucar.edu/software/netcdf/

any MPI library, e.g., MPICH https://www.mpich.org/
