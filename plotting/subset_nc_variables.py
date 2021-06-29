# This script reads in the NetCDF files and reduces the number of variables
# in it before resaving the file. This reduces the size of the files.

# Developer: Lindsay Fitzpatrick 
# Contact: ljob@umich.edu
# Version: 1
# Updated: 06-29-21

import os, time
import netCDF4
from netCDF4 import Dataset
import xarray as xr

# This script reads in WRF-Hydro output files as netcdfs and reduces the number of
# variables before saving as new netcdfs. This reduces the size of the files.

dir = '/glade/scratch/lread/GreatLakes_Retros/td/'

#counter
t = 0
#read in netcdf files
for f in sorted(os.listdir(dir)):
	if f.startswith('2012'):
		if f.endswith('.LDASOUT_DOMAIN1'):
                	start = time.time()
			infile = xr.open_dataset(os.path.join(dir, f))
			ds = infile.copy(deep=True)
			
			#drop the variables that you do NOT want to save
			ds_out = ds.drop(['IVGTYP', 'ISLTYP', 'FVEG', 'LAI', 'SAI', 'COSZ',
                           	        'GRDFLX', 'HFX', 'LH', 'ECAN', 'EDIR', 'ALBEDO', 'ETRAN', 'UGDRNOFF',
                                	'SFCRNOFF', 'CANLIQ', 'CANICE', 'ZWT', 'WA', 'WT', 'ACCECAN', 'ACCEDIR',
                                  	'ACCETRAN', 'SAV', 'TR', 'EVC', 'IRC', 'SHC', 'IRG', 'SHG', 'EVG', 'GHV', 'SAG',
                                  	'SHB', 'EVB', 'GHB', 'TRAD', 'TG', 'TV', 'TAH', 'TGV', 'T2MB',
                                 	'Q2MB', 'EAH', 'FWET', 'ZSNSO_SN', 'SNICE', 'SNLIQ',
                                  	'SNOW_T', 'SOIL_M', 'SNOWH', 'QSNOW', 'ISNOW', 'FSNO', 'ACSNOW',
                                  	'CM', 'CH', 'CHV', 'CHB', 'CHLEAF', 'CHUC', 'CHV2', 'CHB2', 'LFMASS', 'RTMASS',
                                  	'STMASS', 'WOOD', 'STBLCP', 'FASTCP', 'NEE', 'GPP', 'NPP', 'PSN', 'APAR', 'ACCET',
                                  	'CANWAT', 'SNOWT_AVG', 'ALBSND', 'ALBSNI', 'QRAIN'])
                 
			#name the output file
			outfile = f + '_ACCPRCP.nc'
			#save file
			ds_out.to_netcdf(outfile)
			print ("done with " +str(outfile))
                	t = t + 1
                	# close the netcdf file
                	infile.close()
			ds.close()
		
			end = time.time()
			print("Time: {0:.4f}s".format(end - start))




