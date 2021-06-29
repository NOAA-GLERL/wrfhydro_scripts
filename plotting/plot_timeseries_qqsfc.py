# This script reads EOF metadata from the "master list". It then reads in
# the RTOUT files and pulls surface runoff from each timestep to create a 
# continuous record. The script also pulls in the observation text files
# and then plots a time series of the observations and simulated runoff.

# Developer: Lindsay Fitzpatrick 
# Contact: ljob@umich.edu
# Version: 1
# Updated: 06-29-21

from pylab import *
import numpy as np
from netCDF4 import Dataset, num2date
from datetime import datetime
import numpy.ma as ma
from matplotlib import dates
import xarray as xr
from netCDF4 import Dataset as NetCDFFile
import pandas as pd
import os

dir = '/glade/scratch/lfitzpatrick/Retrospec/'
gname='/glade/u/home/lfitzpatrick/Data/eof_i_j.txt'

#read in EOF master list 
for t in range(0,72):
	data=np.loadtxt(gname,skiprows=1,dtype={'names':('line','site','lat','lon','date_start',
					'date_end','i_land','j_land','i_routing','j_routing'),
                                        'formats':('S2','S8',np.float,np.float,'S10','S10',np.int,np.int,np.int,np.int)})
	
	line = data['line'][t]
	lat = data['lat'][t]
	lon = data['lon'][t]
	i=data['i_routing'][t]
	j=data['j_routing'][t]
	i2=data['i_land'][t]
        j2=data['j_land'][t]
	site=data['site'][t]
	date_start=data['date_start'][t]
	date_end=data['date_end'][t]	

	print(line, site, lat, lon, date_start, date_end, i2, j2, i, j)

	#opens observation file
	fname='/glade/u/home/lfitzpatrick/Data/'+site+'.txt'
	data=np.loadtxt(fname,skiprows=3,dtype={'names':('date','runoff'),'formats':('S8',np.float)})
	rn_obs=data['runoff'] *25.4 #convert in to mm

	#create an empty array
	qqsfc_tdon = np.empty(5000)

	#open netcdfs and pull runoff from each file
	h = 0
	for f in sorted(os.listdir(dir)):
		if f.startswith('201'):
			if f.endswith('.RTOUT_DOMAIN1_TDON.nc'):
        			infile = xr.open_dataset(os.path.join(dir, f))
				qqsfc_tdon[h] = infile['qqsfc_acc'][j,i]
				h = h + 1
				infile.close()

	max_sfc = max(qqsfc_tdon)
	max_obs = max(rn_obs)

	if max_obs > max_sfc: max_val = max_obs
	else: max_val = max_sfc
		
	m=np.arange(0,len(rn_obs),1)
	x=num2date(m,units='days since ' +date_start+ ' 00:00',calendar='standard')
	n=np.arange(0,len(qqsfc_tdon),1)
	x2=num2date(n,units='days since 2010-01-01 00:00',calendar='standard')
	
	#create timeseries plot
	fig, ax = plt.subplots(figsize=(15,3))
	
	title(site+" Surface Runoff")
	plot_date(x, rn_obs,color='black',ls='-',linewidth=2,marker=None,label="EOF Obs")
	plot_date(x2, qqsfc_tdon,color='red',ls='-',marker=None,label="TD ON")
	ax.set_ylim([0,60])
	ax.set_ylabel('[mm]')
	ax.set_xlabel('Month/Year')
	ax.set_xlim([datetime(2010,1,1),datetime(2015,1,1)])
	hfmt = dates.DateFormatter('%m/%Y')
	fig.autofmt_xdate()
	ax.xaxis.set_major_formatter(hfmt)
	plt.xticks(rotation='45')
	ax.legend(loc=1,fontsize=14)
	
	#save figure
	savefig(site+'_runoff.png', bbox_inches='tight')
	print('done with ', site)
	plt.close()


