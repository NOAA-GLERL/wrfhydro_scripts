#!/bin/bash
# Download NWM forcing from the google cloud 
# contact: james.kessler@noaa.gov

base_url='gs://national-water-model'

# specify the dimensions and variables to subset
dim_arg='-d x,2250,4250 -d y,1600,3200'  # from geogrid file history for the "GL Bigger" Domain


# set date to grab
#YYYYMMDD=20190319

# or uncomment below and pass date as an argument (a la "./get_nwm_forcing.sh 20190319")
YYYYMMDD=$1



mkdir $YYYYMMDD analysis 2> /dev/null

# get analysis forcing
echo downloading analysis forcing...
gsutil -qm cp ${base_url}/nwm.${YYYYMMDD}/forcing_analysis_assim/nwm.t??z.analysis_assim.forcing.tm00.conus.nc $YYYYMMDD # just analysis

echo clipping...
for f in $(find $YYYYMMDD -type f -iname '*.nc'); do
	ncks -O $dim_arg $f $f
	#ncpdq -O -U $f $f # (dont think this is necessary)
done
echo concatenating and cleaning up
ncrcat $YYYYMMDD/*.nc analysis/$YYYYMMDD.nc
rm -rf $YYYYMMDD



# not fully implemented; get short/medium range forcing
#HH=01
#gsutil -qm cp ${base_url}/nwm.${YYYYMMDD}/forcing_short_range/nwm.t${HH}z.short_range.forcing.f{001..018}.conus.nc $YYYYMMDD

# get medium range term forecast
#echo downloading medium range forcings...
#gsutil -qm cp ${base_url}/nwm.${YYYYMMDD}/forcing_medium_range/nwm.t${HH}z.medium_range.forcing.f{001..240}.conus.nc $YYYYMMDD

## advance date/time (for grabbing multiple FCs in a loop):
#HH=$(( HH+1 ))  
#YYYYMMDD=$(date -d "$YYYYMMDD + 1 day" +%Y%m%d)
