#!/bin/bash

base_url='gs://national-water-model'
#cd /path/to/where/you/want/to/save/data

# define desired feature_ids and variables to process
#                              add -d fid to newline;rm trail ws; add decimal | remove newlines
id_file='/mnt/projects/hpc/kessler/lakes/flow/nwm/more_ids/hecwfs_ids_all.txt'
dim_arg=$(sed -r 's/^/-d feature_id,/g; s/\s+$//g; s/$/.0/g' $id_file | tr '\n' ' ')
var_arg='-v streamflow,velocity,feature_id,time'

t0=20210423        # first record to grab
tf=$(date +%Y%m%d) # last record to grab (default: today)
tf=20210423        # or specified date

anal=1  # get analysis?
shrt=1  # get short fc?
med=0   # get med fc?



YYYYMMDD=$t0
# =====================================================================================================
# ============================   get analysis =========================================================
# =====================================================================================================

if [ $anal == 1 ]; then
echo downloading/processing analysis files... 

# had to implement a new file structure to match Dans so i could use his preprocessing....

rm -rf temp 2> /dev/null
while [[ $YYYYMMDD -le $tf ]]; do
	echo $YYYYMMDD
	YYYY=${YYYYMMDD:0:4}
    [[ -d analysis/$YYYY ]] || mkdir -p  analysis/$YYYY
	MMDD=${YYYYMMDD:4:4}

	#1. download & process streamflows (channel_rt files)
	mkdir temp
	gsutil -qm cp ${base_url}/nwm.${YYYYMMDD}/analysis_assim/nwm.t{00..23}z.analysis_assim.channel_rt.tm00.conus.nc temp
	/home/kessler/.local/bin/cdo --no_warnings copy temp/*channel_rt*.nc temp/${YYYYMMDD}.nc 2> /dev/null # concatenate
	/usr/bin/ncks -O $dim_arg $var_arg temp/${YYYYMMDD}.nc  temp/${YYYYMMDD}.nc         # clip
	mv temp/${YYYYMMDD}.nc analysis/${YYYY}/${MMDD}.nc  
	rm -rf temp
	chmod +r analysis/${YYYY}/${MMDD}.nc		

	YYYYMMDD=$(date -d "$YYYYMMDD + 1 day" +%Y%m%d)
done
fi


# =====================================================================================================
# ============================  get 18 hr forecast (init hours 0, 12) ==================================
# =====================================================================================================

if [ $shrt == 1 ]; then
echo downloading 18 hour FC\'s...
YYYYMMDD=$t0
rm -rf temp 2> /dev/null
while [[ $YYYYMMDD -le $tf ]]; do
	echo -n $YYYYMMDD:
	YYYY=${YYYYMMDD:0:4}
    [[ -d short/$YYYY ]] || mkdir -p short/$YYYY
	MMDD=${YYYYMMDD:4:4}

	for HH in {00,12}; do 
		echo -n ${HH}Z...

		#download
		mkdir temp
		gsutil -qm cp ${base_url}/nwm.${YYYYMMDD}/short_range/nwm.t${HH}z.short_range.channel_rt.f{001..018}.conus.nc temp

		# process (concatenate; clip; unpack)
		/home/kessler/.local/bin/cdo --no_warnings copy temp/*.nc temp/${YYYYMMDD}${HH}.nc 2> /dev/null # concatenate
		/usr/bin/ncks -O $dim_arg $var_arg temp/${YYYYMMDD}${HH}.nc  temp/${YYYYMMDD}${HH}.nc  
		/usr/bin/ncpdq -O -U temp/${YYYYMMDD}${HH}.nc  temp/${YYYYMMDD}${HH}.nc  

		# save final and cleanup
		mv temp/${YYYYMMDD}${HH}.nc short/${YYYY}/${MMDD}${HH}.nc  
		chmod +r short/${YYYY}/${MMDD}${HH}.nc		
		rm -rf temp
	done
	YYYYMMDD=$(date -d "$YYYYMMDD + 1 day" +%Y%m%d)
	echo ' '
done
fi

## =====================================================================================================
## ====================  get medium range forecast (init hours 0, 12) ==================================
## =====================================================================================================

if [ $med == 1 ]; then
echo downloading medium range FC\'s...
# remove any lingering downloaded Forecast files
rm -rf download_m?

YYYYMMDD=$t0
while [[ $YYYYMMDD -le $tf ]]; do # DAY LOOP
	echo $YYYYMMDD...
    YYYY=${YYYYMMDD:0:4}
    MMDD=${YYYYMMDD:4:4}

	for HH in {00,12}; do  # HOUR LOOP
		echo -n ${HH}Z:

		for mem in {1..7}; do # ENSEMBLE MEMBER LOOP
			echo -n m${mem}...
            [[ -d medium/m${mem}/$YYYY ]] || mkdir -p medium/m${mem}/$YYYY 
			# Download
			mkdir download_m${mem} # make temp dir
			gsutil -qm cp ${base_url}/nwm.${YYYYMMDD}/medium_range_mem${mem}/nwm.t${HH}z.medium_range.channel_rt_${mem}.f{003..120..3}.conus.nc download_m${mem}
			wait
		
		    # process (concatenate; unpack; clip;)
			fout=download_m${mem}/${YYYYMMDD}${HH}.nc
			/home/kessler/.local/bin/cdo --no_warnings copy download_m${mem}/*.nc $fout 2> /dev/null 
			/usr/bin/ncpdq -O -U $fout $fout
			/usr/bin/ncks -O $dim_arg $var_arg $fout $fout

			# save final and cleanup
		    mv $fout medium/m${mem}/${YYYY}/${MMDD}${HH}.nc 
			chmod +r medium/m${mem}/${YYYY}/${MMDD}${HH}.nc 
			rm -rf download_m${mem}
		done
		echo ' '
	done
	echo ' ' 
	YYYYMMDD=$(date -d "$YYYYMMDD + 1 day" +%Y%m%d)
done

fi



echo finished at:  `date '+%D %T'` 
1>&2 echo finished at:  `date '+%D %T'` 

