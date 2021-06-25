#!/bin/Rscript
# read in WRFHydro output and then: a) draw hydrographs and/or b) compute skill stats
# contact: james.kessler@noaa.gov 

library(ncdf4)

graphics.off()
rm(list=ls())
load('OBS.Rdata')    # load OBS (previously saved by `get_obs.R`)
source('myutils2.R') # load my functions

# ==========  controls ========================

fn <- 'base_CHANOBS.nc'
#riv_names <- read.table('names.txt',row.names=1) # optional file with river names

plt <- T   # draw hydrographs
prnt <- F  # calculate stats (correlation, bias, NSE, total obserered streamflow)

# define stats output option:
#fout='stats/base_big.txt' # write stats to txt file
fout=''                    # write to stdout

# =============================================


if (prnt) { 
}



ncid <- nc_open(fn)
dts <- as.POSIXct(ncvar_get(ncid,'time')*60, origin='1970-01-01', tz='Z')
fids <- ncvar_get(ncid,'feature_id')
gagetab <- read.table('gages.txt', col.names=c('fid','gageid'),colClasses=c('numeric','character'))  # output from routelink 

#x11(w=15,h=12)
figi <- 1 # figure index
png(sprintf('%i.png',figi), w=2600,h=1600, pointsize=20)
par(mar=c(0,1,1,0), oma=c(5,4,4,2),lwd=2)
layout(matrix(1:20,4,5) )

for (i in 1:nrow(gagetab)){
	fid <- gagetab[i,1] # NWM Feature ID
	gid <- gagetab[i,2] # USGS Gage ID
	# check for obs availability: if none, advance
	if (sum(OBS$site_no==gid) == 0) next

	# read in hydro data at desired stream segment
	sim <- read_hyd(fid)

	# select obs
	obs <- {}
	obs$dts <- OBS[OBS$site_no==gid,]$dateTime
	obs$flo  <- OBS[OBS$site_no==gid,]$X_00060_00000/35.31


	if (plt){ 
		#x11()
		plot(dts, sim, 'l', col='blue', ylab=NA, xlab=NA, axes=F, ylim=range(obs$flo,sim, finite=T))
		axis(2, pretty(c(sim,obs$flow)), lwd=.2)
		lines(obs$dts, obs$flo, lty=3)
		legend('topright', col=c('black', 'blue'), lty=c(3,1), legend=c(gid, fid), inset=c(0.1,0.2))

		# handle the pannels
		if(par()$mfg[1]==c(4)) axis.POSIXct(1,dts) # date axis for bottom row
		if(sum(par()$mfg)==18) {  #panel is full; setup new figure
			dev.off()
			figi <- figi+1
			png(sprintf('%i.png',figi), w=2600,h=1600, pointsize=20)
			par(mar=c(0,1,1,0), oma=c(5,4,4,2),lwd=2); 
			layout(matrix(1:20,4,5)) }
	}

	
	if (prnt){ 
		# print header
		if (i==1) cat('feature\t\tcor\tbias\tnse\ttotal_flow_cms\n', file=fout); 

		# intersect obs and model using dts
		obs_int <- obs$flo[match(dts, obs$dts, nomatch=0)]
		sim_int <- sim[match(obs$dts, dts, nomatch=0)]
		run_stats(fid, obs_int, sim_int)	
	}
}

dev.off()
