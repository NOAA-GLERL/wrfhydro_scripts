#!/usr/bin/R -q
# this script reads in streamflow data from WRFHYDRO output and plots it on the stream network. 
# the CHRTOUT data must be concatenated into a single file with record dimension "time" [minutes since 1970-01-01]) 
# A spatial plot is generated for each desired time step so that an animation can be made.
# JAK Aug '19  (inspired by plot_wrfhydro_flow_vector_network.py)
library(rgdal,q=T)
library(rgeos,q=T) 
library(ncdf4,q=T)
library(chron)          # for datetime (non-essential)
library(colorspace,q=T) # for colorbar (non-essential)

# =================================================================================================
# user settings
# =================================================================================================

clvls <- c(-1,1,2,3,seq(4,10,2),seq(15,50,5),100,250,500,1000)        # colormap levels (needs to be tweaked based on domain, season, etc.)
#cmap <- c('oldlace',sequential_hcl(30,'sunset', rev=T))          # colormap    (best if first value is near white; dry streams barely visible)
cmap <- sequential_hcl(20,'blues', rev=T)          # colormap    (best if first value is near white; dry streams barely visible)
cmap[1] <- 'oldlace'
lk_file <- '/mnt/projects/ipemf/wrfhydro_data/champlain/shp/champlain_basin_100k_lakes.shp'      # shapefile for lakes
str_file <- '/mnt/projects/ipemf/wrfhydro_data/champlain/shp/champlain_basin_100k_flowlines.shp' # shapefile for stream network
flow_file <-'/mnt/projects/hpc/kessler/out/hydro/chrtout_2017.nc'                                # chrtout.nc file
plt_dir   <- './frames'                                                                          # directory to save output

t0 <- 3950
tf <- 4700
dt <- 1
# =================================================================================================
# read in files
# =================================================================================================

# shapefiles
lks <- readOGR(lk_file, v=F)
streams <- readOGR(str_file, v=F)
stream_ids <- unlist(streams@data['ID'])
if ( proj4string(streams) != proj4string(lks) ) { print('inconsistent projections!'); q(status=1) }

# netcdf
ncid <- nc_open(flow_file)
flow <- ncvar_get(ncid, 'streamflow')
flow_ids <- ncvar_get(ncid, 'feature_id')
order <- ncvar_get(ncid, 'order')[,1]
mins <- ncvar_get(ncid, 'time')
dtstr <- chron(mins/60/24, out.format=list(date='y mon d', time='h:m'))
dtstr <- gsub('[()]','',dtstr)

# =================================================================================================
# associate flows with streams
# =================================================================================================

# re-order flows based on IDs to match the streams (neccessary for plotting later)
flow <- flow[match(stream_ids, flow_ids),]
#order <- order[match(stream_ids, flow_ids)] # not use for now

# error checking:  reorder flow_ids & check that they all match now
flow_ids <- flow_ids[match(stream_ids, flow_ids)]         
if ( any(flow_ids != stream_ids) ) { print('reordering failed ~line 50'); q(status=1) }

streams_sl <- SpatialLines(streams@lines)          # grab lines
flow_df <- data.frame(flow)                        # grab data
#colnames(flow_df) <- sprintf('%04i',1:length(flow_df)) # rename columns to 000n format (
rownames(flow_df) <- as.character(stream_ids)      # match rownames (this isn't neccessary but good practice for later)
flow_spdf <- SpatialLinesDataFrame(streams_sl, data=flow_df,match.ID=F) # create SP DF 

# =================================================================================================
# make  plots
# =================================================================================================
for (t in seq(t0,tf,dt)){
	tstr = sprintf('%04i',t)
	print(tstr)
	png(sprintf('%s/%s.png', plt_dir, tstr), width=1200, height=1200)

	# define colorkey attributes
	ckey <- list(labels=as.character(clvls), col=cmap, at=seq(min(clvls),max(clvls),length.out=length(clvls)) )

	# draw lake 
	lk_plot <- list(lks, fill='darkseagreen', col='white', first=F)

	# build title
	my_title <- list(label=paste('Streamflow (CMS)',dtstr[t], sep=' : '), cex=2)

	# overlay streamflow
	print(spplot(flow_spdf, paste('X',t,sep=''), col.regions=cmap, colorkey=ckey, at=clvls, lwd=3, sp.layout=lk_plot, main=my_title))
	invisible(dev.off())
}
