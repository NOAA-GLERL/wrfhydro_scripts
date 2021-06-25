# useful fxns for working with wrfhydro output
# contact: james.kessler@noaa.gov


read_nwm <- function(fin,my_fids){

	        # handle feature_ids IN NETCDF
	        ncid <- nc_open(fin)
	        f_fids <- ncvar_get(ncid,'feature_id') # fids in chanobs file
	        f_sel <- match(my_fids, f_fids) # the indices of the gaged fids in the files
	        f_sel <- f_sel[!is.na(f_sel)] #omit NA's

		#       order <- ncvar_get(ncid,'order')[f_sel] # this is missing from the clipped files but would be nice to have
		        flow <- ncvar_get(ncid, 'streamflow')[f_sel,]
		        rownames(flow) <- f_fids[f_sel]
			        dts <- ncvar_get(ncid,'time')  # this could result in mismatched dates (careful!)
			        dts <- as.POSIXct(dts*60, origin='1970-01-01', 'UTC')
				        return(list(flow=flow, dts=dts))
}



read_hyd <- function(fid){ # read data for a single id (assumes ncid, fids are defined)
	sel <- which(fids==fid)
	ncvar_get(ncid, 'streamflow', start=c(sel,1), count=c(1,-1))
}




start_end_points <- function(spdf){
	sp <- sapply(coordinates(spdf), function(l) l[[1]][1,])
	ep <- sapply(coordinates(spdf), function(l) l[[1]][dim(l[[1]])[1],])
	return(list(startp=sp,endp=ep))
}



bias <- function(obs, sim){ mean(obs-sim) }
nse <- function(obs, sim){ 
	#1 -     SSE        / variance of obs * N
	1- sum((obs-sim)^2)/sum((obs-mean(obs))^2)
}

run_stats<- function(fid, obs, sim){
	
	r <- cor(obs,sim); 
	bi <- bias(obs,sim)
	nash <- nse(obs,sim) 
	sigflo <- sum(obs) # total streamflow during period 
	cat(sprintf('%12s\t%4.2f\t%5.1f\t%5.2f\t%9.0f\n', fid, r, bi, nash, sigflo), file=fout, append=T)
}
