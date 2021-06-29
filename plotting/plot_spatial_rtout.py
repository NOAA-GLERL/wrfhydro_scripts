# This script reads in surface runoff from the output files
# (RTOUT) and creates spatial plots.
# Adapted from a script written by Kevin Sampson (NCAR)

# Developer: Lindsay Fitzpatrick 
# Contact: ljob@umich.edu
# Version: 2
# Updated: 06-29-21

import netCDF4
from netCDF4 import Dataset
import numpy

import pyproj
import xarray as xr
from copy import copy

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

import shapely.geometry as sgeom
import cartopy
from cartopy.feature import NaturalEarthFeature
from cartopy import crs
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from netCDF4 import Dataset, num2date
from datetime import datetime


# -------------- inputs  --------------  #
lat = 47.5
lon = -91.5


print('Reading in files')

dir = '/mnt/projects/hpc/fitzpatrick/Runoff/Retrospec/'
infile = dir + 'DOMAIN/Fulldom_hires.nc'
infile2 = dir + '201206200000.RTOUT_DOMAIN1_DailySum.nc'

ingeogrid = dir + 'DOMAIN/geo_em.d01.nc'


ncfile = Dataset(infile, mode='r', format="NETCDF4")
ncfile2 = Dataset(infile2, mode='r', format="NETCDF4")

# Define grid cell size of input file
cellSize = 250

# Set lat / lon ticks
# longitude
xticks = list(numpy.arange(-90, -80, 0.5))
# latitude
yticks = list(numpy.arange(40, 45, 0.5))



# -------------------------------------  #


def getProjection(geogrid) :  
    # grab projection information from geogrid file
    ds = xr.open_dataset(geogrid)
    
    lccproj = pyproj.Proj(proj='lcc',                                  # projection type: Lambert Conformal Conic
                           lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2,        # Cone intersects with the sphere
                           lat_0=ds.MOAD_CEN_LAT, lon_0=ds.STAND_LON,   # Center point
                           a=6370000, b=6370000)                        # Earth is a perfect sphere
    
    wgs84proj = pyproj.Proj("+init=EPSG:4326")
    
    return lccproj, wgs84proj

def getCornerPnt(routegrid) :

    nc = Dataset(routegrid, mode='r', format="NETCDF4")

    # grab the starting corner point Fulldom_hires script
    # note this is the center point of the 250 m grid cell
    # minimum x value
    x_vals = nc.variables['x']
    xll = x_vals[:].min() - cellSize / 2.0

    # maximum y value
    y_vals = nc.variables['y']
    yul = y_vals[:].max() + cellSize / 2.0

    return xll, yul

def getij (x, y) :
    xp = int((xx - xll) / cellSize )
    yp = int(((yy - yul)) * -1 / cellSize )

    return xp, yp



def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])

def lambert_xticks(ax, ticks, size):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: numpy.vstack((numpy.zeros(n) + t, numpy.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels], size=size)
    
def lambert_yticks_left(ax, ticks, size):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: numpy.vstack((numpy.linspace(b[0], b[1], n), numpy.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels], size=size)

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels


# Get projection of data from geogrid &
# lat/lon projection information
wrf_proj, wgs84_proj = getProjection(ingeogrid)

# convert lat/lon input to x,y coordinate values
xx, yy = pyproj.transform(wgs84_proj, wrf_proj, lon, lat)

# get the starting corner of the routing grid
# which is the upper left corner
xll, yul = getCornerPnt(infile)
#print(xll, yul)

# Find i,j value on routing grid
i,j = getij(xx, yy)
#print(i,j)

print("running")

# https://joehamman.com/2013/10/12/plotting-netCDF-data-with-Python/
lons = ncfile.variables['LONGITUDE'][:]
lats = ncfile.variables['LATITUDE'][:]
topo = ncfile.variables['TOPOGRAPHY'][:]

topo_units = ncfile.variables['TOPOGRAPHY'].units

runoff = ncfile2.variables['qqsfc_acc'][:]

ds = xr.open_dataset(ingeogrid)

# Define the projection
globe = ccrs.Globe(ellipse='sphere', semimajor_axis=6370000, semiminor_axis=6370000)
lcc = ccrs.LambertConformal(globe=globe, # important!
                            central_longitude=ds.STAND_LON, central_latitude=ds.MOAD_CEN_LAT,
                            standard_parallels=(ds.TRUELAT1, ds.TRUELAT2),
                            )

# Download and add the states and lakes
# https://www.naturalearthdata.com/
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
lakes = NaturalEarthFeature(category="physical", scale="50m",
                             facecolor="none",
                             name="lakes")

# WARNING THIS IS TIME CONSUMING!
# matplotlib grabs all the center cells of each
# and makes a contour map. The more grid cells
# the longer it takes to run.

print('Plotting')
# Create a figure
fig = plt.figure(figsize=(10, 10))
# Set the GeoAxes to the projection used by WRF
ax = plt.axes(projection=lcc)

#ax.plot(xx, yy, color="red", marker="o", markersize=8)
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.add_feature(lakes, linewidth=.5, edgecolor="black")

# Make the filled contours 
levels = [x for x in range(0,1000,1)]
plt.contourf(lons, lats, runoff, transform=crs.PlateCarree(), levels=levels, cmap=get_cmap("jet"))
	
# Add a color bar
plt.colorbar(ax=ax, shrink=.70)
#plt.clim(0,0.5)
	
# Build lat/long grid and plot
fig.canvas.draw()
ax.gridlines(xlocs=xticks, ylocs=yticks)
ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
lambert_xticks(ax, xticks, 14)
lambert_yticks_left(ax, yticks, 14)
# Add title and sent font size.  Note the size
# of the title is used to set the size of the lat/lon
# tick marks labels.
#plt.title("Elevation (m)", size=14)
title = "Surface Runoff"
plt.title(title, size=14)

print('Saving') 
	
fig.savefig('runoff.png',bbox_inches='tight')
plt.close(fig)
