# This script reads in the geogrid and plots the elevation
# of the file domain.

# Developer: Lindsay Fitzpatrick 
# Contact: ljob@umich.edu
# Version: 1
# Updated: 06-29-21

# import modules
import os, sys, zipfile, os.path
import netCDF4
from netCDF4 import Dataset
import numpy
from copy import copy

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

import wrf
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

import shapely.geometry as sgeom
import cartopy
from cartopy import crs
import cartopy.crs as ccrs
from cartopy import feature
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# -------------- inputs  --------------  #

# latitude
#lat_pnts = [44.018432,44.018432,44.751452,44.751452,47.326282,47.326282,44.55683,44.477225,44.477225,44.471538,44.471538,45.642755,45.035438,45.035438]
# longitude
#lon_pnts = [-92.708213,-92.708213,-94.373208,-94.373208,-96.151735,-96.151735,-94.801753,-95.264679,-95.264679,-95.269844,-95.269844,-94.943932,-94.058271,-94.058271]
#lat2_pnts = [41.016]
#lon2_pnts = [-84.44] 

infile = "/mnt/projects/hpc/fitzpatrick/Runoff/Retrospec/geo_em.d01.nc"

# longitude
xticks = list(numpy.arange(-90, -80, 0.5))
# latitude
yticks = list(numpy.arange(40, 45, 0.5))

# -------------------------------------  #

# all these functions are necessary only when LCC projection is used.
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


# Download and add the states and lakes
# https://www.naturalearthdata.com/
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
lakes = NaturalEarthFeature(category="physical", scale="50m",
                             facecolor="none",
                             name="lakes")


# Open geogrid (netCDF file)
ncfile = Dataset(infile, mode='r', format="NETCDF4")

# Get the elevation data
elevs = getvar(ncfile, "HGT_M")

# Get the latitude and longitude points
lats, lons = latlon_coords(elevs)

# Get the cartopy mapping object
cart_proj = get_cartopy(elevs)

# Create a figure
fig = plt.figure(figsize=(15,15))
# Set the GeoAxes to the projection used by WRF
ax = plt.axes(projection=cart_proj)

# add lists of latitutde & longitude points to map
#ax.plot([lon_pnts], [lat_pnts], color="red", marker="o",
#            transform=crs.PlateCarree())
#ax.plot([lon2_pnts], [lat2_pnts], color="blue", marker="o",
#            transform=crs.PlateCarree())

# add state & lake outlines
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.add_feature(lakes, linewidth=.5, edgecolor="black")

# Make the filled contours of the grid layer
plt.contourf(to_np(lons), to_np(lats), to_np(elevs), 60,
             transform=crs.PlateCarree(),
             cmap=get_cmap("terrain"))
# Add a color bar
plt.colorbar(ax=ax, shrink=.60, label="Elevation [m]")

# Set the map bounds
ax.set_xlim(cartopy_xlim(elevs))
ax.set_ylim(cartopy_ylim(elevs))

# Build lat/long grid and plot
fig.canvas.draw()
ax.gridlines(xlocs=xticks, ylocs=yticks)
ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
lambert_xticks(ax, xticks, 11)
lambert_yticks_left(ax, yticks, 11)

# Add title and sent font size.  Note the size
# of the title is used to set the size of the lat/lon
# tick marks labels.
plt.title("Elevation", size=14)
fig.savefig('elevation.png', bbox_inches='tight')
# Display the figure
#plt.show()

