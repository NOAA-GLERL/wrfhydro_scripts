# This script reads in precipitation data from the forcing files
# (LDASIN) and creates spatial plots.
# Adapted from a script written by Kevin Sampson (NCAR)

# Developer: Lindsay Fitzpatrick 
# Contact: ljob@umich.edu
# Version: 2
# Updated: 06-29-21

# import python core modules
import os
import sys
import time
from copy import copy

# import third-party modules
import numpy
import xarray as xr

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import wrf
from wrf import (getvar, get_cartopy, latlon_coords, projection)
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import shapely.geometry as sgeom
# --- Global Variables --- #

# infile is the geogrid for the geo information. infile2 is the actual data
inDir = '/glade/p/ral/hap/lread/GreatLakes/scratch_aorc_forcing/'

infile2 = os.path.join(inDir, '2007100921.LDASIN_DOMAIN1')
#infile2 = r'C:\Users\ksampson\Desktop\Lindsay_GLERL_Test\2021_03_16\model_output\2012061900.LDASIN_DOMAIN1'

# Output directory (to save image)
outDir = r'./'
outFName = 'rain_2007100921.png'

LDASIN = True
LDASOUT = False

# Variable to plot from infile2
varToPlot = "RAINRATE"
conversion = 3600 #convert from mm/s to in/hr
units = 'mm'
clim = [0, 100]                                 # Color-ramp limits. Essentially the data limits for vizualization
mask_vals = True                                # switch for masking a particular cell value
maskwhere = 0                                   # A value to set as transparent in the output

# Plot options
pTitle = "Precipitation"                        # Plot title
cmap = 'jet'                                    # Color bar: 'jet', 'terrain'
useNaturalEarth = True                          # Option to use downloadable Natural Earth data to supplement the political map.
NEscale = '50m'                                 # Scale of data to use from Natural Earth. Options = ['10m', '50m', '110m']
plotsize = (14, 10)
draw_grid_labels = True                         # Switch to draw graticule labels (slows down the script)

# Set the plot extent.
# Either specify one of the extents below or keep them commented to plot the entire domain
expand = False                                          # Plot an area larger than the grid. Overridden by specified extents below
#plot_extent_proj = [500000, 1000000, 200000, 700000]    # Specify plot limits in projected coordinates [Xmin, Xmax, Ymin, Ymax]
plot_extent_latlon = [-95, -85, 40, 50]                # Specify plot limits in geocentric coordinates [lon min, lon max, lat min, lat max]
# -a_ullr -54000.000 1281000.000 1947000.000 -320000.000

# Graticule
xlocs = range(-180,181,5)
ylocs = range(-90,91,5)


# --- Do not edit below this line --- #

# Projection parameters specific to NWM
projparams ={'CEN_LAT': 40.0,
             'CEN_LON': -97.0,
             'DX': 1000,
             'DY': 1000,
             'KNOWN_X': 0.0,
             'KNOWN_Y': 0.0,
             'MAP_PROJ': 1,
             'MOAD_CEN_LAT': 40.0,
             'POLE_LAT': 90.0,
             'POLE_LON': 0.0,
             'REF_LAT': 40.0,
             'REF_LON': -97.0,
             'REF_X': 0.0,
             'REF_Y': 0.0,
             'STAND_LON': -97.0,
             'TRUELAT1': 30.0,
             'TRUELAT2': 60.0}

# Switch for using wrf-python. Will not probably work with NWM outputs.
wrf_python_utils = False

# Option to flip the y-axis
flipY = True

# Dimension names to be used to identify certain known dimensions
tVarName = 'valid_time'                       # Time variable
yDims = ['south_north', 'y']
xDims = ['west_east', 'x']
timeDim = ['Time', 'time', 'valid_time']

interp_method = 'none'          # None, 'none', 'nearest', 'bilinear'
# --- End Global Variables --- #

# --- Functions --- #
def subset_ncVar(ncVar, dimensions, times=slice(None), flipIdx=-1):
    '''
    7/7/2020:
    This function will accept an array object, and will attempt
    to identify the time, x, and y dimensions. If requested, a dimension can be
    specified that will be reversed. This is typically "south_north" dimesnion.
    Also, a time index or slice may be provided. The time dimension size may
    only be greater than 1 if the variable has only 3 dimensions (time, x, y),
    for example.
    '''

    # Ensure x and y dimensions are in the dimensions of the input dataset
    assert all([any([dim in dimensions for dim in dims]) for dims in [yDims, xDims]])

    # Construct slice to index entire array in original order
    ind = [slice(None)] * len(dimensions)
    ind[flipIdx] = slice(None,None,-1)
    print("    Reversing order of dimension '{0}'".format(dimensions[flipIdx]))

    # Find the index for the time dimension
    timeDimIdx = [dimensions.index(dim) for dim in dimensions if dim in timeDim]
    if len(timeDimIdx) > 0:
        timeDimIdx = timeDimIdx[0]
        print("    Time dimension found: '{0}'.".format(dimensions[timeDimIdx]))
        ind[timeDimIdx] = times

        # Set any additional dimensions to 0
        additionals = [dim for dim in dimensions if dim not in set(yDims + xDims + timeDim)]
        for extraDim in additionals:
            print("    Selecting '{}' = 0.".format(extraDim))
            ind[dimensions.index(extraDim)] = 0
    else:
        print('    No time dimension found.')

    # Read the array as requested, reversing y if necessary, and subsetting in time
    print('    Dimensions and indices or slices on those dimensions:')
    for dim,indslice in zip(list(dimensions),ind):
        print('        {0}: {1}'.format(dim,indslice))
    ncArr = ncVar[ind]
    assert len(ncArr.shape) <= 3
    del dimensions, timeDimIdx, ind
    return ncArr

def find_corner(DX,DY,nrows,ncols,divFac,xOffset=0.0,yOffset=0.0):
    '''
    This function will take a grid (defined by rows and columns and cellsize) and
    produce the corner coordinates in projected space. Assumes that the 0,0 point
    is the center of the projected coordinate system. Allows for false eastings and
    northings (xOffset,yOffset) in projected coordinates.
    '''
    xdist = (float(DX)*float(ncols)/divFac)
    ydist = (float(abs(DY))*float(nrows)/divFac)
    minX = xOffset-xdist
    minY = yOffset-ydist
    maxX = xOffset+xdist
    maxY = yOffset+ydist
    return minX, maxX, minY, maxY

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

# --- End Functions --- #

# --- Main Codeblock --- #
if __name__ == '__main__':
    tic = time.time()
    print('Process initiated at {0}'.format(time.ctime()))

    ncfile2 = xr.open_mfdataset(infile2, concat_dim='Time', combine='nested')

    # Isolate variable of interest
    variables = ncfile2.variables
    ncVar = variables[varToPlot]
    varDims = ncVar.dims

    # Deal with time
    tVar = variables[tVarName]
    datevals = tVar.values
    datestr = ncfile2[tVarName].dt.strftime(r"%b %d %Y %X").values
    del tVar, datevals

    # Find the index for the y dimension
    xDimName = [dim for dim in varDims if dim in xDims][0]
    xDimIdx = varDims.index(xDimName)

    # Find the index for the y dimension
    yDimName = [dim for dim in varDims if dim in yDims][0]
    yDimIdx = varDims.index(yDimName)

    tDimName = [dim for dim in varDims if dim in timeDim][0]
    tDimIdx = varDims.index(tDimName)

    # Get arrays shape
    nc_cols = ncfile2.dims[xDimName]
    nc_rows = ncfile2.dims[yDimName]

    ncVarArr = ncfile2[varToPlot].values
    if flipY:
        ncVarArr = subset_ncVar(ncVarArr, varDims, flipIdx=yDimIdx)

    # This will sum over all times
    ncVarArr = numpy.sum(ncVarArr, axis=tDimIdx) * conversion

    # Use the wrf-python funcitonality to define the coordinate system given the input parameters
    inproj = projection.getproj(**projparams)
    cart_proj = inproj.cartopy()  # Create a cartopy projection object

    # Obtain grid resolution
    dx = float(projparams['DX'])
    dy = float(projparams['DY'])

    # Get the grid extent in projected coordinates.
    # Only works if you have a coordinate variables ['x', 'y']
    if LDASIN:
        minX = -54000.000
        maxY = 1281000.000
    if LDASOUT:
        minX = variables[xDimName][:].min() - (dx/2.)
        maxY = variables[yDimName][:].max() + (dy/2.)
    maxX = minX + (float(nc_cols)*dx)
    minY = maxY - (float(nc_rows)*dy)

    # Calculate the domain center as the midpoint of the x and y coordinate variable values
    d_center = [maxX-(maxX-minX)/2., maxY-(maxY-minY)/2.]
    expanded_extent = find_corner(dx, dy, nc_rows, nc_cols, 1.75,
        xOffset=d_center[0], yOffset=d_center[1])

    # Build the bounding box based on nest lower-left corner x,y and dx,dy and nrows,ncols
    if expand:
        xbbox = [minX, minX, maxX, maxX, minX]    # Used to draw a box in ax.plot()
        ybbox = [minY, maxY, maxY, minY, minY]    # Used to draw a box in ax.plot()
    domain_extent = [minX, maxX, minY, maxY]
    print('  Domain extent [minX, maxX, minY, maxY]: %s' %domain_extent)

    # Create a figure
    fig = plt.figure(figsize=plotsize)
    #ax = plt.axes(projection=cart_proj)                # Set the GeoAxes to the projection used by WRF
    ax = fig.add_subplot(1, 1, 1, projection=cart_proj)
    figsize = fig.get_size_inches()*fig.dpi # Figure size in pixels

    # To alter the extent of the plot, specify an axis limit in projected coordinates
    if expand:
        plot_extent = expanded_extent
        extentproj = cart_proj
    elif 'plot_extent_proj' in locals():
        plot_extent = plot_extent_proj
        extentproj = cart_proj
    elif 'plot_extent_latlon' in locals():
        plot_extent = plot_extent_latlon
        extentproj = ccrs.PlateCarree()
    else:
        plot_extent = domain_extent
        extentproj = cart_proj
    ax.set_extent(plot_extent, crs=extentproj)

    # Download and add the states and lakes: https://www.naturalearthdata.com/
    if useNaturalEarth:
        land = cartopy.feature.NaturalEarthFeature(
            category='physical',
            name='land',
            scale=NEscale,
            edgecolor='black',
            facecolor=cartopy.feature.COLORS['land'])           # edgecolor='face'
        land_outline = cartopy.feature.NaturalEarthFeature(
            category='physical',
            name='ocean',
            scale=NEscale,
            edgecolor='face',
            facecolor=cartopy.feature.COLORS['water'])
        lake_outline = cartopy.feature.NaturalEarthFeature(
            category='physical',
            name='lakes',
            scale=NEscale,
            facecolor=cartopy.feature.COLORS['water'])

        # Create a feature for States/Admin 1 regions from Natural Earth
        states_provinces = cartopy.feature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale=NEscale,
            facecolor='none')
        international = cartopy.feature.NaturalEarthFeature(
            category='cultural',
            name='admin_0_boundary_lines_land',
            scale=NEscale,
            facecolor='none')
        ax.add_feature(land_outline, zorder=0)
        ax.add_feature(states_provinces, edgecolor='grey', facecolor='none', linewidth=0.5, zorder=2)
        ax.add_feature(lake_outline, zorder=3, facecolor='grey')
        ax.add_feature(international, edgecolor='black', linewidth=1.0, zorder=4)
    else:
        # Using feature data
        ax.add_feature(cartopy.feature.BORDERS)

    # Make the filled contours of the grid layer
    if mask_vals:
        ncVarArr = numpy.ma.masked_where(ncVarArr == maskwhere, ncVarArr)
    plt.imshow(ncVarArr, cmap=cmap, extent=domain_extent, zorder=1, interpolation=interp_method)

    # Add the polygon boundary extent
    if expand:
        ax.plot(xbbox, ybbox, color='red', transform=cart_proj, zorder=5)

    # Add a color bar
    plt.colorbar(ax=ax, shrink=.65)
    plt.clim(clim[0], clim[1])

    # Grid Lines
    gl = ax.gridlines(color='black', linewidth=0.5, zorder=6, draw_labels=False, xlocs=xlocs, ylocs=ylocs) #, linestyle='--')
    gl.n_steps = 90                     # Smoothes the gridlines for curved domains like polar stereographic

    # Grid line labels
    if draw_grid_labels:
        fig.canvas.draw()
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        lambert_xticks(ax, list(xlocs), 14)
        lambert_yticks_left(ax, list(ylocs), 14)

    # Add title and set font size.
    plt.title(pTitle + '\n Variable: {0} Units: {1}\n {2}'.format(varToPlot, units, datestr[0]), size=14)
    #ax.text(plot_extent[1], plot_extent[3], 'Min={0: 3.2f} Max={1: 3.2f}'.format(ncVarArr.min(), ncVarArr.max()),
     #       horizontalalignment='right',
     #       verticalalignment='bottom',
     #       transform=cart_proj, size=12)
    #ax.text(plot_extent[0], plot_extent[3], '{0}'.format(datestr[0]),
     #       horizontalalignment='left',
     #       verticalalignment='bottom',
     #       transform=cart_proj, size=12)

    # Show the figure and save a copy when the displayed image is closed
    #plt.savefig(os.path.join(outDir, outFName), bbox_inches='tight')

    # Faster than the savefig method
    canvas = FigureCanvas(fig)
    canvas.print_figure(os.path.join(outDir, outFName),bbox_inches='tight')

    # Clean up
    ncfile2.close()
    print('Process completed in {0:3.2f}s'.format(time.time() - tic))
# --- End Main Codeblock --- #
