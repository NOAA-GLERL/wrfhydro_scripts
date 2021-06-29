# This script reads in individual hourly files from WRF-Hydro and 
# aggregates them to daily values.

# Developer: Lindsay Fitzpatrick 
# Contact: ljob@umich.edu
# Version: 1
# Updated: 05-17-21

# import python core modules
import os
import sys
import time
from copy import copy

# import third-party modules
import numpy as np
import xarray as xr

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

### Cutout the GL domain ###

## Directory to the hourly output ##
# These are the 1km land files
inland = '/glade/scratch/lfitzpatrick/Scripts/Tmp/*.LDASOUT_DOMAIN1'
# These are the 250m routing files
inroute = '/glade/scratch/lfitzpatrick/Scripts/Tmp/*.RTOUT_DOMAIN1'

# Track run time 
tic = time.time()

# Reads in all the files and makes one long record according to 'time'
ncland = xr.open_mfdataset(inland, combine='by_coords', parallel=True)
variables = ncland.variables
ncroute = xr.open_mfdataset(inroute, combine='by_coords', parallel=True)

# Define known dimensions
tVarName = 'time'                       # Time variable
yDims = ['south_north', 'y']
xDims = ['west_east', 'x']
timeDim = ['Time', 'time', 'valid_time']

# Deal with time
tVar = variables[tVarName]
datevals = tVar.values
datestr = ncland[tVarName].dt.strftime(r"%b %d %Y %X").values
del tVar

# Aggregate entire dataset to daily values
## variables that take daily accumulations (sum)
SWFORC = ncland.SWFORC.resample(time='D').sum() # W/m2
LWFORC = ncland.LWFORC.resample(time='D').sum() # W/m2
FSA = ncland.FSA.resample(time='D').sum() # W/m2
FIRA = ncland.FIRA.resample(time='D').sum() # W/m2
IRB = ncland.IRB.resample(time='D').sum() # W/m2
QQSFC = ncroute.qqsfc_acc.resample(time='D').sum() # mm
QQSUB = ncroute.qqsub_acc.resample(time='D').sum() # mm

## variables that take the daily mean
TGB = ncland.TGB.resample(time='D').mean() # K
T2MV = ncland.T2MV.resample(time='D').mean() # K
Q2MV = ncland.Q2MV.resample(time='D').mean() # kg/kg
SOIL_M = ncroute.SOIL_M.resample(time='D').mean() # m3/m3
SOIL_T = ncland.SOIL_T.resample(time='D').mean() # K
SOIL_W = ncland.SOIL_W.resample(time='D').mean() # m3/m3
SOILICE = ncland.SOILICE.resample(time='D').mean() # fraction
SOILSAT = ncland.SOILSAT.resample(time='D').mean() # fraction
ZWATTABLRT = ncroute.zwattablrt.resample(time='D').mean() # m

## variables that are accumulated where the daily must be calculated
## by taking the last value and subtracting the first
PRECIP = ncland.ACCPRCP.resample(time='D',closed='right').last(skipna=False)-ncland.ACCPRCP.resample(time='D').first(skipna=False) # mm
ACSNOM = ncland.ACSNOM.resample(time='D',closed='right').last(skipna=False)-ncland.ACSNOM.resample(time='D').first(skipna=False) # mm
QBDRYRT = ncroute.QBDRYRT.resample(time='D',closed='right').last(skipna=False)-ncroute.QBDRYRT.resample(time='D').first(skipna=False) # mm

## variables that take the mean of all non zero hourly values ##
EMISS = ncland['EMISS']
# this step replaces 0 with nans, takes the mean, turns nans back to 0 
EMISS = EMISS.where(EMISS != 0).resample(time='D').mean(skipna=True).fillna(0)
SFCHEADSUBRT = ncroute['sfcheadsubrt']
SFCHEADSUBRT = SFCHEADSUBRT.where(SFCHEADSUBRT != 0).resample(time='D').mean(skipna=True).fillna(0) # mm