#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 13:28:42 2021

@author: donglaiyang
    Analyzing ISMIP6 data
    Getting to use NetCDF4 and Xarray
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import cftime
import nc_time_axis

datadir = '../AWI/ISSM1/historical/'
AllFiles = os.listdir(datadir)

# using xarray to import the dataset
ds = xr.open_dataset(datadir + AllFiles[0],
                       engine='netcdf4')

ds.data_vars
ds_use = ds['yvelmean']
SLS = ds_use.isel(x = 150, y = 150)
SLS.plot()


fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(14,4))
ds_use.lon.plot(ax=ax1)
ds_use.lat.plot(ax=ax2)

# visualizing on the logical grid (x, y)
plt.figure()
ds.yvelmean[0].plot()

# visualizing over the projection grid
plt.figure(figsize=(14,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ds.yvelmean[0].plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), x='lon', y='lat', add_colorbar=False)
ax.coastlines()
ax.set_ylim([0,90]);



# choose a projection
proj=ccrs.NorthPolarStereo()


### a number of common comands
ds.info()
ds.coords # 
ds.dims # here, time is dimension coordinates, lon and lat are not (non-dim coords)

# output ndarrays
ds_use.lon.values
ds_use.lat.values

