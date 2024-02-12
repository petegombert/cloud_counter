"""
This code is my first delve into the more complicated dataset.
Probably will eventually be used as a get data.


Leave notes: Need to make a fix for file doesn't exist.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from pyhdf.HDF import *
from pyhdf.VS import *
import netCDF4 as cdf
import cartopy.crs as ccrs
import cartopy.feature as feature
import julian
from datetime import datetime
import sys
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
import subprocess
from warnings import filterwarnings
import get_data
import mast_plot
import go_stats
import coarsify

#first lets look at cloud_column count & total column count.

data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_rlgeoprof_cbtz_dz2'

def find_files(data_id, yearly_flag):
    """
    This function pulls the files from the complex (non-zonal) dataset.
    Yearly flag = 1: Get all months of a specific year
    Yearly flag = 0: Get all years of a specific month
    Yearly flag = 2: Get specific month & specific year
    Returns:
        lat, lon, cbtz_bin, cbtz_dbin, dz_bin, dz_dbin, single_layer_occ, tot_count, year
    """
    print(str(data_id))
# The flags determine whether the data is grabbed by month or year
#    data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_Freq_RL'
    data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_rlgeoprof_cbtz_dz2'
    os.chdir(data_path)
    for filename in os.listdir(data_path):
        if filename.split('_')[-1][0:4] == str(data_id) and yearly_flag == 1:
            data_hold = grab_data(filename)
            try:
                data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3],data_hold[4],data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object)
        elif filename.split('_')[-1][4:6] == str(data_id) and yearly_flag == 0:
            data_hold = grab_data(filename)
            try:
                data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3],data_hold[4], data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object)
        elif filename.split('_')[-1][0:6] == str(data_id) and yearly_flag == 2:
            data_hold = grab_data(filename)
            try:
                data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3],data_hold[4], data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object)

    if yearly_flag == 'all':
        for filename in os.listdir(data_path):
            if filename.split('.')[-1] != 'txt':
                data_hold = grab_data(filename)
                try:
#                    print(filename)
                    data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object), axis=1)
                except UnboundLocalError:
#                    print(filename)
                    data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3], data_hold[4], data_hold[5], data_hold[6], data_hold[7], data_hold[8]], dtype=object)

    try:
        data_arr = np.rot90(data_arr)
    except UnboundLocalError:
        return 0
#    data_arr = data_arr[data_arr[:,0].argsort()]
    return data_arr

def grab_data(filename):
    os.chdir(data_path)
    year = filename.split('_')[-1][0:4]
    print(filename)
    f = cdf.Dataset(filename)
    lat_get = f.variables['Latitude']
    lon_get = f.variables['Longitude']
    cbtz_bin_get = f.variables['cbtz_bin_mid']
    cbtz_dbin_get = f.variables['cbtz_dbin']
    dz_bin_get = f.variables['dz_bin_mid']
    dz_dbin_get = f.variables['dz_dbin']
    single_layer_occ_get = f.variables['cbz_dz1_grid']
#    col_count_get = f.variables['cld_column_count']
    tot_count_get = f.variables['total_column_count']

    lat = lat_get[:]
    lon = lon_get[:]
    cbtz_bin = cbtz_bin_get[:]
    cbtz_dbin = cbtz_dbin_get[:]
    dz_bin = dz_bin_get[:]
    dz_dbin = dz_dbin_get[:]
    single_layer_occ = single_layer_occ_get[:]
#    col_count = col_count_get[:]
    tot_count = tot_count_get[:]

    data_arr = np.array([[lat], [lon], [cbtz_bin], [cbtz_dbin], [dz_bin], [dz_dbin], [single_layer_occ], [tot_count], [year]], dtype=object)
    return data_arr


#data = grab_data('CPR_Monthly_rlgeoprof_cbtz_dz2_200904.cdf')
#print(data[1,0].shape)


"""
lat = data[0]
lon = data[1]
col_count = data[2]
col_count = col_count[1,:,:]
col_count = np.rot90(col_count)
col_count_coarse = coarsify.make_me_coarse_2d(col_count, lat, lon)

plt.imshow(col_count)
plt.show()
plt.imshow(col_count_coarse)
plt.show()
print(col_count.shape)
print(col_count_coarse.shape)
exit()
tot_count = data[3]
tot_count = tot_count[1,:,:]
tot_count = np.rot90(tot_count)

proj = ccrs.PlateCarree()
fig, ax = plt.subplots(1, 2, figsize=(10,10))

#ax1 = plt.subplot(1, 2, 1, projection=proj)
ax1 = plt.subplot2grid(shape=(3,3), loc=(0,0), colspan = 2, rowspan=2, projection=proj)
g1 = ax1.gridlines(draw_labels = False, linestyle = 'dotted')
g1.top_labels = True
g1.right_labels = True
ax1.coastlines()

p1 = ax1.imshow(col_count, origin='lower', aspect='equal', extent=[-180,180,-90,90])
cb1 = plt.colorbar(p1, shrink=0.5, ax=ax1)

ax2 = plt.subplot(1, 2, 2, projection=proj)
g2 = ax2.gridlines(draw_labels = False, linestyle = 'dotted')
g2.top_labels = True
g2.right_labels = True
ax2.coastlines()

p2 = ax2.imshow(col_count_coarse, origin='lower', aspect='equal', extent=[-180,180,-90,90])
cb2 = plt.colorbar(p1, shrink=0.5, ax=ax2)

plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/new_data_climate_test.png')
plt.close()

"""
