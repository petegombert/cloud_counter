#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

*This is the first script of the new climate project. This will grab a year of data and plot cloudyness fraction averages for the year.
leave notes: Not sure what cbz_dz or ctz_dz or cbtz is. Might have to ask Jay or look how he uses in IDL.
             	--It appears that there are 2 data arrays. One that the statistics are generated from and another that show global mean statistics. Not sure where the first is going to
			come into play.
	     Data array is going to be very confusing. It's a 5-6 dimensional array. Will need to really think through this one and create a really nice key for myself.
To do:
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

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

year = 2010

def get_files():
#    data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_rlgeoprof_cbtz_dz2/'
    data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_Freq_RL'
    os.chdir(data_path)
    for filename in os.listdir(data_path):
        if filename.split('_')[-1][0:4] == str(year):
            month = filename.split('.')[0][-2:]

            f = cdf.Dataset(filename)
            lat_get = f.variables['Latitude']
            lon_get = f.variables['Longitude']
            z_get = f.variables['Height_MSL']
            dbz_get = f.variables['DbZ']
            dt_flag_get = f.variables['Daytime_flag']
            cld_occ_get = f.variables['Cld_Occurrence']

            lat = np.array(lat_get[:])
            lon = np.array(lon_get[:])
            z = z_get[:]
            dbz = dbz_get[:]
            dt_flag = dt_flag_get[:]
            cld_occ = cld_occ_get[:]

#* Array will be (month, lat, lon, height, dbz, cld_occ)
#* Array dictionary: [month,data]
#* Daytime flag will be excluded, -1 if not daytime, 1 if daytime.

            try:
                data_arr = np.append(data_arr, np.array([[month], [lat], [lon], [z], [dbz], [cld_occ]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([[month],[lat],[lon],[z],[dbz],[cld_occ]], dtype=object)

    print(data_arr.shape)
    data_arr = np.rot90(data_arr)
    data_arr = data_arr[data_arr[:,0].argsort()]
    print(data_arr.shape)
    exit()
    return data_arr

def sum_along_axis(month_index, axis, axis_hold):
    data_sel = get_files_dict[month_index,:]
    cloud_occ = data_sel[-1]
    cloud_occ = cloud_occ[1,:] #assumes we always want daytime values. I think this is the case.
    tstep = np.arange(0, cloud_occ.shape[axis-1])
    for t in tstep:
        if axis == 1:
            if t == 0:
                sum = cloud_occ[0, axis_hold, :, :]
                continue
            hold = cloud_occ[t, axis_hold, :, :]
            sum = sum+hold
        else:
            if t == 0:
                sum = cloud_occ[axis_hold, 0, :, :]
                continue
            hold = cloud_occ[axis_hold, t, :, :]
            sum = sum+hold

    return sum

def test_plot():
    data_sel = get_files_dict[0,:]
    month = data_sel[0]
    month_dt = datetime.strptime(str(month), '%m')
    month_name = datetime.strftime(month_dt, '%b')
    lat = data_sel[1]
    lon = data_sel[2]
    z = data_sel[3]
    dbz = data_sel[4]
#    cloud_occ = data_sel[-1]
#    cloud_occ = cloud_occ[1, 2, 3, :, :]
    cloud_occ = z_sum
    cloud_occ = np.rot90(cloud_occ)

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(1, 1, figsize = (10, 10))

    axs = plt.subplot(1, 1, 1, projection=proj)
    axs.coastlines()
    gs = axs.gridlines(draw_labels = True, linestyle = 'dotted')
    gs.top_labels = False
    gs.right_labels = False

#    axs.contourf(lon, lat, cloud_occ)
    axs.imshow(cloud_occ, transform=ccrs.PlateCarree())
    axs.set_title('Monthly Cloud Occurance: ' + month_name + ', ' + str(year) + '\n Height: ' + str(z[3]) + '. dBZ: ' + str(dbz[2]))

#    ax=plt.gca() #get the current axes
#    PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
#    plt.colorbar(PCM, ax=ax)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/climate_test_zsum_imshow.png')

get_files_dict = get_files()
dbz_sum = sum_along_axis(0, 1, 3)
z_sum = sum_along_axis(0, 2, 3)
test_plot()
