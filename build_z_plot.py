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

year = 2008
daytime_filter = 0

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

    data_arr = np.rot90(data_arr)
    data_arr = data_arr[data_arr[:,0].argsort()]

    return data_arr

def sum_along_dbz(data):
    cloud_occ = data #assumes we always want daytime values. I think this is the case.
    tstep = np.arange(0, cloud_occ.shape[0])
    for t in tstep:
        if t == 0:
            sum = cloud_occ[0, :, :, :]
            continue
        hold = cloud_occ[t, :, :, :]
        sum = sum+hold

    return sum

def sum_along_lon(data):
    cloud_occ = data
    tstep = np.arange(0, cloud_occ.shape[1])
    for t in tstep:
        if t == 0:
            sum = cloud_occ[:, 0, :]
            continue
        hold = cloud_occ[:, t, :]
        sum = sum+hold

    return sum

def make_plot():
    data = lon_sum
    print(data.shape)
    z = get_files_dict[0][3]
    lat = get_files_dict[0][1]

    fig, ax = plt.subplots(1, 1, figsize = (10, 10))

    ax.contourf(lat, z, data)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/climate_height.png')

get_files_dict = get_files()
data_send = get_files_dict[0][-1]
data_send = np.rot90(data_send)
dbz_sum = sum_along_dbz(data_send[:,1,:,:,:])
lon_sum = sum_along_lon(dbz_sum)
make_plot()

