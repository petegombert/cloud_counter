#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

*This script is designed to recreate Jay's idl code that plots longitutally summed data over the year and plots cloud occurence frequencies for the year. Based on z_no_sum.py.
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
min_dbz = -3

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

def test_monthly_data():
#tests to makes sure that the dbz/z bins are the same. Summation will not work if they aren't all the same.
    data = get_files_dict
    monthly_dbz = data[:,4]
    monthly_z = data[:,3]
    tstep = np.arange(0, monthly_dbz.shape[0])
    for t in tstep:
        if t == 0:
            hold_dbz = monthly_dbz[t]
            hold_z = monthly_z[t]
            continue
        test_arr_dbz = monthly_dbz[t]-hold_dbz
        test_arr_z = monthly_z[t]-hold_z
        if np.sum(test_arr_dbz) != 0:
            print('dbz values are not the same')
            exit()

        if np.sum(test_arr_z) != 0:
            print('z values are not the same')
            exit()
    print('dbz/z values are the same between months')

def sum_months():
    data = get_files_dict
    tstep = np.arange(0, data.shape[0])
    for t in tstep:
        if t == 0:
            mast_arr = data[t][-1]
            continue
        hold = data[t][-1]
        mast_arr = mast_arr+hold

    return mast_arr

def sum_along_dbz(data):
#Month_sum will only be okay if we have the same values for z and dbz for all months. Should probably make it error if this is not the case.
    cloud_occ = data #assumes we always want daytime values. I think this is the case.
    dbz = get_files_dict[0][4]
    if min_dbz != 1000:
        start_index = np.where(dbz==min_dbz)
    else:
        start_index = np.where(dbz==0) #finds where dbz equals 0, this will be the starting spot for the summation.
    tstep = np.arange(int(start_index[0]), dbz.shape[0])
    for t in tstep:
        if t == int(start_index[0]):
            sum = cloud_occ[t, :, :]
            continue
        hold = cloud_occ[t, :, :]
        sum = sum+hold

    return sum

def make_plot2():
    data = data_send
    data = np.sum(data, axis=2)
    data = sum_along_dbz(data)

    z = get_files_dict[0][3]
    lat = get_files_dict[0][1]

    fig, ax = plt.subplots(1, 1, figsize = (10, 10))

    ax.imshow(data)
    ax.set_ylim([0, 15])
    ax.set_ylabel('Height (km)')
    ax.set_xlabel('Latitude')
    ax.set_title('200801-200812 Mean Cloudiness. No min dBZ')
#    ax=plt.gca() #get the current axes
#    PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
#    plt.colorbar(PCM, ax=ax)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/climate_height_sum_min_3.png')

get_files_dict = get_files()
test_monthly_data()
month_sum_data = sum_months()
data_send = month_sum_data[1,:,:,:,:]
make_plot2()

