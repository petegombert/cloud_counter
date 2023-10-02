#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

* This script gets the monthly data from the cdfs, will be called from other programs.
leave notes:

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
import subprocess
from warnings import filterwarnings

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

# data_id is either year or month depending on what the yearly_flag is set to.
# if yearly_flag = 0, it will get monthly data, if it is set to 1, it will get yearly data.
try:
    data_id = sys.argv[1]
except IndexError:
    data_id = '01'
try:
    yearly_flag = int(sys.argv[2])
except IndexError:
    yearly_flag = 0

def find_files(data_id, yearly_flag):
# The flags determine whether the data is grabbed by month or year
    data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_Freq_RL'
    os.chdir(data_path)
    for filename in os.listdir(data_path):
        if filename.split('_')[-1][0:4] == str(data_id) and yearly_flag == 1:
            data_hold = grab_data(filename)
            try:
                data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3],data_hold[4],data_hold[5]], dtype=object)
        elif filename.split('_')[-1][4:6] == str(data_id) and yearly_flag == 0:
            print(filename)
            data_hold = grab_data(filename)
            try:
                data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3],data_hold[4],data_hold[5]], dtype=object)
        elif filename.split('_')[-1][0:6] == str(data_id) and yearly_flag == 2:
            print(filename)
            data_hold = grab_data(filename)
            try:
                data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3],data_hold[4],data_hold[5]], dtype=object)
        #This is for grabbing months that correspond with the global dataset (exclude 2011 and years later than 2016.)
        elif filename.split('_')[-1][4:6] == str(data_id) and yearly_flag ==3 and filename.split('_')[-1][0:4] != '2011' and int(filename.split('_')[-1][0:4]) <= 2016:
            data_hold = grab_data(filename)
            if data_id == '06':
                print(filename.split('_')[-1][0:4])
            try:
                data_arr = np.append(data_arr, np.array([data_hold[0], data_hold[1], data_hold[2], data_hold[3], data_hold[4], data_hold[5]], dtype=object), axis=1)
            except UnboundLocalError:
                data_arr = np.array([data_hold[0],data_hold[1],data_hold[2],data_hold[3],data_hold[4],data_hold[5]], dtype=object)

    try:
        data_arr = np.rot90(data_arr)
        data_arr = data_arr[data_arr[:,0].argsort()]
        return data_arr
    except UnboundLocalError:
        return 0

def grab_data(filename):
    data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_Freq_RL'
    os.chdir(data_path)
    month = filename.split('.')[0][-2:]
    year = filename.split('_')[-1][0:4]

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

    if yearly_flag == 1:
        return np.array([[month],[lat],[lon],[z],[dbz],[cld_occ]], dtype=object)
    else:
        return np.array([[year],[lat],[lon],[z],[dbz],[cld_occ]], dtype=object)

def test_monthly_data(data_send):
#tests to makes sure that the dbz/z bins are the same. Summation will not work if they aren't all the same.
    data = data_send
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

#* Array will be (month, lat, lon, height, dbz, cld_occ)
#* Array dictionary: if yearly_flag = 1: [month,data], if yearly_flag = 0: [year, data]
#* Daytime flag will be excluded, -1 if not daytime, 1 if daytime.

if __name__ == '__main__':
    data = find_files(data_id, yearly_flag)
