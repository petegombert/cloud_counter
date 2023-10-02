#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

* This script will build anomaly data/plots for each month from the average.
leave notes: getting data from the import run. Need to figure out how to return the data back here.

To do:
	Need to figure out how to aggrigate multiple months. Should just be able to add them up, not that clean though...
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

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

year_start = 2006
daytime_filter = 0
min_dbz = -39.

win_months = ['12', '01', '02']
som_months = ['06', '07', '08']
seas_sel = som_months

def get_monthly_data(month):
    os.chdir('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj')
    print('fetching data')
    data = get_data.find_files(month, 3)
    get_data.test_monthly_data(data)
    z = data[0,3]
    lat = data[0,1]
    dbz = data[0,4]
    data_send = go_stats.sum_months(data)

    return np.array([data_send, z, lat, dbz], dtype=object)

def do_stats():
    data = get_monthly_data_dict
#calls go_stats.py to build cloud_occ stats
    cld_occ = go_stats.build_stats(data, 0)

    return cld_occ

def make_cdf():
    ncdf_save_path = '/uufs/chpc.utah.edu/common/home/u1113223/grad_research/ncdfs'
    os.chdir(ncdf_save_path)
    ncfile = cdf.Dataset('jja_avg_lt2017.cdf', mode='w', format='NETCDF4')
    height = ncfile.createDimension('height')
    latitude = ncfile.createDimension('latitude')
    total_mean = ncfile.createVariable('total_mean', np.float32, ('height','latitude'))
    total_mean[:] = cld_occ[:]

if __name__ == '__main__':
    for t in seas_sel:
        print('starting month '+str(t))
        if t =='06' or t=='12':
            get_monthly_data_dict = get_monthly_data(t)
            continue
        get_monthly_data_dict[0] = get_monthly_data_dict[0]+get_monthly_data(t)[0]
    print('making cloud occurrence stats')
    cld_occ = go_stats.build_stats_presum(get_monthly_data_dict[0], get_monthly_data_dict[3], get_monthly_data_dict[2], 1)
    print('making ncdf')
    make_cdf()
    mast_plot.make_seas_avg_plot(cld_occ, 'JJA')
