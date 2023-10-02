#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

* This script will build anomaly data/plots for each month from the average.
leave notes: getting data from the import run. Need to figure out how to return the data back here.

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
import get_data
import mast_plot
import go_stats
import build_ann_avg

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

years = [2006, 2007, 2008, 2009, 2010, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
daytime_filter = 0
min_dbz = -39.

def pull_data(month, year):
    data_path = '/uufs/chpc.utah.edu/common/home/mace-group4/qzhang/occ/P1_R05/CPR_Monthly_Freq_RL'
    os.chdir(data_path)
    for filename in os.listdir(data_path):
        if filename.split('_')[-1][0:4] == str(year) and filename.split('_')[-1][4:6] == str(month):
            data = get_data.grab_data(filename)
            break
    try:
        return data
    except NameError:
        print('No Data')
        return 'no data'
#        exit()

if __name__ == '__main__':
    for year in years:
        print('starting year: '+str(year))
        try:
            os.mkdir('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/monthly_annom_plots_full/'+str(year))
        except FileExistsError:
            print('ok')
        for month in months:
            print('starting month: '+month)
            data_dict = pull_data(month, year)
            if data_dict[0] == 'n':
                print('skipping month: '+month+'/'+str(year))
                continue
            data_sel = go_stats.build_stats(np.squeeze(data_dict, axis=1), 1)
            z = data_dict[3]
            lat = data_dict[1]
            data_fetch_dict = build_ann_avg.get_monthly_data(int(month))
            data_fetch = data_fetch_dict[0]

#compare data
            anom = data_sel-data_fetch
#plot_data
#        mast_plot.make_zonal_anom_plot(anom, lat, z, str(month)+'/'+str(year)+' anomaly plot', month, str(year))
            mast_plot.make_mast_anom_plot(data_fetch, anom, data_sel, lat, z, str(month), str(year))
