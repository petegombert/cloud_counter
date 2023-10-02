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
import build_ann_anom

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

daytime_filter = 0
min_dbz = -39.


win_months = ['12', '01', '02']
som_months = ['06', '07', '08']
seas_sel = win_months

def get_ncdf_data():
    os.chdir('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/ncdfs')
    f = cdf.Dataset('jja_avg_lt2017.cdf')
    data_hold = f.variables['total_mean']
    data = data_hold[:]

    return data

if __name__ == '__main__':
    flag = 0 #flag is used to make dataset if there is no data in 06 or 12 month.
    years = np.arange(2007, 2019)
    for y in years:
        for t in seas_sel:
            print('starting month '+str(t)+' / '+str(y))
            if t =='06' or t=='12' or flag==1:
                print('making new dataset')
                data_hold = build_ann_anom.pull_data(t, y)
                if data_hold[0] == 'n':
                    flag=1
                    continue
                get_monthly_data_dict = np.array([data_hold[-1,0], data_hold[3,0], data_hold[1,0], data_hold[4,0]], dtype=object)
                flag = 0
                continue
            data_hold = build_ann_anom.pull_data(t, y)
            if data_hold[0] == 'n':
                continue
            get_monthly_data_dict[0] = get_monthly_data_dict[0]+data_hold[-1,0]
        print('making cloud occurrence stats')
        cld_occ = go_stats.build_stats_presum(get_monthly_data_dict[0], get_monthly_data_dict[3], get_monthly_data_dict[2], 1)
        avg = get_ncdf_data()
        anom = cld_occ-avg
        mast_plot.make_mast_anom_plot(avg, anom, cld_occ, get_monthly_data_dict[2], get_monthly_data_dict[1], 'jja', y)

