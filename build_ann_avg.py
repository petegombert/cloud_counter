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

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

year_start = 2006
daytime_filter = 0
min_dbz = -39.

def get_monthly_data(month):
#must change month to have leading 0 if it is less that 10.
    if month < 10:
        month = '0'+str(month)
    else:
        month = str(month)
    os.chdir('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj')
    print('fetching data')
    data = get_data.find_files(month, 0)
    get_data.test_monthly_data(data)
    z = data[0,3]
    lat = data[0,1]

#calls go_stats.py to build cloud_occ stats
    cld_occ = go_stats.build_stats(data, 0)

    return cld_occ, lat, z

if __name__ == '__main__':
    tstep = np.arange(1, 13)
    for t in tstep:
        print('starting month '+str(t))
        monthly_avg_dict = get_monthly_data(t)
        mast_plot.make_zonal_mean_plot(monthly_avg_dict[0], monthly_avg_dict[1], monthly_avg_dict[2], 'Zonal averaged cloud occurrence. Month: '+str(t)+'/ 2006-2019', str(t))

