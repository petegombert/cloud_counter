"""
This code recreates low cloud figure on the data website. It's just another test with the new data. Should be a nice figure to put into my prospectus for ENSO.

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
import diveinnewdata
import coarsify

def get_monthly_data(month):
    if month < 10:
        month = '0'+str(month)
    else:
        month = str(month)
    os.chdir('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj')
    data = diveintonewdata.find_files(month, 0)

    return data

def fix_tot_data(data):
    hits = np.argwhere(data == 0)
    hitstep = np.arange(0, hits.shape[0])
    for h in hitstep:
        lat = hits[h][0]
        lon = hits[h][1]
        data[lat,lon] = 1

    return data

if __name__ == '__main__':
    data = diveinnewdata.find_files(0, 'all')
    test = data[5,-2]
    ref = data[5, -1]
    lat = data[0,0]
    lon = data[0,1]
    cld_occ = test[:,0:2,:]
    cld_occ = np.sum(cld_occ, axis=0)
    cld_occ = np.sum(cld_occ, axis=0)
    cld_occ = np.sum(cld_occ, axis=0)
    ref = np.sum(ref, axis=0)
    cld_occ = np.rot90(cld_occ)
    ref = np.rot90(ref)
    ref = fix_tot_data(ref)
    crs_arr = coarsify.make_me_coarse_2d(cld_occ, lat, lon)
    ref_crs = coarsify.make_me_coarse_2d(ref, lat, lon)
    lat_vec = crs_arr[1]
    lon_vec = crs_arr[2]
    crs_occ = crs_arr[0]
    crs_occ_frq = crs_occ/ref_crs[0]
    mast_plot.low_cloud(lat_vec, lon_vec, crs_occ_frq)

