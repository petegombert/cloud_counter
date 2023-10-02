"""
This code analyses the merra data. First step is comparing the high and low ENSO year 500mb heights to one another.

To do:
    For the correlation analysis, I'm going to play around with keeping the monthly data. 
    I think that it might be useful for being able to look at time periods within the annual mean data.
    Clean up build_ds function. All the parts are still there, I think it might be pretty clunky.

Leave notes:
    In making the merra analysis for the correlation between enso & sam I obliterated the build_ds function. 
    

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
from warnings import filterwarnings
import jja_global_avg
from pydap.client import open_url
from pydap.cas.urs import setup_session
import requests
import merra_first
import indexes
import mast_plot

#merra_path = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/ENSO_SO'
#merra_path = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/FD_ENSO'
merra_path = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/SAM_GLB/'

def grab_data(filename,cat,year,month_sel):
    if month_sel == 1:
        os.chdir(os.path.join(merra_path,cat,year))
    else:
        os.chdir(os.path.join(merra_path,cat))
    print(filename)
    f = cdf.Dataset(filename)
    lat_get = f.variables['lat']
    lon_get = f.variables['lon']
    mb500_get = f.variables['mb500']
    slp_get = f.variables['slp']
    t2m_get = f.variables['t2m']
    pbltop_get = f.variables['pbltop']
    q500_get = f.variables['q500']

    lat = lat_get[:]
    lon = lon_get[:]
    mb500 = mb500_get[:]
    slp = slp_get[:]
    t2m = t2m_get[:]
    pbltop = pbltop_get[:]
    q500 = q500_get[:]
    f.close()

    return lat, lon, mb500, slp, t2m, pbltop, q500

def build_ds(year, cat, months):
    print('making mast array')

    #* This solves the data generation changes.
    if year == 'all' and months == 'all':
        for filename in os.listdir(os.path.join(merra_path,cat)):
            if filename[0] == 'k':
                print(filename)
                lat, lon, mb500, slp, t2m, pbltop, q500 = grab_data(filename,cat)
                try:
                    mast_arr = np.append(mast_arr, np.array([[lat, lon, mb500, slp, t2m, pbltop, q500]], dtype=object), axis=0)
                except UnboundLocalError:
                    mast_arr = np.array([[lat, lon, mb500, slp, t2m, pbltop, q500]], dtype=object)

        return mast_arr

    elif len(year) > 1:
        for years in year:
            for filename in os.listdir(os.path.join(merra_path,cat,years)):
                if filename.split('.')[0][4:] == months[0] or filename.split('.')[0][4:] == months[1] or filename.split('.')[0][4:] == months[2]:
                    print(filename)
                    lat, lon, mb500, slp, t2m, pbltop, q500 = grab_data(filename,cat,years,1)
                    try:
                        mast_arr = np.append(mast_arr, np.array([[lat, lon, mb500, slp, t2m, pbltop, q500]], dtype=object), axis=0)
                    except UnboundLocalError:
                        mast_arr = np.array([[lat, lon, mb500, slp, t2m, pbltop, q500]], dtype=object)
        return mast_arr

    #os.chdir(os.path.join(merra_path,cat))
    #os.mkdir(year)
    for filename in os.listdir(os.path.join(merra_path,cat)):
        #if filename[0:4] == year and filename[0] != 'k' and filename[0] != 'p' and filename != year:
         if filename[0] == 'k' and filename[0:4] == year:   
            print(filename)
            lat, lon, mb500, slp, t2m, pbltop, q500 = grab_data(filename,cat)
            #*Only turn on when making the merra mean dataset initially* 
            #os.remove(filename)
            #os.rename(filename, year+'/'+filename)
            try:
                mast_arr = np.append(mast_arr, np.array([[lat, lon, mb500, slp, t2m, pbltop, q500]], dtype=object), axis=0)
            except UnboundLocalError:
                mast_arr = np.array([[lat, lon, mb500, slp, t2m, pbltop, q500]], dtype=object)

    return mast_arr

def f_lat_lon(lat, lon, lat_s, lon_s):
    lat_hold = abs(lat-lat_s)
    lon_hold = abs(lon-lon_s)
    lat_hit = np.argwhere(lat_hold == np.min(lat_hold))
    lon_hit = np.argwhere(lon_hold == np.min(lon_hold))

    return int(lat_hit[0]), int(lon_hit[0])

def build_ts(years, cat, lat_s, lon_s):
    #**Years needs to be a np array
    os.chdir(os.path.join(merra_path,cat))
    tstep = np.arange(0, years.shape[0])
    for t in tstep:
        year = str(years[t])
        for filename in os.listdir(os.path.join(merra_path,cat)):
            if filename.split('.')[0][4:] == year:
                print('filename')
                lat, lon, mb500, slp, t2m, pbltop, q500 = grab_data(filename,cat)
                lat_indx, lon_indx = f_lat_lon(lat, lon, lat_s, lon_s)
                if t == 0:
                    data_ts = np.array([mb500[lat_indx, lon_indx]])
                    continue
                data_ts = np.append(data_ts, np.array([mb500[lat_indx, lon_indx]]))
    
    return data_ts

def get_mean_vals(data):
    mb500 = data[:,2]
    m_mb500 = np.mean(mb500)

    slp = data[:,3]
    m_slp = np.mean(slp)

    t2m = data[:,4]
    m_t2m = np.mean(t2m)

    pbltop = data[:,5]
    m_pbltop = np.mean(pbltop)

    q500 = data[:,6]
    m_q500 = np.mean(q500)

    return m_mb500, m_slp, m_t2m, m_pbltop, m_q500

if __name__ == '__main__':
    low_arr = build_ds('low')
    lat = low_arr[0,0]
    lon = low_arr[0,1]
    high_arr = build_ds('high')
    l_mb500, l_slp, l_t2m = get_mean_vals(low_arr)
    h_mb500, h_slp, h_t2m = get_mean_vals(high_arr)
    mast_plot.merra_1st_look(lat, lon, h_mb500, l_mb500)





