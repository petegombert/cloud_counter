"""
This code will build the average seasonal data with the new dataset. Should be very similar to the code to build the zonal average data.
*Issue is that the cld_tot data isn't perfect because it's daytime only now so we need to run the replacement function.

Leave notes: I should add in an averaging code that will build averages to let me get an overview.
!!!Not going to run right rn. Need to change the region_arr function to return the mean. I changed it for 1d_corr.py!!!

To do:
Should rework region_arr3d. I had to play with it a bit to get it to work for the levels and im not sure
if it'll work that good moving on. Code should error if it doesn't work though so it won't be an undiagnosed problem.
The problem in the region arrays have to do with the order of the bound/lat/lon arrays. Can either fix by adjusting the
latlon_bounds array that initalizes the whole sequence. Or test the way I have it now where it shouldn't
matter which way the data goes in. Either or. 

Still not sure what's going on with the region_arr function. Why do I need to add 1 sometimes but not 
others!!! I think it might be something to do with half steps.
About the above issues with the region_arr functions. Part of it definitely has to do with the order of
the lat and lon arrays & bounds. I think I've fixed that part but there is still the issue that I've jerry
riged by adding the one on the end, but ultimately has to do with half steps I think.
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
import get_data
import mast_plot
import go_stats
import build_ann_anom
import coarsify
import diveinnewdata


months = ['06','07','08']
#months = ['06']
years = [2006, 2007, 2008, 2009, 2010, 2012, 2013, 2014, 2015, 2016]
z_bot = 0
z_top = 1000
dz_bot = 0
dz_top = 1500

def get_files(year):
    for month in months:
        data_id = str(year)+month
        if month == '06':
            data_arr = diveinnewdata.find_files(data_id, 2)
            continue
        data_hold = diveinnewdata.find_files(data_id, 2)
        data_arr = np.append(data_arr, data_hold, axis=0)

    return data_arr

def get_arr_index(dz_index, dz_t, dz_b, z_index, z_t, z_b, dz_bin, z_bin):
    """
    This function gets the indexes in the array that mark the bounds for the dz & cbtz bins
    dz_dbin tells us the width of the bins, while dz_bin tells the center value.
    """
#needs some work. The np.argwhere part doesn't like the and...erg
    dz_index_bot = dz_index - (dz_bin)
    dz_index_top = dz_index + (dz_bin)
    z_index_bot = z_index - (z_bin)
    z_index_top = z_index + (z_bin)

    print('Getting dz: ' +str(dz_b)+'-'+str(dz_t) +' Getting z: ' +str(z_b)+'-'+str(z_t))
    dz_bot = np.argwhere(dz_index_bot >= dz_b)
    dz_top = np.argwhere(dz_index_top <= dz_t)
    dz_bot = int(dz_bot[0])
    dz_top = int(dz_top[-1])

    z_bot = np.argwhere(z_index_bot >= z_b)
    z_top = np.argwhere(z_index_top <= z_t)
    z_bot = int(z_bot[0])
    z_top = int(z_top[-1])

    data_arr = np.array([dz_bot, dz_top, z_bot, z_top])

    return data_arr

def test_data(data_send):
    print('testing')
    hits = np.argwhere(np.isnan(data_send))
    print(hits.shape)
    if hits.shape[0] > 0:
        tstep = np.arange(0, hits.shape[0])
        print(tstep.shape)
        for t in tstep:
            hit_lat = hits[t][0]
            hit_lon = hits[t][1]
            data_send[hit_lat, hit_lon] = 0
        return data_send
    else:
        print('nerp')
        return data_send

def restrict_domain(lat, lon, lat_top, lat_bottom, lon_right, lon_left):
    hld_arr = abs(lat-lat_top)
    lat_hit_top = np.where(hld_arr == np.min(hld_arr))
    hld_arr = abs(lat-lat_bottom)
    lat_hit_bot = np.where(hld_arr == np.min(hld_arr))
    hld_arr = abs(lon-lon_right)
    lon_hit_right = np.where(hld_arr == np.min(hld_arr))
    hld_arr = abs(lon-lon_left)
    lon_hit_left = np.argwhere(hld_arr == np.min(hld_arr))

    return np.array([lat_hit_top, lat_hit_bot, lon_hit_right, lon_hit_left], dtype=object)

def region_arr(data_send, lat_send, lon_send, bounds_arr):
    t_bnd = int(bounds_arr[0])
    b_bnd = int(bounds_arr[1])
    r_bnd = int(bounds_arr[2])
    l_bnd = int(bounds_arr[3])
    
    if t_bnd>b_bnd:
        if r_bnd>l_bnd:
            data = data_send[b_bnd:t_bnd,l_bnd:r_bnd]
            lats = lat_send[b_bnd:t_bnd]
            lons = lon_send[l_bnd:r_bnd]
        else:
            data = data_send[b_bnd:t_bnd,r_bnd:l_bnd]
            lats = lat_send[b_bnd:t_bnd]
            lons = lon_send[r_bnd:l_bnd]
    elif b_bnd>t_bnd:
        if r_bnd>l_bnd:
            data = data_send[t_bnd:b_bnd+1,l_bnd:r_bnd+1]
            lats = lat_send[t_bnd:b_bnd+1]
            lons = lon_send[l_bnd:r_bnd+1]
        else:
            data = data_send[t_bnd:b_bnd+1,r_bnd:l_bnd+1]
            lats = lat_send[t_bnd:b_bnd+1]
            lons = lon_send[r_bnd:l_bnd+1]
#    mast_plot.low_cloud(lats, lons, data)
#    exit()
    #data_flt = data.flatten()
    #mean = np.mean(data_flt)

    return data, lats, lons

def region_arr3d(data_send, lat_send, lon_send, bounds_arr):
    t_bnd = int(bounds_arr[0])
    b_bnd = int(bounds_arr[1])
    r_bnd = int(bounds_arr[2])
    l_bnd = int(bounds_arr[3])
    if b_bnd > t_bnd:
        data = data_send[:,t_bnd:b_bnd+1,l_bnd:r_bnd+1]
        lats = lat_send[t_bnd:b_bnd+1]
        lons = lon_send[l_bnd:r_bnd+1]
    elif t_bnd > b_bnd:
        data = data_send[:,b_bnd:t_bnd,l_bnd:r_bnd]
        lats = lat_send[b_bnd:t_bnd]
        lons = lon_send[l_bnd:r_bnd]

    return data, lats, lons

def get_index():
    indx = np.genfromtxt('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/index_values/meiv2.txt', delimiter='	')
    print(indx.shape)

def go_stats_func(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr):
    """
    This function takes the cloud count data and returns the cloud occurrence statistics. It also filters the data by 
    dz & cbtz. 
    Returns:
        array: cld_occ_send, lat_vector, lon_vector
    """
    tstep = np.arange(0, cld_cnt_raw.shape[0])
    data_inx = get_arr_index(restrict_arr[0], restrict_arr[1], restrict_arr[2], restrict_arr[3],
                                restrict_arr[4], restrict_arr[5], restrict_arr[6], restrict_arr[7])
    for t in tstep:
        cld_cnt_raw[t] = cld_cnt_raw[t][data_inx[0]:data_inx[1]+1,data_inx[2]:data_inx[-1]+1,1]
        cld_cnt_raw[t] = np.sum(cld_cnt_raw[t], axis=0)
        cld_cnt_raw[t] = np.sum(cld_cnt_raw[t], axis=0)
        cld_tot_raw[t] = cld_tot_raw[t][1]

    cld_cnt = np.sum(cld_cnt_raw, axis=0)
    cld_tot = np.sum(cld_tot_raw, axis=0)
    cld_cnt = np.rot90(cld_cnt)
    cld_tot = np.rot90(cld_tot)

    cld_cnt_crs_arr = coarsify.make_me_coarse_2d(cld_cnt, lat, lon)
    cld_tot_crs_arr = coarsify.make_me_coarse_2d(cld_tot, lat, lon)

    lat_vector = cld_cnt_crs_arr[1]
    lon_vector = cld_cnt_crs_arr[2]

    cld_cnt_crs = cld_cnt_crs_arr[0]
    cld_tot_crs = cld_tot_crs_arr[0]

    cld_occ_send = cld_cnt_crs/cld_tot_crs
    cld_occ_send = np.flip(cld_occ_send[0], axis=0)
    cld_occ_send = test_data(cld_occ_send)
    
    return np.array([cld_occ_send, lat_vector, lon_vector], dtype=object)

def go_stats_func_wm(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr):
    """
    This function does the transformations that returns cloud occurrence frequency. 
    This function will return the data each month, the one above is the annual mean values.
    """
    tstep = np.arange(0, cld_cnt_raw.shape[0])
    data_inx = get_arr_index(restrict_arr[0], restrict_arr[1], restrict_arr[2], restrict_arr[3],
                                restrict_arr[4], restrict_arr[5], restrict_arr[6], restrict_arr[7])
    for t in tstep:
        cld_cnt_raw[t] = cld_cnt_raw[t][data_inx[0]:data_inx[1]+1,data_inx[2]:data_inx[-1]+1,1]
        cld_cnt_raw[t] = np.sum(cld_cnt_raw[t], axis=0)
        cld_cnt_raw[t] = np.sum(cld_cnt_raw[t], axis=0)

        cld_tot_raw[t] = cld_tot_raw[t][1]

        cld_cnt_raw[t] = np.rot90(cld_cnt_raw[t])
        cld_tot_raw[t] = np.rot90(cld_tot_raw[t])

        cld_cnt_crs_arr = coarsify.make_me_coarse_2d(cld_cnt_raw[t], lat, lon)
        cld_tot_crs_arr = coarsify.make_me_coarse_2d(cld_tot_raw[t], lat, lon)

        lat_vector = cld_cnt_crs_arr[1]
        lon_vector = cld_cnt_crs_arr[2]

        cld_cnt_crs = cld_cnt_crs_arr[0]
        cld_tot_crs = cld_tot_crs_arr[0]

        cld_occ_hold = cld_cnt_crs/cld_tot_crs
        cld_occ_hold = np.flip(cld_occ_hold[0], axis=0)
        cld_occ_hold = test_data(cld_occ_hold)

        if t == 0:
            cld_occ_send = np.array([cld_occ_hold])
            continue
        cld_occ_send = np.append(cld_occ_send, np.array([cld_occ_hold]), axis=0)
    
    return np.array([cld_occ_send, lat_vector, lon_vector], dtype=object)

if __name__ == '__main__':
    for year in years:
        print(year)
        raw_data = get_files(year)
        lat = raw_data[0,0]
        lon = raw_data[0,1]
        cbtz_bin = raw_data[0,2]
        cbtz_dbin = raw_data[0,3]
        dz_bin = raw_data[0,4]

        dz_dbin = raw_data[0,5]

        restrict_arr = np.array([dz_bin, dz_top, dz_bot, cbtz_bin, z_top, z_bot, dz_dbin, cbtz_dbin])

        cld_cnt_raw = raw_data[:,6]
        cld_tot_raw = raw_data[:,-2]

        go_stats_arr = go_stats_func(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr)
        lat_vector = go_stats_arr[1]
        lon_vector = go_stats_arr[2]
        cld_occ = go_stats_arr[0]

        bounds_arr = restrict_domain(lat_vector[0], lon_vector[0], 37, 10, -165, -116)
        region_arr(cld_occ, lat_vector[0], lon_vector[0], bounds_arr)
        if year == 2006:
            mean_arr = np.array([region_arr(cld_occ, lat_vector[0], lon_vector[0], bounds_arr)])
            year_arr = np.array([year])
            continue
        mean_hold = region_arr(cld_occ, lat_vector[0], lon_vector[0], bounds_arr)
        mean_arr = np.append(mean_arr, mean_hold)
        year_arr = np.append(year_arr, year)
    index = np.array([0.4, -0.7, -0.95, 0.55, -2.5, 0.2, -0.65, 0.3, 1.8, -0.35])
    mast_plot.corr(year_arr, index, mean_arr)
    exit()
    mast_plot.time_series(year_arr, index, mean_arr)
    exit()

    plt.plot(year_arr, mean_arr)
    plt.plot(year_arr, index)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/test.png')

    mast_plot.low_cloud_wregion(lat_vector, lon_vector, cld_occ, bounds_arr)
