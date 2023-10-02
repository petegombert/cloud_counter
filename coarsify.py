#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

* This code takes count data arrays and builds occurance statistics using the methodology from jays code.
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

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

def sum_data(data, start_index, end_index):
    tstep = np.arange(start_index, end_index)
    for t in tstep:
        if t ==start_index:
            data_coarse = data[:,t]
            continue

        data_hold = data[:,t]
        data_coarse = data_coarse+data_hold

    return data_coarse


def make_me_coarse(data, lat):
#make lat vector
    avg_factor=6
    lat_vector = np.arange(lat[0], lat[-1], avg_factor)
#    print(lat_vector)
#    exit()
    tstep = np.arange(0, lat_vector.shape[0]-1)
    for t in tstep:
        if t == 0:
            hits = np.argwhere(lat <= lat_vector[t+1])
            start_index = hits[0]
            end_index = hits[-1]
            index_hold = hits[-1]
            data_coarse = np.array([sum_data(data, start_index, end_index)])
            continue
        hits = np.argwhere(lat <= lat_vector[t+1])
        start_index = index_hold+1
        end_index = hits[-1]
        index_hold = end_index
        data_hold = sum_data(data, start_index, end_index)
        data_coarse = np.append(data_coarse, np.array([data_hold]), axis =0)

    data_coarse = np.flip(data_coarse, axis=0)
    data_coarse = np.rot90(data_coarse, k=3)
    return data_coarse
#    return np.rot90(data_coarse, k=3)

def go_sum_2d(data, start_index_lon, end_index_lon, start_index_lat, end_index_lat):
    col_sum = np.sum(data[int(start_index_lat):int(end_index_lat),int(start_index_lon):int(end_index_lon)])

    return col_sum

def make_me_coarse_2d(data, lat, lon):
    print(lon.shape)
    avg_factor = 8
    lat_vector = np.arange(lat[0], lat[-1], avg_factor)
    lon_vector = np.arange(lon[0], lon[-1], avg_factor)
    lat_step = np.arange(0, lat_vector.shape[0]-1)
    lon_step = np.arange(0, lon_vector.shape[0]-1)
    for lo in lon_step:
        if lo == 0:
            hits_lon = np.argwhere(lon <= lon_vector[lo+1])
            start_index_lon = hits_lon[0]
            end_index_lon = hits_lon[-1]
            index_hold_lon = hits_lon[-1]
        hits_lon = np.argwhere(lon <= lon_vector[lo+1])
        start_index_lon = index_hold_lon+1
        end_index_lon = hits_lon[-1]
        index_hold_lon = end_index_lon
        for la in lat_step:
            if la == 0:
                hits_lat = np.argwhere(lat <= lat_vector[la+1])
                start_index_lat = hits_lat[0]
                end_index_lat = hits_lat[-1]
                index_hold_lat = hits_lat[-1]
                sum_hold = np.array([go_sum_2d(data, start_index_lon, end_index_lon, start_index_lat, end_index_lat)])
                continue
            hits_lat = np.argwhere(lat <= lat_vector[la+1])
            start_index_lat = index_hold_lat+1
            end_index_lat = hits_lat[-1]
            index_hold_lat = end_index_lat
            sum_hold = np.append(sum_hold, np.array([go_sum_2d(data, start_index_lon, end_index_lon, start_index_lat, end_index_lat)]))
        sum_hold_send = np.array([sum_hold])
        sum_hold_send = np.rot90(sum_hold_send)
        if lo == 0:
            data_sum = sum_hold_send
            continue
        data_sum = np.append(data_sum, sum_hold_send, axis=1)
    data_sum = np.flip(data_sum, axis=0)
#    data_sum = np.rot90(data_sum, k=3)

    return np.array([[data_sum], [lat_vector], [lon_vector]], dtype=object)

def make_me_coarse_2d_na_return(data, lat, lon):
    print(lon.shape)
    print(lat.shape)
    avg_factor = 1.5
    lat_vector = np.arange(lat[0], lat[-1], avg_factor)
    lon_vector = np.arange(lon[0], lon[-1], avg_factor)
    print(lat_vector.shape, lon_vector.shape)
    lat_step = np.arange(0, lat_vector.shape[0]-1)
    lon_step = np.arange(0, lon_vector.shape[0]-1)
    for lo in lon_step:
        if lo == 0:
            hits_lon = np.argwhere(lon <= lon_vector[lo+1])
            start_index_lon = hits_lon[0]
            end_index_lon = hits_lon[-1]
            index_hold_lon = hits_lon[-1]
        hits_lon = np.argwhere(lon <= lon_vector[lo+1])
        start_index_lon = index_hold_lon+1
        end_index_lon = hits_lon[-1]
        index_hold_lon = end_index_lon
        for la in lat_step:
            if la == 0:
                hits_lat = np.argwhere(lat <= lat_vector[la+1])
                start_index_lat = hits_lat[0]
                end_index_lat = hits_lat[-1]
                index_hold_lat = hits_lat[-1]
                data_arr = np.mean(data[int(start_index_lat):int(end_index_lat),int(start_index_lon):int(end_index_lon)])
                sum_hold = np.array([data_arr])
                #sum_hold = np.array([go_sum_2d(data, start_index_lon, end_index_lon, start_index_lat, end_index_lat)])
                continue
            hits_lat = np.argwhere(lat <= lat_vector[la+1])
            start_index_lat = index_hold_lat+1
            end_index_lat = hits_lat[-1]
            index_hold_lat = end_index_lat
            sum_hold = np.append(sum_hold, np.array([data_arr]))
            #sum_hold = np.append(sum_hold, np.array([go_sum_2d(data, start_index_lon, end_index_lon, start_index_lat, end_index_lat)]))
        sum_hold_send = np.array([sum_hold])
        #sum_hold_send = np.rot90(sum_hold_send)
        if lo == 0:
            data_sum = sum_hold_send
            continue
        data_sum = np.append(data_sum, sum_hold_send, axis=0)
    #data_sum = np.flip(data_sum, axis=0)
    #data_sum = np.rot90(data_sum, k=3)

    return data_sum, lat_vector, lon_vector

if __name__ == '__main__':
    data_sel = get_data.find_files('08', 0)
    data_send = data_sel[0,-1]
    data_send = data_send[0,3,:,10,:]
    test = np.argwhere(data_send != 0)
    data = make_me_coarse(data_send, data_sel[0, 1])

