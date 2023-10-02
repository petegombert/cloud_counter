"""
*This code will combine the output arrays from the make_clim_bam.py scripts which produced multiple arrays
of seasonal means. Due to time limitations and submission limitations, they had to be run in chunks.
This code will import all the chunks and create seasonal mean files. 

To do:

Leave notes:
    Need to import the data for each season. 
    The son_2003-2013 data didn't run.

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
import bam_merra_calc
import go_stats
import make_bam_index
import make_clim_bam_smal

season = 'son'

bam_file_path = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bam/'

def merge_data(u_send, v_send, t_send, init):
    if init == 1:
        global u_arr
        u_arr = np.array([u_send])
        global v_arr
        v_arr = np.array([v_send])
        global t_arr
        t_arr = np.array([t_send])

    else:
        u_arr = np.append(u_arr, np.array([u_send]), axis=0)
        v_arr = np.append(v_arr, np.array([v_send]), axis=0)
        t_arr = np.append(t_arr, np.array([t_send]), axis=0)


if __name__ == '__main__':
    os.chdir(bam_file_path)
    t = 0
    for filename in os.listdir(bam_file_path):
        if season == filename.split('_')[0]:
            years = filename.split('_')[-2]+'_'+filename.split('_')[-1].split('.')[0]
            lat_h, lev_h, u_h, v_h, t_h = make_bam_index.get_seasonal_mean_yrs(season, years)
            if t == 0:
                print('making new data arrays')
                merge_data(u_h, v_h, t_h, 1)
            else:
                print('adding to existing data arrays')
                merge_data(u_h, v_h, t_h, 0)

            t=t+1
 
    mean_u = np.mean(u_arr, axis=0)
    mean_v = np.mean(v_arr, axis=0)
    mean_t = np.mean(t_arr, axis=0)

    fname = season+'_full_climatological_mean.nc4'
    make_clim_bam_smal.save_nc(mean_u, mean_v, mean_t, lat_h, lev_h, fname)




            


        


