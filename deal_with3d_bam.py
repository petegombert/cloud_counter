"""
*This code will plot and handle the data from the 3d (time,lon,height,lat) data and makes figures 
with it 

To do:

Leave notes:

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
import matplotlib
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
import bam_merra_calc_smal
import make_bam_index_smal
import go_stats
import mast_plot
import new_make_raw_bam_anim
import make_raw_bam_anim
import eof_analysis
import bam_analysis

years_arr = np.arange(1986,1987)
years = go_stats.make_np_list(years_arr)
months = ['06','07','08']
#months = ['01','02','03','04','05','06','07','08','09','10','11','12']
days_arr = np.arange(1, 32)
days_send = go_stats.make_np_list(days_arr)
days = go_stats.fix_days(days_send)

latlon_bounds = [-20,-70,180,-180]
lat_height_arr = [-70,-20,900,200]
mix_years = 0

bam_smal_data = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bamtest'
bam_data = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bam/'

def fiveh_mb_height_get(data_s, level_s):
    min_arr = abs(level_s-500)
    hit_arr = np.argwhere(min_arr==np.min(min_arr))
    hit_arr = int(hit_arr[0])
    
    data_s_500 = data_s[:,:,hit_arr,:]

    return data_s_500


if __name__ == '__main__':
    level, lat, lon, u_anom, u3d, v_anom, v3d, eke_data, eke3d = new_make_raw_bam_anim.import_data3d()

    data_s_500 = fiveh_mb_height_get(eke3d, level)

    print(eke_data.shape)
    mean_eke = np.mean(eke_data, axis=0)
    #print(mean_eke.shape)
    #plt.imshow(mean_eke[:19,:])
    
    """
    !!! If using data from the following netcdf4s the u & v data must be weighted, otherwise
    it's already weighted.
        u_anom_test3d_1986_2.nc4 
        u_anom_test3d_1986_3.nc4
        u_anom_test3d_2008.nc4
    """
    #lat_weight_arr = go_stats.latitudinal_weight2(lat, u3d[4,10,:,:], 'cos')
    #u3d_w, v3d_w = make_raw_bam_anim.apply_weights(u3d, v3d, lat_weight_arr)
    mass_weight_arr = go_stats.mass_weight_3(level, eke3d[10,5,:,:])
    eke3d_m = np.mean(eke3d, axis=1)
    u3d_m = np.mean(u3d, axis=1)
    eke3d_w, u3d_w = make_raw_bam_anim.apply_weights(eke3d_m, u3d_m, mass_weight_arr)

    time_step, lev_step, lat_step = eke3d_w.shape
    eke3d_send = np.reshape(eke3d_w, (time_step,(lev_step*lat_step)))

    eof1, pc1 = eof_analysis.myeof(eke3d_send)
    corr_arr = bam_analysis.build_cov_arr_1st_dem(eke3d_w, 'reg', pc1)

    p = plt.imshow(corr_arr, cmap='jet')
    plt.colorbar(p)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testinger.png')
    exit()

    tstep = np.arange(u3d.shape[0])
    xstep = np.arange(u3d.shape[1])
    for t in tstep:
        print(t)
        eke_day = np.mean(eke3d[t,:,:,:], axis=0)
        eke_day = mass_weight_arr*eke_day
        print(eke_day.shape)
        if t == 0:
            #eke_mean_day = np.mean(u3d_w[t,:,:,:], axis=0)
            eke_mean_day = eke_day
        elif t !=0:
            #eke_mean_day = eke_data[t,:,:]
            eke_mean_day_hold = np.mean(eke3d[:t,:,:,:], axis=1)
            print(eke_mean_day_hold.shape)
            eke_mean_day = np.mean(eke_mean_day_hold, axis=0)
            eke_mean_day = mass_weight_arr*eke_mean_day

        mast_plot.eke_plot_loop(eke_day, eke_mean_day, data_s_500[t,:,:], level, lat, lon, t)