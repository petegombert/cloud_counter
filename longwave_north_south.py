"""
*This code is going to look into something I thought of a few weeks ago. The idea is this, 
how is the longwave pattern in the northern hemisphere related to the longwave pattern in the southern
hemisphere could one tell us about the other?

To to:

Leave notes: 

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
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
import jja_global_avg
import overview
import atmos_6040_final
import indexes
import oned_corr

years_n = go_stats.make_np_list(np.arange(2007,2018))
twen011 = years_n.pop(4)

years_s = go_stats.make_np_list(np.arange(2006,2017))
twen011 = years_s.pop(5)

months_low = ['06','07','08']
months_high = ['12','01','02']

def twod_corr(data1, data2, cor_sel):
    ystep = np.arange(data1.shape[0])
    xstep = np.arange(data1.shape[1])
    for y in ystep:
        for x in xstep:
            if cor_sel == 'corr':
                cov_m_hold = np.corrcoef(data1[y,x,:], data2[y,x,:])
                cov_hold = cov_m_hold[0,1]
            if x == 0:
                cov_x_arr = np.array([cov_hold])
                continue
            cov_x_arr = np.append(cov_x_arr, np.array([cov_hold]), axis=0)
        if y == 0:
            cov_arr = np.array([cov_x_arr])
            continue
        cov_arr = np.append(cov_arr, np.array([cov_x_arr]), axis =0)
    
    return cov_arr

if __name__ == '__main__':
    nh_win, lat, lon = oned_corr.eof_call_main(years_n, months_high, 1)
    nh_win = np.flipud(nh_win)

    sh_win, lat, lon = oned_corr.eof_call_main(years_s, months_low, 0)
    sh_win = np.flipud(sh_win)

    print(nh_win.shape)
    nh_win = np.mean(nh_win,axis=0)
    nh_win = np.mean(nh_win,axis=0)
    print(nh_win.shape)
    corr_arr = oned_corr.build_cov_arr(sh_win, 2, nh_win)
    #corr_arr = twod_corr(nh_win, sh_win, 'corr')
    print(corr_arr.shape)

    # std_test = np.std(corr_arr.flatten(), ddof=1)
    # std_2 = std_test*2
    # stat_arr = oned_corr.spatial_sig(corr_arr, std_2)

    # sig = np.argwhere(stat_arr == 1)
    # nsig = np.argwhere(stat_arr == 0)
    #lat_arr, lon_arr = oned_corr.get_lat_lon(corr_arr, sig, lat, lon)

    #mast_plot.oned_corr(lat, lon, corr_arr, lat_arr, lon_arr, 'JJA', 'DJF', '0-3000', '3000-6000', 'longwave_test')

    norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=corr_arr.min(), vmax=corr_arr.max())

    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1,projection=proj)
    ax.coastlines()
    ax.imshow(corr_arr, origin='lower', extent = [-180,180,-90,90,], aspect='equal', cmap = 'RdBu')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/tester.png')



