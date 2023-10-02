"""
* This script does the zonal mean cloud occurrence correlations with the BAM index.

To do:

Leave notes:

"""
"""
*This function is designed to compare the correlation between ENSO and BAM. The reasoning behind this
is because there 
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
#from concurrent.futures.process import _MAX_WINDOWS_WORKERS
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import scipy.stats as sp
from scipy.optimize import curve_fit
from pyEOF import *
import jja_global_avg
import requests
#import merra_first
import indexes
import go_stats
import mast_plot
import zonal_enso
import get_data
import matplotlib
import import_eke_data
import make_raw_bam_anim
import bam_analysis

months_arr = np.arange(1,13)
months = go_stats.make_np_list(months_arr)
months = go_stats.fix_days(months)

years_arr = np.arange(2006,2010)
years = go_stats.make_np_list(years_arr)

def remove_nan_vals(data_s):
    """
    This function loops through the time dimension of a 3d array and removes the nan values at 
    each timestep. Input must have shape (time, lat, lon). RETURNS: Array of same shape
    """
    tstep = np.arange(0,data_s.shape[0])
    for t in tstep:
        data_nan_r = jja_global_avg.test_data(data_s[t,:,:])
        if t == 0:
            data_nnan_arr = np.array([data_nan_r])
            continue
        data_nnan_arr = np.append(data_nnan_arr, np.array([data_nan_r]), axis=0)

    return data_nnan_arr


if __name__ == '__main__':
    for year in years:
        data_arr = zonal_enso.data_call(year, months, 0)
        lat = data_arr[0,1]
        lon = data_arr[0,2]
        z = data_arr[0,3]
        dbz = data_arr[0,4]

        cld_cnt_raw = data_arr[:,-1]
        cld_occ = go_stats.build_stats_wm(data_arr, 1)

        if year == years[0]:
            cld_occ_arr = cld_occ
            continue
        cld_occ_arr = np.append(cld_occ_arr, cld_occ, axis=0)

    # p = plt.imshow(np.mean(cld_occ_arr, axis=0))
    # plt.colorbar(p)
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')

    cld_occ_anom = go_stats.remove_mean(cld_occ_arr, 1500,0,1000,0, 1)
    # p = plt.imshow(np.mean(cld_occ_anom, axis=0))
    # plt.colorbar(p)
    #plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')
    cld_occ_anom = remove_nan_vals(cld_occ_anom)
    # p = plt.imshow(np.mean(cld_occ_anom, axis=0))
    # plt.colorbar(p)
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/tester.png')

    eof1, bam_index = bam_analysis.run_main_bam()

    # mcld_occ_anom = np.mean(cld_occ_anom, axis=1)
    # mcld_occ_anom_s = np.mean(mcld_occ_anom, axis=1)
    # mcld_occ_arr = np.mean(cld_occ_arr, axis=1)
    # mcld_occ_arr_s = np.mean(mcld_occ_arr, axis=1)
    # print(np.mean(mcld_occ_anom_s))
    # exit()

    # p = plt.imshow(np.mean(cld_occ_anom, axis=0))
    # plt.colorbar(p)
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')

    # print(mcld_occ_anom_s)
    # exit()

    # mast_plot.seasonal_cycle_test(mcld_occ_anom_s, mcld_occ_arr_s)
    # exit()

    corr_arr = bam_analysis.build_cov_arr_1st_dem(cld_occ_anom, 2, bam_index, z, lat, 12, 60)
    reg_intc = bam_analysis.build_cov_arr_1st_dem(cld_occ_anom, 'reg_itc', bam_index, z, lat, 10, 20)
    reg_slp = bam_analysis.build_cov_arr_1st_dem(cld_occ_anom, 'reg_slp', bam_index, z, lat, 10, 20)

    mast_plot.bam_cld_corr_zonal(corr_arr, reg_intc, reg_slp)



    # cmap_send = matplotlib.cm.get_cmap('RdBu')
    # shifted_cmap1 = mast_plot.shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')

    # p = plt.imshow(corr_arr, cmap = shifted_cmap1)
    # plt.colorbar(p)
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')

    


