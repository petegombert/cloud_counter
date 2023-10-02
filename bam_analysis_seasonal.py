"""
*This code has a similar purpose to the bam_analysis.py code, but there will be no seasonality
removal and thus will only look at individulal seasons.

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
from pyhdf.HDF import *
from pyhdf.VS import *
import datetime as dt
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
import merra_analize
import indexes
import go_stats
import mast_plot
import overview
import oned_corr
import eof_analysis
import matplotlib
import import_eke_data
import make_raw_bam_anim
import bam_analysis

months_arr = np.arange(6,9)
months = go_stats.make_np_list(months_arr)
months = go_stats.fix_days(months)
#months = ['12','01','02']
season = 'JJA'
years_arr = np.arange(2006,2017)
years_arr_delete = np.where(years_arr == 2011)
years_arr = np.delete(years_arr, years_arr_delete)
years = go_stats.make_np_list(years_arr)

def rest_season_bam(bam_mast, months, mix_years, bamlag):
    """
    This function will restrict the season of the bam index to whatever season is selected.
    Season is selected by the months inputed
    """
    #*Should just make a check using the months and this method:
    #This is the probably best way of doing it (getting the exact months)
    # for month in months:
    #     hit_arr = np.argwhere(bam_mast[1,:] == month)
    #     hstep = np.arange(0,hit_arr.shape[0])
        #for h in hstep:
            #bam_index_h = 

    #This is the easier way of doing it by doing the math on it
    hold = 12
    tstep = np.arange(0,years_arr.shape[0])
    # if mix_years == 1:
    #     tstep = tstep[:-1]
    for t in tstep:
        for month in months:
            indx = int(month)%12
            if indx == 0:
                indx = 12-bamlag
            indx = (indx-1-bamlag)+hold
            if mix_years == 1 and month == '01' or month == '02':
                indx=indx+12
            bam_arr_h = bam_mast[0,indx]
            #This code will save seasonal monthly values.
            # if t == 0 and month == months[0]:
            #     bam_arr = np.array([bam_arr_h])
            #     continue
            # bam_arr = np.append(bam_arr, bam_arr_h)
        #hold = hold+12

            #This code will save seasonal means.
            if month == months[0]:
                season_hold = np.array([bam_arr_h])
                continue
            season_hold = np.append(season_hold, bam_arr_h)
        hold = hold+12
        season_val = np.mean(season_val)
        if t == 0:
            bam_arr = np.array([season_val])
            continue
        bam_arr = np.append(bam_arr, season_val)

    print(bam_arr.shape)
    exit()
        

    return bam_arr

def rest_season_bam_2(bam_mast, years, months, mix_years, bamlag, season_mean):
    """
    This function will restrict the season of the bam index to whatever season is selected.
    Season is selected by the months inputed
    Must pass arrays in months and years.
    """
    #*Should just make a check using the months and this method:
    #This is the probably best way of doing it (getting the exact months)
    # for month in months:
    #     hit_arr = np.argwhere(bam_mast[1,:] == month)
    #     hstep = np.arange(0,hit_arr.shape[0])
        #for h in hstep:
            #bam_index_h = 

    #This is the easier way of doing it by doing the math on it
    s = 0
    tstep = np.arange(0, bam_mast.shape[1])
    for t in tstep:
        print(bam_mast[1,t])
        print(type(bam_mast[1,t]))
        year = bam_mast[1,t].year
        month = bam_mast[1,t].month
        if s < years.shape[0]:
            if year == years[s] and month == months[0]:
                print(year)
                print(month)
                print('getting bam vals!')
                season_hold = np.mean([bam_mast[0,t:t+3]])
                if s == 0:
                    season_arr = np.array([season_hold])
                    s = s+1
                    continue
                season_arr = np.append(season_arr, season_hold)
                s = s+1

    return season_arr

if __name__ == '__main__':
    cld_occ, lat, lon = oned_corr.eof_call_main(years, months, 0, 0)
    # cld_occ_zonal_mean = np.mean(cld_occ, axis=1)
    # cld_occ_area_mean = np.mean(cld_occ_zonal_mean, axis=1)
    # zeros = np.zeros(cld_occ_area_mean.shape[0])
    # mast_plot.seasonal_cycle_test(zeros, cld_occ_area_mean)
    # exit()
    # dates_occ = indexes.make_dt_arr(np.arange(1,2), 11, 3, '06-01-2006', 2011,)
    # print(dates_occ.shape)
    # exit()
    # cld_occ, bd_arr = go_stats.test_cld_occ(cld_occ_test, dates_occ)
    print(cld_occ.shape)
    bd_arr = np.array([dt.date(1999,11,25)])
    eof1, bam_index_mast = bam_analysis.run_main_bam(bd_arr)
    bam_arr = rest_season_bam_2(bam_index_mast, years_arr, months_arr, 0, 1, 0)
    plt.plot(np.arange(0,bam_arr.shape[0]), bam_arr)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')
    exit()

    # std_arr = bam_analysis.build_cov_arr_1st_dem(cld_occ, 'reg_stderr', bam_arr, lat, lon, 20,20)
    # mast_plot.std_err_bam(std_arr)
    # exit()
    corr_arr = bam_analysis.build_cov_arr_1st_dem(cld_occ, 2, bam_arr, lat, lon, 10,-115, 0.35)
    reg_inc = bam_analysis.build_cov_arr_1st_dem(cld_occ, 'reg_itc', bam_arr, lat, lon, 20,20, 0.35)
    reg_slp = bam_analysis.build_cov_arr_1st_dem(cld_occ, 'reg_slp', bam_arr, lat, lon, 20,20, 0.35)

    mast_plot.bam_cld_corr(corr_arr, reg_inc, reg_slp)



    