"""
This code builds the zonal-mean correlation with ENSO. This is the atmospheric cross section code. 

Leave notes:
    Trying to get a zonal sam look to see if we can see the shift in the jet stream.

To do:
    Need to make data_call func resilient to missing data.
    Need to test the mix years part of the data_call.
    Need to test the comp function by adding anomalies and comparing it to the mean.
    Hard to make anything out of the data.
    When oned_corr.build_cov_arr is called, it is using the index from oned_corr.
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
import jja_global_avg
import overview
import oned_corr
import indexes
import build_seas_annom
import bam_analysis

# months = ['12','01','02']
# inx_months = ['11','12','01']
months_arr = np.arange(1,13)
months = go_stats.make_np_list(months_arr)
months = go_stats.fix_days(months)
months_inx_arr = np.arange(8,11)
months_inx = go_stats.make_np_list(months_inx_arr)
months_inx = go_stats.fix_days(months_inx)
#years = ['2006','2007','2008','2009','2010','2012','2013','2014','2015','2016']
years_arr = np.arange(2007,2017)
delete = np.argwhere(years_arr == 2011)
years_arr = np.delete(years_arr, delete)
# delete = np.argwhere(years_arr == 2012)
# years_arr = np.delete(years_arr, delete)

years = go_stats.make_np_list(years_arr)

#index = np.array([0.4, -0.7, -0.95, 0.55, -2.5, 0.2, -0.65, 0.3, 1.8, -0.35])
#index = indexes.enso(months)
index = indexes.sam(months_inx, years_arr, 0)


def data_call(year, months, mix_years):
    for month in months:
        if mix_years == 1 and month == '01':
            year = int(year)+1
        print(str(month))
        data_id = str(year)+str(month)
        print(data_id)
        data_hold = get_data.find_files(data_id, 2)
        if data_hold == 0:
            print('Ommitting: '+str(year)+' - '+month)
            continue
        try:
            data_arr = np.append(data_arr, data_hold, axis=0)
        except UnboundLocalError:
            data_arr = data_hold
    return data_arr

def get_yearly_occ_matrix(year):
    flag =0
    print('starting year: ' + str(year))
    for month in months:
        if month == '06' or flag == 1:
            data_hold = build_ann_anom.pull_data(month, year)
            if data_hold[0] == 'n':
                flag =1
                continue
            data_arr = np.array([data_hold[-1,0], data_hold[3,0], data_hold[1,0], data_hold[4,0]], dtype=object)
            flag=0
            continue
        data_hold = build_ann_anom.pull_data(month, year)
        if data_hold[0] == 'n': #not sure why I put this in there. The code is from build_seas_annom.py
            continue
        data_arr[0] = data_arr[0]+data_hold[-1,0]

    return data_arr

if __name__ == '__main__':

#** This code builds the correlation plots with ENSO
    for year in years:
        print('Starting year: '+ str(year))
        data_arr = data_call(year, months, 0)
        lat = data_arr[0,1]
        lon = data_arr[0,2]
        z = data_arr[0,3]
        dbz = data_arr[0,4]

        cld_cnt_raw = data_arr[:,-1]
        cld_occ = go_stats.build_stats(data_arr, 1)
        #cld_occ = np.flip(cld_occ, axis=0)

        if year == years[0]:
            cld_occ_arr = np.array([cld_occ])
            continue
        cld_occ_arr = np.append(cld_occ_arr, np.array([cld_occ]), axis=0)

    mean_state = np.mean(cld_occ_arr, axis=0)

    anom_cld_occ = cld_occ_arr-mean_state

        # if year == years[0]:
        #     cld_occ_hold = cld_occ
        #     continue
        # if year == years[1]:
        #     hold_test = oned_corr.go_append(cld_occ_hold, cld_occ, 1)
        #     continue
        # hold_test = oned_corr.go_append(hold_test, cld_occ, 0)

    # print(index)
    # cov_arr = oned_corr.build_cov_arr(hold_test, 2, index)
    # cov_arr = jja_global_avg.test_data(cov_arr)
    corr_arr, reg_slp, reg_inc = bam_analysis.build_cov_arr_1st_dem(cld_occ_arr, 'all', index, lat, lon, -40, 80, 0.06, 'noname')
    mast_plot.oned_zonal(reg_slp, corr_arr, reg_inc)

    # step = 0
    # indeces = indexes.enso(['07','08']) 
    # for index in indeces:
    #     if index >= 0.35:
    #         h_data_hold = get_yearly_occ_matrix(years[step])
    #         try:
    #             h_data[0] = h_data[0] + h_data_hold[0]
    #         except NameError:
    #             h_data = h_data_hold
            

    #     elif index <= -0.35:
    #         l_data_hold = get_yearly_occ_matrix(years[step])
    #         try:
    #             l_data[0] = l_data[0] + l_data_hold[0]
    #         except NameError:
    #             l_data = l_data_hold

    #     else:
    #         n_data_hold = get_yearly_occ_matrix(years[step])
    #         try:
    #             n_data[0] = n_data[0] = n_data_hold[0]
    #         except NameError:
    #             n_data = n_data_hold

            
    #     #n_data_hold = get_yearly_occ_matrix(years[step])
        
    #     step=step+1

    # h_cld_occ = go_stats.build_stats_presum(h_data[0], h_data[3], h_data[2], 1)
    # l_cld_occ = go_stats.build_stats_presum(l_data[0], l_data[3], l_data[2], 1)
    # n_cld_occ = go_stats.build_stats_presum(n_data[0], n_data[3], n_data[2], 1)
    # t_cld_occ = go_stats.build_stats_presum((h_data[0]+l_data[0]+n_data[0]), n_data[3], n_data[2], 1)

    # avg = build_seas_annom.get_ncdf_data()

    # mast_plot.zonal_enso_comps(avg[:12,:], h_cld_occ[:12,:], l_cld_occ[:12,:], h_data[1][:12], h_data[2])




        








