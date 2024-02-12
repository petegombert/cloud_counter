"""
* This code is used for exploring the different raw cloud occurrence statistics.

To do:
    Not sure why but when I put a list of the all the months, the get_data function wasn't working 
    exactly right.

Leave notes:

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
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

#months = ['06','07','08']
#months = ['12','01','02']
#months = 'all'
#years = ['2010']
months = ['03','04','05']
# years_arr = np.arange(2006,2010)
# years = go_stats.make_np_list(years_arr)
years = ['2007','2008','2009','2010','2013','2014','2015','2016']
#High ENSO Years
#years = ['2006','2009','2012','2015']
#Low ENSO Years
#years = ['2007','2008','2009','2010','2012','2013','2014','2015','2016']
#years = ['2007']
z_bot = 0
z_top = 14000
dz_bot = 10000
dz_top = 14000

restrict_domain_call = 'no'
#restrict_domain_call = 'no'
#domain_restriction = np.array([-5, -30, -100, -70])
#domain_restriction = np.array([37, 10, -140, -116])
#domain_restriction = [48,25,-124,-113]
domain_restriction = [-40,-80,-78,-180]

mix_years =0

def get_data(years, months, mix_years):
    if len(years) == 1:
        for month in months:
            if years == 'all':
                data_hold = diveinnewdata.find_files(month, 0)
            else:
                if mix_years == 1 and month == '01':
                    years[0] = int(years[0])+1
                data_id = str(years[0])+month
                data_hold = diveinnewdata.find_files(data_id, 2)
                if data_hold == 0:
                        print('omitting: '+ str(years[0])+ '-'+ month)
                        continue
            try:
                data_arr = np.append(data_arr, data_hold, axis=0)
            except UnboundLocalError:
                data_arr = data_hold
        
        return data_arr

    else:
        for year in years:
            if months == 'all':
                data_hold = diveinnewdata.find_files(year, 1)
                try:
                    data_arr = np.append(data_arr, data_hold, axis=0)
                except UnboundLocalError:
                    data_arr = data_hold
            else:
                for month in months:
                    if mix_years == 1 and month == '01':
                        year = int(year)+1
                    data_id = str(year)+month
                    print(data_id)
                    data_hold = diveinnewdata.find_files(data_id, 2)
                    if data_hold == 0:
                        print('omitting: '+ str(year)+ '-'+ month)
                        continue
                    try:
                        data_arr = np.append(data_arr, data_hold, axis=0)
                    except UnboundLocalError:
                        data_arr = data_hold
           
        return data_arr 
    
def call_main(years, months, mix_years, restrict_domain_call, domain_restriction):
    print(years)
    raw_data = get_data(years, months, mix_years)
    lat = raw_data[0,0]
    lon = raw_data[0,1]
    cbtz_bin = raw_data[0,2]
    cbtz_dbin = raw_data[0,3]
    dz_bin = raw_data[0,4]
    dz_dbin = raw_data[0,5]
    cld_cnt_raw = raw_data[:,6]
    cld_tot_raw = raw_data[:,-2]

    indx_arr = jja_global_avg.get_arr_index(dz_bin, dz_top, dz_bot, cbtz_bin, z_top, z_bot, dz_dbin, cbtz_dbin)
    tstep = np.arange(0, cld_cnt_raw.shape[0])
    for t in tstep:
        cld_cnt_raw[t] = cld_cnt_raw[t][indx_arr[0]:indx_arr[1]+1,indx_arr[2]:indx_arr[-1]+1,1]
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
    cld_occ = cld_cnt_crs/cld_tot_crs
    cld_occ = np.flip(cld_occ[0], axis=0)

    if restrict_domain_call == 'yes':
        bounds_arr = jja_global_avg.restrict_domain(lat_vector[0], lon_vector[0], domain_restriction[0], domain_restriction[1],
                                                domain_restriction[2], domain_restriction[3])
        cld_occ, lat_vector, lon_vector = jja_global_avg.region_arr(cld_occ, lat_vector[0], lon_vector[0], bounds_arr)
    
    
    cld_occ = jja_global_avg.test_data(cld_occ)

    return lat_vector, lon_vector, cld_occ 

if __name__ == '__main__':
    raw_data = get_data(years, months, mix_years)
    lat = raw_data[0,0]
    lon = raw_data[0,1]
    cbtz_bin = raw_data[0,2]
    cbtz_dbin = raw_data[0,3]
    dz_bin = raw_data[0,4]
    dz_dbin = raw_data[0,5]
    cld_cnt_raw = raw_data[:,6]
    cld_tot_raw = raw_data[:,-2]

    indx_arr = jja_global_avg.get_arr_index(dz_bin, dz_top, dz_bot, cbtz_bin, z_top, z_bot, dz_dbin, cbtz_dbin)
    tstep = np.arange(0, cld_cnt_raw.shape[0])
    for t in tstep:
        cld_cnt_raw[t] = cld_cnt_raw[t][indx_arr[0]:indx_arr[1]+1,indx_arr[2]:indx_arr[-1]+1,1]
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

    if restrict_domain_call == 'yes':
        bounds_arr = jja_global_avg.restrict_domain(lat_vector[0], lon_vector[0], domain_restriction[0], domain_restriction[1],
                                                domain_restriction[2], domain_restriction[3])

    cld_cnt_crs = cld_cnt_crs_arr[0]
    cld_tot_crs = cld_tot_crs_arr[0]
    cld_occ = cld_cnt_crs/cld_tot_crs
    cld_occ = np.flip(cld_occ[0], axis=0)
    cld_occ = jja_global_avg.test_data(cld_occ)

    if restrict_domain_call == 'yes':
        print(bounds_arr)
        mast_plot.overview_plt_rd(lat_vector, lon_vector, cld_occ, dz_bot, dz_top, 
                                z_bot, z_top, bounds_arr, 'reg_anal_test')
    else:
        mast_plot.overview_plt(lat_vector, lon_vector, cld_occ, dz_bot, dz_top, z_bot, z_top, 'test_mam_bam_conv')



