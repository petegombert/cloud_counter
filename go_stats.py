#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

* This code takes count data arrays and builds occurance statistics using the methodology from jays code.

To do:

Leave notes:
    Still thinking about how to build mass weighting functions
    Not sure what T should be or how to calculate it. 
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
import julian
from datetime import datetime
import datetime as dt
import sys
import subprocess
from warnings import filterwarnings
import get_data
import mast_plot
import coarsify
import math
import indexes
import diveinnewdata
import jja_global_avg

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

def sum_months(data):
#sums all cloud_occ values. Returns array of cloud_occ with all months summed.
#checks for specific months
    if len(data.shape) < 2:
        return data[-1]

    tstep = np.arange(0, data.shape[0])
    for t in tstep:
        if t == 0:
            mast_arr = data[t][-1]
            continue
        hold = data[t][-1]
        mast_arr = mast_arr+hold

    return mast_arr

def sum_along_dbz(data, dbz, min_dbz):
#Month_sum will only be okay if we have the same values for z and dbz for all months. Should probably make it error if this is not the case.
    cloud_occ = data #assumes we always want daytime values. I think this is the case.
    if min_dbz != 1000:
        start_index = np.where(dbz==min_dbz)
        print(start_index)
    else:
        start_index = np.where(dbz==0) #finds where dbz equals 0, this will be the starting spot for the summation.
    tstep = np.arange(int(start_index[0]), dbz.shape[0])
    for t in tstep:
        if t == int(start_index[0]):
            sum = cloud_occ[t, :, :]
            continue
        hold = cloud_occ[t, :, :]
        sum = sum+hold

    return sum

def build_stats(data_send, daytime_flag):
#Data must come in as the whole dataset. Need to be able to pull the dbz bins for the summation.
#Daytime flag key: daytime_flag = 1: daytime only. daytime_flag = 0: nighttime only. daytime_flag = 2: all data.
    print('building stats')
    data = sum_months(data_send)
    if daytime_flag ==1:
        data = data[1,:,:,:,:]
    elif daytime_flag==0:
        data = data[0,:,:,:,:]
    else:
        data = np.sum(data, axis=0)
#sum along longitude:
    data = np.sum(data, axis=2)

#sum to get total array:
    data_tot = np.sum(data, axis=0)
#fix issue with multiple months/vs single month - annom sends single month
    if len(data_send.shape) < 2:
        dbz = data_send[4]
        lat = data_send[1]
    else:
        dbz = data_send[0,4]
        lat = data_send[0,1]
    data_rel = sum_along_dbz(data, dbz, -36.)
    data_rel = coarsify.make_me_coarse(data_rel, lat)
    data_tot = coarsify.make_me_coarse(data_tot, lat)

    data_occ = data_rel/data_tot

    return data_occ

def build_stats_wm(data_send, daytime_flag):
#Data must come in as the whole dataset. Need to be able to pull the dbz bins for the summation.
#Daytime flag key: daytime_flag = 1: daytime only. daytime_flag = 0: nighttime only. daytime_flag = 2: all data.
    print('building stats')
    tstep = np.arange(0, data_send.shape[0])
    for t in tstep:
        data = data_send[t,-1]

        if daytime_flag ==1:
            data = data[1,:,:,:,:]
        elif daytime_flag==0:
            data = data[0,:,:,:,:]
        else:
            data = np.sum(data, axis=0)
    #sum along longitude:
        data = np.sum(data, axis=2)

    #sum to get total array:
        data_tot = np.sum(data, axis=0)
    #fix issue with multiple months/vs single month - annom sends single month
        if len(data_send.shape) < 2:
            dbz = data_send[4]
            lat = data_send[1]
        else:
            dbz = data_send[0,4]
            lat = data_send[0,1]
        data_rel = sum_along_dbz(data, dbz, -36.)
        data_rel = coarsify.make_me_coarse(data_rel, lat)
        data_tot = coarsify.make_me_coarse(data_tot, lat)

        data_occ = data_rel/data_tot
        data_occ = np.flipud(data_occ)
        data_occ = jja_global_avg.test_data(data_occ)

        if t == 0:
            data_occ_arr = np.array([data_occ])
            continue
        data_occ_arr = np.append(data_occ_arr, np.array([data_occ]), axis=0)

    return data_occ_arr

def build_stats_presum(data_send, dbz, lat, daytime_flag):
    data = data_send
    if daytime_flag ==1:
        data = data[1,:,:,:,:]
    elif daytime_flag==0:
        data = data[0,:,:,:,:]
    else:
        data = np.sum(data, axis=0)

#sum along longitude:
    data = np.sum(data, axis=2)

#sum to get total array:
    data_tot = np.sum(data, axis=0)
    data_rel = sum_along_dbz(data, dbz, -36.)
    data_rel = coarsify.make_me_coarse(data_rel, lat)
    data_tot = coarsify.make_me_coarse(data_tot, lat)

    data_occ = data_rel/data_tot

    return data_occ

def merge_data(data, index):
    low_hit = np.argwhere(index < -0.65)
    high_hit = np.argwhere(index > 0.3)

    lstep = np.arange(0, low_hit.shape[0])
    for l in lstep:
        l_data_hold = data[:,:,int(low_hit[l])]
        try:
            l_data = np.append(l_data, np.array([l_data_hold]), axis=0)
        except UnboundLocalError:
            l_data = np.array([l_data_hold])
    
    hstep = np.arange(0, high_hit.shape[0])
    for h in hstep:
        h_data_hold = data[:,:,int(high_hit[h])]
        try:
            h_data = np.append(h_data, np.array([h_data_hold]), axis=0)
        except UnboundLocalError:
            h_data = np.array([h_data_hold])

    l_data_avg = np.sum(l_data, axis=0)/3
    h_data_avg = np.sum(h_data, axis=0)/3

    return l_data_avg, h_data_avg

#def t_test(data, index, ddof):
def add_zeros(num_pass, n_zeros):
    new_num = num_pass.zfill(n_zeros + len(num_pass))
    
    return new_num

def p2c(data):
    rho = np.linalg.norm(data)
    print(rho.shape)
    phi = np.arctan2(data[:,0],data[0,:])
    print(phi.shape)

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    print(rho.shape)
    print(phi.shape)
    exit()
    return(rho, phi)
  
# def make_rolling_corr(data1, data2):
# #* Data must bee time series of the same length
#     tstep = np.arange(0, data1.shape[0], 5)
#     print(tstep)
#     for t in tstep:
#         if t == 0:
#             continue
#         corr_hold = np.corrcoef(data1[t-5:t], data2[t-5:t])
#         corr = corr_hold[0,1]
#         if t == 5:
#             corr_arr = np.array([corr])
#             continue
#         corr_arr = np.append(corr_arr, corr)

#     return corr_arr

def make_rolling_corr(data1, data2):
    tstep = np.arange(0, data1.shape[0])
    for t in tstep:
        # if t == 0:
        #     continue
        corr_hold = np.corrcoef(data1[0:t], data2[0:t])
        corr = corr_hold[0,1]
        if t == 0:
            corr_arr = np.array([corr])
            continue
        corr_arr = np.append(corr_arr, corr)

    return corr_arr

def make_np_list(list_send):
    list_send = list_send.tolist()
    list_str = map(str, list_send)
    list_fin = list(list_str)

    return list_fin    

def fix_days(days):
    tstep = np.arange(0,len(days))
    for t in tstep:
        day = days[t]
        if int(day) < 10:
            print(day)
            days[t] = add_zeros(day, 1)
            
    return days

def latitudinal_weight(lats, data_s, weight_select):
    if weight_select == 'cos':
        xstep = np.arange(data_s.shape[1])
        for x in xstep:
            lat_app = (np.cos(math.radians(-1*lats[x])))**0.5
            data_app = data_s[:,x]*lat_app
            if x == 0:
                data_arr = np.array([data_app])
                continue
            data_arr = np.append(data_arr, np.array([data_app]), axis=0)
        data_arr = np.rot90(data_arr, k=-1)
        data_arr = np.fliplr(data_arr)

        return data_arr
    
def latitudinal_weight3d(lats, data_s, weight_select):
    #*This doesn't work right now.
    if weight_select == 'cos':
        ystep = np.arange(data_s.shape[1])
        for y in ystep:
            lat_app = (np.cos(math.radians(-1*lats[y])))**0.5
            data_app = data_s[:,y,:]*lat_app
            if y == 0:
                data_arr = np.array([data_app])
                continue
            data_arr = np.append(data_arr, np.array([data_app]), axis=0)
        data_arr = np.rot90(data_arr, k=-1)
        data_arr = np.fliplr(data_arr)

        return data_arr
        
def latitudinal_weight2(lats, data_s, weight_select):
    if weight_select == 'cos':
        xstep = np.arange(data_s.shape[1])
        for x in xstep:
            lat_app = (np.cos(math.radians(-1*lats[x])))**0.5
            #data_app = data_s[:,x]*lat_app
            data_app_h = np.ones(data_s[:,x].shape[0])
            data_app = data_app_h*lat_app
            if x == 0:
                data_arr = np.array([data_app])
                continue
            data_arr = np.append(data_arr, np.array([data_app]), axis=0)
        data_arr = np.rot90(data_arr, k=-1)
        data_arr = np.fliplr(data_arr)

        return data_arr
    
def make_weight(rho_col_send):
    sum_rho = np.sum(rho_col_send)
    rho_col = rho_col_send/sum_rho

    return rho_col

def mass_weight(levels, data_s, data_t):
    ystep = np.arange(0,data_s.shape[0])
    xstep = np.arange(0,data_s.shape[1])
    for x in xstep:
        for y in ystep:
            temp = data_t[y,x]
            pres = levels[y]*100
            r = 286
            rho_0 = 9.8
            if temp == 0:
                temp = 0.0000001
            rho = pres/(temp*r)
            # print(data_s[y,x])
            # print(temp)
            # print(pres)
            # print(rho)
            #w = (rho_0-rho)/rho_0
            if y == 0:
                y_arr = np.array([rho])
                continue
            y_arr = np.append(y_arr, np.array([rho]), axis=0)
        y_arr = make_weight(y_arr)
        if x == 0:
            w_arr = np.array([y_arr])
            continue
        w_arr = np.append(w_arr, np.array([y_arr]), axis=0)
        # print('sum!!!')
        # print(np.sum(y_arr))
    w_arr = np.rot90(w_arr, k=-1)
    w_arr = np.fliplr(w_arr)

    data_w = data_s*w_arr

    return data_w

def mass_weight2(levels, data_s, data_t):
    ystep = np.arange(0,data_s.shape[0])
    xstep = np.arange(0,data_s.shape[1])
    for x in xstep:
        for y in ystep:
            temp = data_t[y,x]
            pres = levels[y]*100
            r = 286
            rho_0 = 9.8
            if temp == 0:
                temp = 0.0000001
            rho = pres/(temp*r)
            # print(data_s[y,x])
            # print(temp)
            # print(pres)
            # print(rho)
            #w = (rho_0-rho)/rho_0
            if y == 0:
                y_arr = np.array([rho])
                continue
            y_arr = np.append(y_arr, np.array([rho]), axis=0)
        y_arr = make_weight(y_arr)
        if x == 0:
            w_arr = np.array([y_arr])
            continue
        w_arr = np.append(w_arr, np.array([y_arr]), axis=0)
        # print('sum!!!')
        # print(np.sum(y_arr))
    w_arr = np.rot90(w_arr, k=-1)
    w_arr = np.fliplr(w_arr)

    data_w = data_s*w_arr

    return w_arr
        
def mass_weight_3(level_s, data_s):
    xstep = np.arange(0,data_s.shape[1])
    ystep = np.arange(0,data_s.shape[0])
    for x in xstep:
        for y in ystep:
            mass_w_hold = (level_s[y]**0.5)/(1000**0.5)

            if y == 0:
                y_arr = np.array([mass_w_hold])
                continue
            y_arr = np.append(y_arr, np.array([mass_w_hold]), axis=0)
        if x == 0:
            w_arr = np.array([y_arr])
            continue
        w_arr = np.append(w_arr, np.array([y_arr]), axis=0)

    w_arr = np.rot90(w_arr, k=-1)
    w_arr = np.fliplr(w_arr)
    return w_arr

def make_monthly_mean_without_time_dem(bam_send, day_step, n_years, n_months):
    """
    This function was originginally designed to take the BAM index (pc1 of EKE) which has data
    at the time interval specified by the code that made the data (make_raw_bam_anim.py for testing)
    and it will make it into monthly mean data. This function will call the month_id function in indexes.py
    Input:
        1. Raw bam index, 2. Day step, 3. Number of years, 4. Number of months
    Output:
        monthly mean bam index, mast_arr. Mast_arr is an array that includes the values and the months.
    """
    day_arr = np.arange(1,32,day_step)
    dt_arr = indexes.make_dt_arr(day_arr, n_years, n_months)
    if dt_arr.shape[0] != bam_send.shape[0]:
        raise Exception('Datetime array shape does not equal bam array shape')

    tstep = np.arange(0,bam_send.shape[0])
    index_start = 0
    for t in tstep:
        if t == 0:
            month_test = dt_arr[t].month
            continue
        month = dt_arr[t].month
        print(month)
        if month != month_test:
            print('!!')
            print(month_test)
            print(month)
            m_data_hold = np.mean(bam_send[index_start:t])
            index_start = t

            try:
                m_data_arr = np.append(m_data_arr, m_data_hold)
                month_arr = np.append(month_arr, month_test)
            except UnboundLocalError:
                m_data_arr = np.array([m_data_hold])
                month_arr = np.array([month_test])
            month_test = month

        #Needed to include this for the end of the loop (there is no next month so month == month_test)
        elif t == tstep[-1]:
            h_data_hold = np.mean(bam_send[index_start:])
            m_data_arr = np.append(m_data_arr, m_data_hold)
            month_arr = np.append(month_arr, month_test)

    mast_arr = np.array([m_data_arr, month_arr])
        
    return m_data_arr, mast_arr

def make_monthly_mean(bam_send, time, bad_months):
    """
    This function was originginally designed to take the BAM index (pc1 of EKE) which has data
    at the time interval specified by the code that made the data (make_raw_bam_anim.py for testing)
    and it will make it into monthly mean data. This function will call the month_id function in indexes.py
    This is a copy, after I made the larger dataset with the accurate time dimension.
    Input:
        1. Raw bam index, 2. Day step, 3. Number of years, 4. Number of months
    Output:
        monthly mean bam index, mast_arr. Mast_arr is an array that includes the values and the months.
   
    Notes:
        Not sure I've managed to get the bad_months integrated in the best way.
    """
    # day_arr = np.arange(1,32,day_step)
    # dt_arr = indexes.make_dt_arr(day_arr, n_years, n_months)
    # if dt_arr.shape[0] != bam_send.shape[0]:
    #     raise Exception('Datetime array shape does not equal bam array shape')
    x = 0
    tstep = np.arange(0,bam_send.shape[0])
    index_start = 0
    for t in tstep:
        if t == 0:
            month_test = time[t].month
            continue
        month = time[t].month
        year = time[t].year
        print(month)      

        if month != month_test:
            print('!!')
            print(month_test)
            print(month)
            m_data_hold = np.mean(bam_send[index_start:t])
            index_start = t
            if x < bad_months.shape[0]:
                if month_test == bad_months[x].month and year == bad_months[x].year:
                    print('Omitting: '+ dt.datetime.strftime(bad_months[x], '%m-%d-%Y'))
                    month_test = month
                    x = x+1
                    continue
            try:
                m_data_arr = np.append(m_data_arr, m_data_hold)
                month_arr = np.append(month_arr, time[t-1])
            except UnboundLocalError:
                m_data_arr = np.array([m_data_hold])
                month_arr = np.array([time[t-1]])
            month_test = month
            print(m_data_arr.shape)

        #Needed to include this for the end of the loop (there is no next month so month == month_test)
        elif t == tstep[-1]:
            h_data_hold = np.mean(bam_send[index_start:])
            m_data_arr = np.append(m_data_arr, m_data_hold)
            month_arr = np.append(month_arr, time[t])

    mast_arr = np.array([m_data_arr, month_arr])
        
    return m_data_arr, mast_arr

def make_monthly_cs_means(dz_top, dz_bot, z_top, z_bot):
    """
    This function will use diveinnewdata.py to build monthly cloud occurrence means.
    """
    month_arr = np.arange(1,13)
    months = make_np_list(month_arr)
    months = fix_days(months)

    for month in months:
        data_arr = diveinnewdata.find_files(month, 0)
        cld_cnt_raw = data_arr[:,6]
        cld_tot_raw = data_arr[:,7]
        lat = data_arr[0,0]
        lon = data_arr[0,1]
        dz_bin = data_arr[0,4]
        dz_dbin = data_arr[0,5]
        cbtz_bin = data_arr[0,2]
        cbtz_dbin = data_arr[0,3]

        restrict_arr = np.array([dz_bin, dz_top, dz_bot, cbtz_bin, z_top, z_bot, dz_dbin, cbtz_dbin], dtype=object)
        go_stats_arr = jja_global_avg.go_stats_func(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr)
        cld_occ = go_stats_arr[0]
        lat_c = go_stats_arr[1]
        lon_c = go_stats_arr[2]

        if month == months[0]:
            monthly_mean_arr = np.array([cld_occ])
            continue
        monthly_mean_arr = np.append(monthly_mean_arr, np.array([cld_occ]), axis=0)

    return monthly_mean_arr

def make_monthly_zonal_cs_means():
    """
    This function is very similar to make_monthly_cs_means, but it does the zonal data.
    """
    month_arr = np.arange(1,13)
    months = make_np_list(month_arr)
    months = fix_days(months)
    for month in months:
        data_arr = get_data.find_files(month, 0)
        lat = data_arr[0,1]
        lon = data_arr[0,2]
        z = data_arr[0,3]
        dbz = data_arr[0,4]

        cld_cnt_raw = data_arr[:,-1]
        cld_occ = build_stats(data_arr, 1)

        if month == months[0]:
            monthly_mean_arr = np.array([cld_occ])
            continue
        monthly_mean_arr = np.append(monthly_mean_arr, np.array([cld_occ]), axis=0)

    return monthly_mean_arr

def find_lat_lon(lat_s, lon_s, lat_sel, lon_sel):
    """
    This function will return the index values of specific lat/lon
    """
    min_arr_lat = abs(lat_s-lat_sel)
    hit_arr_lat = np.argwhere(min_arr_lat == np.min(min_arr_lat))
    hit_inx_lat = int(hit_arr_lat[0])

    min_arr_lon = abs(lon_s-lon_sel)
    hit_arr_lon = np.argwhere(min_arr_lon == np.min(min_arr_lon))
    hit_inx_lon = int(hit_arr_lon[0])

    return hit_arr_lat, hit_arr_lon

def remove_mean(data, dz_top,dz_bot,z_top,z_bot, zonal_select):
    """
    This function removes the monthly cloudsat cloud occurrence mean from the data.
    Returns: data_arr 
    """
    if zonal_select == 1:
        monthly_means = make_monthly_zonal_cs_means()
    else:
        monthly_means = make_monthly_cs_means(dz_top,dz_bot,z_top,z_bot)
    tstep = np.arange(0, data.shape[0])
    for t in tstep:
        month_data = (t+5)%12

        month_mean = monthly_means[month_data-1,:,:]
        cld_occ_anom_h = data[t,:,:] - month_mean
        
        if t == 0:
            cld_occ_arr = np.array([cld_occ_anom_h])
            continue
        cld_occ_arr = np.append(cld_occ_arr, np.array([cld_occ_anom_h]), axis=0)

    #print(cld_occ_arr.shape)
    return cld_occ_arr

def apply_sin_model(model_s, data_s):
    """
    This function was made for applying a sin model to the seasonality of the bam to try to
    remove it. It will return a time series of values that the model predicts.
    """
    c1 = model_s[0]
    c2 = model_s[1]
    c3 = model_s[2]
    c4 = model_s[3]
    c5 = model_s[4]
    #c6 = model_s[5]
    dstep = np.arange(0,data_s.shape[0])
    for d in dstep:
        y = (d**4)*c1 + (d**3)*c2 + (d**2)*c3 + d*c4 + c5
        #y = (d**5)*c1 + (d**4)*c2 + (d**3)*c3 + (d**2)*c4 + c5*d + c6
        if d == 0:
            y_arr = np.array([y])
            continue
        y_arr = np.append(y_arr, y)
    
    return y_arr

def make_daily_mean_vals_test(bam_data_s):
    """This function will loop through the "daily" bam_index and remove the mean to hopefully
    remove the seasonality of the data
    The fact that Feburary doesn't have a 29 is messing it up.
    """
    tstep = np.arange(0,97)
    ystep = np.arange(3)
    for t in tstep:
        if t == 16:
            continue
        print(t)
        for y in ystep:
            d_hold = bam_data_s[t+(y*96)]
            if y == 0:
                y_arr = np.array([d_hold])
                continue
            y_arr = np.append(y_arr, d_hold)
        d_mean = np.mean(y_arr)
        if t == 0:
            d_mean_arr = np.array(d_mean)
            continue
        d_mean_arr = np.append(d_mean_arr, d_mean)

    return d_mean_arr

def make_daily_mean_vals(bam_data_s, day_skip, years):
    """This function will loop through the "daily" bam_index and remove the mean to hopefully
    remove the seasonality of the data
    The fact that Feburary doesn't have a 29 is messing it up.
    """
    tstep = np.arange(0,365,day_skip)
    tstep = tstep/day_skip
    ystep = np.arange(0,years)
    for t in tstep:
        if t == 16:
            continue
        print(t)
        for y in ystep:
            print(y)
            print(int(y)*int(tstep.shape[0]))
            d_hold = bam_data_s[int(t)+(int(y)*int(tstep.shape[0]))]
            if y == 0:
                y_arr = np.array([d_hold])
                continue
            y_arr = np.append(y_arr, d_hold)
        d_mean = np.mean(y_arr)
        if t == 0:
            d_mean_arr = np.array(d_mean)
            continue
        d_mean_arr = np.append(d_mean_arr, d_mean)

    return d_mean_arr

def remove_mean_bam(data_send, daily_mean):
    """
    This function will take the daily mean values and subtract them from the bam index. Hopefully removing
    the seasonality
    """
    tstep = np.arange(0,data_send.shape[0])
    for t in tstep:
        m = t%daily_mean.shape[0]
        print(m)
        d_anom = data_send[t]-daily_mean[m]
        if t == 0:
            d_anom_arr = np.array([d_anom])
            continue
        d_anom_arr = np.append(d_anom_arr, d_anom)

    return d_anom_arr

def cmpt_cross_cov_mtx(d1_s, d2_s):
    """
    This function will compute the cross covariance matrix of two values. The values need to be sent as
    flattened matrices, or matrixes with shape: (time, (lat*lon)). The output will be a matrix with shape: (N1*N2)
    where N = lat*lon observations.
    """
    xstep = np.arange(0,d1_s.shape[-1])
    ystep = np.arange(0,d2_s.shape[-1])
    if d1_s.shape[-1] != d2_s.shape[-1]:
        raise Exception('inputs do not have same shape')
    for x in xstep:
        print('!!!!x'+str(x))
        for y in ystep:
            cov_hold = (d1_s[:,x]*d2_s[:,y])/d1_s.shape[0]
            if y == 0:
                y_arr = np.array([cov_hold])
                continue
            y_arr = np.append(y_arr, np.array([cov_hold]), axis=0)
        if x == 0:
            cov_arr = np.array([y_arr])
            continue
        cov_arr = np.append(cov_arr, np.array([y_arr]), axis=0)

    return cov_arr

def corr_test(std_err, std_err_ctoff):
    """
    This function will test the BAM correlation values to see if they should be included or exluded.
    Likely could be used for other things. The number of 0 values are creating  misleading correlation 
    coefficient values.
    Currently using standard error as the determinate.
    """

def test_cld_occ(data_s, dates_s):
    """
    This function will loop through the cloud occurrence data and remove years and months that
    have a mean cloud occurrence value of 0 (no data).
    Returns:
        data_s without bad years/months
        data_s_remove: the years and months it took out.
    """
    tstep = np.arange(0,data_s.shape[0])
    for t in tstep:
        date_hold = dates_s[t]
        data_hold = data_s[t,:,:]
        test = np.mean(data_hold)
        if test == 0.0:
            try:
                bd_data = np.append(bd_data, np.array([date_hold]), axis=0)
            except UnboundLocalError:
                bd_data = np.array([date_hold])
            continue
        if t == 0:
            cld_occ = np.array([data_hold])
            continue
        cld_occ = np.append(cld_occ, np.array([data_hold]), axis=0)

    return cld_occ, bd_data

if __name__ == '__main__':
    make_daily_mean_vals(np.array([1, 2,3,4]), 32, 3)

if __name__ == '__main__':
    #make_rolling_corr(np.arange(1,20), np.arange(1,20))
    remove_mean(np.arange(0,43), 1500,0,1000,0)
    #make_monthly_cs_means(1500,0,1000,0)
    exit()
    data_get_dict = get_data.find_files('08', 0)
    data_occ = build_stats(data_get_dict, 2)
    mast_plot.make_zonal_mean_plot(data_occ, data_get_dict[0,1], data_get_dict[0,3], 'More testing 08', 'imshow_test_08', '200000')
