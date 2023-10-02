"""
*This script is where some of the BAM EOF analysis will take place.

Script notes:
    The first 3 months of 2012 and 2-2016 don't have data, however, they return 0.0 all across the globe.
    This skews the anomalies to negative values.

To do:
    Clean up main call. 
    Clean up integration with eof_analysis
    Make a better monthly mean maker function.
    I still think the seasonality remove that I've made doesn't work that well.

Leave notes:
    Something is still weird with the make monthly anom function.

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
from statsmodels.tsa.seasonal import seasonal_decompose

def build_cov_arr_1st_dem(data, corsel, index, lat, lon, lat_sel, lon_sel, std_err_ctoff):
    # model, intercept, rval, pval, stderr = sp.linregress(index, data[:,12,150])
    # cov_m_hold = np.corrcoef(index, data[:,12,150])
    # print(lat[lat_get],lon[lon_get])
    # mast_plot.one_point_anal_bsc(index, data[:,y,x], model, intercept, cov_m_hold[0,1], pval, lat[lat_get], lon[lon_get])
    # mast_plot.make_ts_bam_cld_occ(index, data[:,y,x], model, intercept, cov_m_hold[0,1], pval, lat[lat_get], lon[lon_get])
    # exit()
#* This is the same code as oned_corr.build_cov_arr(). I just needed to rework it because for the 
#cloudsat data, the time dimension is the last dimension, but for this data it is the first dimension.
    xstep = np.arange(0,data.shape[-1])
    ystep = np.arange(0,data.shape[1])
    print(index)
    lat_get, lon_get = go_stats.find_lat_lon(lat, lon, lat_sel, lon_sel)
    if lat_get.shape[0] > 1:
        lat_get = lat_get[0]
    if lon_get.shape[0] > 1:
        lon_get = lon_get[0]
    print(lat_get, lon_get)
    for y in ystep:
        for x in xstep:
            if corsel == 1:
                cov_m_hold = np.correlate(data[:,y,x], index)
                cov_hold = cov_m_hold
            elif corsel == 2:
                model, intercept, rval, pval, stderr = sp.linregress(index, data[:,y,x])
                cov_m_hold = np.corrcoef(index, data[:,y,x])
                # print(y)
                # print(lat_get)
                # if y == lat_get and x == lon_get:
                #     print(lat[lat_get],lon[lon_get])
                #     mast_plot.one_point_anal_bsc(index, data[:,y,x], model, intercept, cov_m_hold[0,1], pval, lat[lat_get], lon[lon_get])
                #     mast_plot.make_ts_bam_cld_occ(index, data[:,y,x], model, intercept, cov_m_hold[0,1], pval, lat[lat_get], lon[lon_get])
                #     exit()
                cov_hold = cov_m_hold[0,1]
                # if stderr >= std_err_ctoff:
                #     cov_hold = 0.0
            elif corsel == 'reg_itc':
                model, intercept, rval, pval, stderr = sp.linregress(index, data[:,y,x])
                cov_hold = intercept
                # if stderr >= std_err_ctoff:
                #     cov_hold = 0.0
            elif corsel == 'reg_slp':
                model, intercept, rval, pval, stderr = sp.linregress(index, data[:,y,x])
                cov_hold = model
                # if stderr >= std_err_ctoff:
                #     cov_hold = 0.0
            elif corsel == 'reg_rval':
                model, intercept, rval, pval, stderr = sp.linregress(index, data[:,y,x])
                cov_hold = rval**2
            elif corsel == 'reg_pval':
                model, intercept, rval, pval, stderr = sp.linregress(index, data[:,y,x])
                cov_hold = pval
            elif corsel == 'reg_stderr':
                model, intercept, rval, pval, stderr = sp.linregress(index, data[:,y,x])
                cov_hold = stderr
            elif corsel == 'onep_anal':
                model, intercept, rval, pval, stderr = sp.linregress(index, data[:,y,x])
                cov_m_hold = np.corrcoef(data[:,y,x], index)
                cov_hold = cov_m_hold[0,1]
                print(lat[y], lon[x])
                mast_plot.one_point_anal(index, data[:,y,x], model, intercept, lat[y], lon[x], pval, cov_hold)
            else:
                cov_m_hold = np.cov(data[:,y,x], index)
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

def sam_test_monthly_mean(pc1):
#*This will take the daily data and make it into a mean to compare it to the sam index
    dec_days = 31
    feb_days = 28
    tstep = np.arange(0,pc1.shape[0])
    dec_mean = np.mean(pc1[:dec_days])
    mean_arr = np.array([dec_mean])
    jan_mean = np.mean(pc1[dec_days:dec_days*2])
    mean_arr = np.append(mean_arr, np.array([jan_mean]))
    feb_mean = np.mean(pc1[dec_days*2:dec_days*2+feb_days])
    mean_arr = np.append(mean_arr, np.array([feb_mean]))

    year1 = dec_days*2+feb_days
    dec2_mean = np.mean(pc1[year1:year1+dec_days])
    mean_arr = np.append(mean_arr, np.array([dec2_mean]))
    jan2_mean = np.mean(pc1[year1+dec_days:year1+dec_days*2])
    mean_arr = np.append(mean_arr, np.array([jan2_mean]))
    feb2_mean = np.mean(pc1[year1+dec_days*2:])
    mean_arr = np.append(mean_arr, np.array([feb2_mean]))

    print(mean_arr.shape)
    return mean_arr

def project(pc, data_send):
    data_reshape = data_send.reshape(data_send.shape[0], -1)
    pc_mean = np.mean(pc)
    pc_anom = pc - pc_mean

    projection = np.dot(data_reshape.T, pc_anom)
    
    twod_projection = projection.reshape(42,100)

    return twod_projection

def projection_arr(eof1, data_send):
    tstep = np.arange(0,data_send.shape[0])
    for t in tstep:
        data = data_send[t,:,:]*eof1
        if t==0:
            projection_arr = np.array([data])
            continue
        projection_arr = np.append(projection_arr, np.array([data]), axis=0)

    return projection_arr

def eof_ts(eof1, pc1):
    tstep = np.arange(0,pc1.shape[0])
    for t in tstep:
        data = eof1*pc1[t]
        if t==0:
            projection_arr = np.array([data])
            continue
        projection_arr = np.append(projection_arr, np.array([data]), axis=0)

    return projection_arr

def make_eke_filter_test(eke_data_s, lat_s, levels_s, temp_s):
    tstep = np.arange(eke_data_s.shape[0])
    for t in tstep:
        #eke_data_hold = go_stats.mass_weight(levels_s, eke_data_s[t,:,:], temp_s)
        if t == 0:
            mass_weight_arr = go_stats.mass_weight2(levels_s, eke_data_s[t,:,:], temp_s)
            lat_weight_arr = go_stats.latitudinal_weight2(lat_s, eke_data_s[t,:,:], 'cos')
        if t == 0:
            eke_data = np.array([eke_data_s[t,:,:]*mass_weight_arr*lat_weight_arr])
            continue
        eke_data = np.append(eke_data, np.array([eke_data_s[t,:,:]*mass_weight_arr*lat_weight_arr]), axis=0)

    return eke_data

def run_main_bam(bd_months):
    """
    This is a callable function that will return the bam_index in monthly mean form.
    Returns:
        EOF1 of eke, PC1 of EKE (monthly mean)&zscored.
    """
    #level, lat_b, lon_b, u_anom, u3d, v_anom, v3d, eke_data, eke3d = import_eke_data.import_data3d('u_anom_test3d_2006_2009.nc4', 0)
    level, lat_b, lon_b, time, u_anom, u3d, v_anom, v3d, eke_data, eke3d = import_eke_data.merge_data_3d()

    #Not sure why I do this.
    eke3d = eke3d*2
    eke3d_m = np.mean(eke3d, axis=1)
    eke3d_m = eke3d_m/2
    u3d_m = np.mean(u3d, axis=1)
    mass_weights = go_stats.mass_weight_3(level, eke3d[10,5,:,:])
    eke3d_w, uanom_w = make_raw_bam_anim.apply_weights(eke3d_m, u3d_m, mass_weights)

    eofs, pcs, explained_var, pres_shape, lat_shape = eof_analysis.eof_main(eke3d_w, 6)

    bam_index = sp.zscore(pcs[:,0])
    bam_index_nsr, mast_arr_nsr = go_stats.make_monthly_mean(bam_index, time, bd_months)

    dailymean = go_stats.make_daily_mean_vals_test(bam_index)
    d_anom_arr = go_stats.remove_mean_bam(bam_index, dailymean)
    
    #!!Took out for NAO comparision
    bam_index_mm, mast_arr = go_stats.make_monthly_mean(d_anom_arr, time, bd_months)
    bam_index_mm_send = bam_index_mm[4:-1]

    #bam_index_mm_send = go_stats.make_monthly_mean(bam_index, 4, 4, 12)
  
    eof1 = eofs[0].reshape(pres_shape, lat_shape)

    return eof1, mast_arr_nsr


if __name__ == '__main__':
    years_arr = np.arange(2006,2017)
    hit2011 = np.where(years_arr==2011)
    years_arr = np.delete(years_arr, hit2011)
    years = go_stats.make_np_list(years_arr)
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']

    cld_occ_dates = indexes.make_dt_arr(np.arange(1,2),11,12,'01-01-2006',2011)
    cld_occ_dates = cld_occ_dates[5:]
    cld_occ, lat, lon = oned_corr.eof_call_main(years, months, 0, 0)
    cld_occ_gd, bd_data = go_stats.test_cld_occ(cld_occ, cld_occ_dates)

    cld_occ_anom = go_stats.remove_mean(cld_occ_gd, 1500,0,1000,0,0)

    # cld_anom_mean = np.mean(cld_occ_anom, axis=1)
    # cld_anom_mean = np.mean(cld_anom_mean, axis=1)
    # cld_occ_mean = np.mean(cld_occ, axis=1)
    # cld_occ_mean = np.mean(cld_occ_mean, axis=1)

    # pr`int(cld_occ_mean.shape)
    # print(cld_anom_mean.shape)

    # mast_plot.seasonal_cycle_test(cld_anom_mean, cld_occ_mean)
    # exit()

    #level, lat_b, lon_b, time, u_anom, u3d, v_anom, v3d, eke_data, eke3d = import_eke_data.import_data3d('u_anom_test3d_2010_2016.nc4', 1)
    level, lat_b, lon_b, time, u_anom, u3d, v_anom, v3d, eke_data, eke3d = import_eke_data.merge_data_3d()

    eke3d = eke3d*2
    eke3d_m = np.mean(eke3d, axis=1)
    eke3d_m = eke3d_m/2
    u3d_m = np.mean(u3d, axis=1)
    mass_weights = go_stats.mass_weight_3(level, eke3d[10,5,:,:])
    eke3d_w, uanom_w = make_raw_bam_anim.apply_weights(eke3d_m, u3d_m, mass_weights)

    eofs, pcs, explained_var, pres_shape, lat_shape = eof_analysis.eof_main(eke3d_w, 6)

    bam_index = sp.zscore(pcs[:,0], ddof=1)
    print(bam_index.shape)

    #bam_index_raw_mm, mast_arr_raw = go_stats.make_monthly_mean(bam_index, time, bd_data)
    #bam_index_raw_mm_send = bam_index_raw_mm[5:]

    dailymean = go_stats.make_daily_mean_vals(bam_index, 4, 3)
    dailymean_test = go_stats.make_daily_mean_vals_test(bam_index)

    d_anom_arr = go_stats.remove_mean_bam(bam_index, dailymean)
    d_anom_arr_test = go_stats.remove_mean_bam(bam_index, dailymean_test)
    print(d_anom_arr_test.shape)

    #*Seasonality reduction function 2
    # bam_index_mm = go_stats.make_monthly_mean(d_anom_arr, 4, 4, 12)
    # bam_index_mm_send = bam_index_mm[7:]

    #*Seasonality reduction function 1
    bam_index_mm, mast_arr_anom = go_stats.make_monthly_mean(d_anom_arr_test, time, bd_data)
    #print(mast_arr_anom)


    bam_index_mm_send = bam_index_mm[4:-1]
    cld_occ_anom_send = cld_occ_anom[:]

    eof1 = eofs[0].reshape(pres_shape, lat_shape)

    #*This code is for making the EOF projected on eke.
    # eke_corr = build_cov_arr_1st_dem(eke3d_w, 2, bam_index, level, lat_b, 800,-40, 0.4)
    # eke_reg_inc = build_cov_arr_1st_dem(eke3d_w, 'reg_itc', bam_index, level, lat_b, 800,-40, 0.4)
    # eke_reg_slp = build_cov_arr_1st_dem(eke3d_w, 'reg_slp', bam_index, level, lat_b, 800,-40, 0.4)

    # mast_plot.bam_eke_corr(eke_corr, eke_reg_inc, eke_reg_slp)
    # exit()

    print('Corr vals!!')
    print(bam_index_mm_send.shape)
    print(cld_occ.shape)

    #test_corr = build_cov_arr_1st_dem(eke3d_w, 'reg_itc', bam_index)
    #onep_anal = build_cov_arr_1st_dem(cld_occ_anom_send, 'onep_anal', bam_index_mm_send, lat, lon, -90, -90)
    # exit()
    cld_corr = build_cov_arr_1st_dem(cld_occ_anom_send, 2, bam_index_mm_send, lat, lon, 41, 29, 1)
    cld_reg_int = build_cov_arr_1st_dem(cld_occ_anom_send, 'reg_itc', bam_index_mm_send, lat, lon, -41, -18, 1)
    cld_reg_slp = build_cov_arr_1st_dem(cld_occ_anom_send, 'reg_slp', bam_index_mm_send, lat, lon, -41, 18, 1)
    # rval_arr = build_cov_arr_1st_dem(cld_occ_anom, 'reg_rval', bam_index_mm_send, lat, lon, 45, 21)
    # pval_arr = build_cov_arr_1st_dem(cld_occ_anom, 'reg_pval', bam_index_mm_send, lat, lon, 45, 21)

    # mast_plot.bam_stats_package(pval_arr, rval_arr)
    # exit()

    print(np.min(cld_reg_int))
    #cld_corr = np.flipud(cld_corr)

    mast_plot.bam_cld_corr(cld_corr, cld_reg_int, cld_reg_slp)
    exit()

    lat_lon_bounds = [-20,-90,180,-180]
    bounds_arr = jja_global_avg.restrict_domain(lat, lon, lat_lon_bounds[0],lat_lon_bounds[1], lat_lon_bounds[2], lat_lon_bounds[3],
                                                lat_lon_bounds[-1])

    exit()

    #norm = mast_plot.norm_maker(cld_corr)
    p = plt.imshow(cld_corr, cmap = 'RdBu')
    plt.colorbar(p)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')
    


    

