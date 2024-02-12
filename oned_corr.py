"""
This code builds the 1 grid spacing correlation with the ENSO index.
If subdomain is selected, change domain call to 1 & enter latlon_bounds.
For other correlation codes that use this, the same timeframe must be selected.
**Notes on go_stats_func & go_stats_func_wm. This function will get the statistics associated 
with occurrence from the raw data arrays. go_stats_func will average the occurrence over the season
that is selected (annual average). go_stats_func_wm will get the statistics for each month and return an array
of as many months as are in the season. **To switch between the two functions just change the function call
on ~225/273 depending on where you are running the code. ALSO, the append method will need to be changed at the
bottom of the loop.

Leave notes: Pearson isn't working well. Not sure why there is correlation over the continents. 
    Should probably look at the raw data for each year and see if there are some years with cloud occurrence over
    the continents. 

To do:
    Need to fix the statistical significance stipling, if there is no statistical significance it errors out.
    Should make the index sorting non-manual.
    

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
import datetime as dt
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
import bam_analysis

#months = ['09','10','11']
# months_arr = np.arange(1,13)
# months = go_stats.make_np_list(months_arr)
# months = go_stats.fix_days(months)
# index_months = ['09','10','11']
# years = ['2006','2007','2008','2009','2010','2012','2013','2014','2015','2016']
# high_enso_years = ['2006','2009','2012','2015']
# low_enso_years = ['2013','2014','2016']
# mean = ['2006','2007','2008','2009','2010','2012','2013','2014','2015','2016']
z_bot = 0
z_top = 1000
dz_bot = 0
dz_top = 1500

#index = np.array([0.4, -0.7, -0.95, 0.55, -2.5, 0.2, -0.65, 0.3, 1.8, -0.35])
#index = indexes.enso(months)
#0.8 was taken out of the index array because 2006 was taken out above
#-0.3 was taken out of the end index array because 2016 was taken out above
#index = np.array([0.8, -1.2, -1.1, 1.2, -2.1, -0.2, -0.2, 0.3, 1.9, -0.3])

restrict_domain_call = 0
#**Current domain options**
#SEP
# bounds_arr = jja_global_avg.restrict_domain(lat, lon, -5, -40, -100, -70)
#latlon_bounds = [-5,-40,-100,-70]
#NEP
# bounds_arr = jja_global_avg.restrict_domain(lat, lon, 37, 10, -140, -116)
#latlon_bounds = [37,10,-140,-116]
#NEA
#latlon_bounds = [40,20,-30,-10]
#SOP
#latlon_bounds = [-43,-64,-140,-80]
#SA
#latlon_bounds = [-43,-64,-50,2]
#TEST
#latlon_bounds = [0,-20,-65,-55]
#Southern ocean projection
#latlon_bounds = [-40,-90,-180,180]
#Longwave test 
latlon_bounds = [48,25,-124,-113]


def go_append(x1, x2, init_flag):
#!!Need to write out how I did this....soon
    if init_flag == 1:
        print(x1.shape)
        ystep = np.arange(0, x1.shape[0])
        xstep = np.arange(0, x1.shape[1])
        for y in ystep:
            for x in xstep:
                if x == 0:
                    hold_x = np.array([[x1[y,x], x2[y,x]]])
                    continue
                hold_x = np.append(hold_x, np.array([[x1[y,x], x2[y,x]]]), axis=0)
            if y == 0:
                hold_y = np.array([hold_x])
                continue
            hold_y = np.append(hold_y, np.array([hold_x]), axis =0)
    else:
        ystep = np.arange(0, x1.shape[0])
        xstep = np.arange(0, x1.shape[1])
        for y in ystep:
            for x in xstep:
                if x == 0:
                    hold_x = np.array([np.insert(x1[y,x,:], x1.shape[-1], x2[y,x])])
                    continue
                hold_x = np.append(hold_x, np.array([np.insert(x1[y,x,:], x1.shape[-1], x2[y,x])]), axis=0)
            if y == 0:
                hold_y = np.array([hold_x])
                continue
            hold_y = np.append(hold_y, np.array([hold_x]), axis =0)

    return hold_y

def build_cov_arr(data, corsel, index):
    xstep = np.arange(0,data.shape[1])
    ystep = np.arange(0,data.shape[0])

    for y in ystep:
        for x in xstep:
            if corsel == 1:
                cov_m_hold = np.correlate(data[y,x,:], index)
                cov_hold = cov_m_hold
            elif corsel == 2:
                cov_m_hold = np.corrcoef(data[y,x,:], index)
                cov_hold = cov_m_hold[0,1]
            else:
                cov_m_hold = np.cov(data[y,x,:], index)
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

def spatial_sig(data, std_2):
    ystep = np.arange(0, data.shape[0])
    xstep = np.arange(0, data.shape[1])
    for y in ystep:
        for x in xstep:
            if abs(data[y,x]) > std_2:
                hold = 1
            else:
                hold = 0
            if x == 0:
                stat_x_arr = np.array([hold])
                continue
            stat_x_arr = np.append(stat_x_arr, np.array([hold]), axis=0)
        if y == 0:
            stat_arr = np.array([stat_x_arr])
            continue
        stat_arr = np.append(stat_arr, np.array([stat_x_arr]), axis=0)

    return stat_arr

def get_lat_lon(data, sig, lat, lon):
    tstep = np.arange(0, sig.shape[0])
    for t in tstep:
        y = int(sig[t][0])
        print(type(y))
        x = int(sig[t][1])
        lat_hold = lat[y]
        lon_hold = lon[x]
        if t == 0:
            lat_arr = np.array([lat_hold])
            lon_arr = np.array([lon_hold])
            continue
        lat_arr = np.append(lat_arr, np.array([lat_hold]), axis=0)
        lon_arr = np.append(lon_arr, np.array([lon_hold]), axis=0)

    return lat_arr, lon_arr

def my_pearson(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    n = x.shape[0]+1

    tstep = np.arange(0, x.shape[0])
    for t in tstep:
        x_hold = (x[t]-x_mean)**2
        y_hold = (y[t]-y_mean)**2
        t_hold = (x[t]-x_mean)*(y[t]-y_mean)
        if t == 0:
            x_sum = x_hold
            y_sum = y_hold
            t_sum = t_hold
            continue
        x_sum = x_sum + x_hold
        y_sum = y_sum + y_hold
        t_sum = t_sum + t_hold
    num = t_sum/(n-1)
    denom = ((x_sum/(n-1))**0.5)*((y_sum/(n-1))**0.5)
    r = num/denom
    print('r:')
    print(r)

def eof_call_main(years, months, restrict_domain_call, mix_years, latlon_bounds):
    for year in years:
        data_arr = overview.get_data([year], months, mix_years)
        lat = data_arr[0,0]
        lon = data_arr[0,1]
        cbtz_bin = data_arr[0,2]
        cbtz_dbin = data_arr[0,3]
        dz_bin = data_arr[0,4]
        dz_dbin = data_arr[0,5]

        restrict_arr = np.array([dz_bin, dz_top, dz_bot, cbtz_bin, z_top, z_bot, dz_dbin, cbtz_dbin], dtype=object)

        cld_cnt_raw = data_arr[:,6]
        #print(cld_cnt_raw.shape)
        cld_tot_raw = data_arr[:,-2]

        go_stats_arr = jja_global_avg.go_stats_func(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr)
        cld_occ = go_stats_arr[0]
        lat = go_stats_arr[1][0]
        lon = go_stats_arr[2][0]

        if restrict_domain_call == 1:
            bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0], latlon_bounds[1],
                                                        latlon_bounds[2], latlon_bounds[3])
            print(bounds_arr)
            # if using go stats function wm, use:
            # cld_occ, lat, lon = jja_global_avg.region_arr3d(cld_occ, lat, lon, bounds_arr)
            # else:
            cld_occ, lat, lon = jja_global_avg.region_arr(cld_occ, lat, lon, bounds_arr)
        #*This is the way to append together if using go_stats_func_wm
        # if year == years[0]:
        #     cld_occ_arr = cld_occ
        #     continue
        # cld_occ_arr = np.append(cld_occ_arr, cld_occ, axis=0)

        #*This is the way to append together if using go_stats_func
        if year == years[0]:
            cld_occ_arr = np.array([cld_occ])
            continue
        cld_occ_arr = np.append(cld_occ_arr, np.array([cld_occ]), axis=0)

        #*This was the old way it was done. Not sure if it will work that well when using go_stats_func_wm
        # if year == years[0]:
        #     cld_occ_hold = cld_occ
        #     continue
        # if year == years[1]:
        #     hold_test = go_append(cld_occ_hold, cld_occ, 1)
        #     continue
        # hold_test = go_append(hold_test, cld_occ, 0)

    return cld_occ_arr, lat, lon

def eof_call_main_2(years, months, restrict_domain_call, mix_years, latlon_bounds):
    print(years)
    for year in years:
        data_arr = overview.get_data([year], months, mix_years)
        lat = data_arr[0,0]
        lon = data_arr[0,1]
        cbtz_bin = data_arr[0,2]
        cbtz_dbin = data_arr[0,3]
        dz_bin = data_arr[0,4]
        dz_dbin = data_arr[0,5]

        restrict_arr = np.array([dz_bin, dz_top, dz_bot, cbtz_bin, z_top, z_bot, dz_dbin, cbtz_dbin], dtype=object)

        cld_cnt_raw = data_arr[:,6]
        #print(cld_cnt_raw.shape)
        cld_tot_raw = data_arr[:,-2]

        go_stats_arr = jja_global_avg.go_stats_func(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr)
        cld_occ = go_stats_arr[0]
        lat = go_stats_arr[1][0]
        lon = go_stats_arr[2][0]

        if restrict_domain_call == 1:
            bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0], latlon_bounds[1],
                                                        latlon_bounds[2], latlon_bounds[3])
            print(bounds_arr)
            # if using go stats function wm, use:
            # cld_occ, lat, lon = jja_global_avg.region_arr3d(cld_occ, lat, lon, bounds_arr)
            # else:
            cld_occ, lat, lon = jja_global_avg.region_arr(cld_occ, lat, lon, bounds_arr)
        #*This is the way to append together if using go_stats_func_wm
        # if year == years[0]:
        #     cld_occ_arr = cld_occ
        #     continue
        # cld_occ_arr = np.append(cld_occ_arr, cld_occ, axis=0)

        #*This is the way to append together if using go_stats_func
        if year == years[0]:
            cld_occ_arr = np.array([cld_occ])
            continue
        cld_occ_arr = np.append(cld_occ_arr, np.array([cld_occ]), axis=0)

        #*This was the old way it was done. Not sure if it will work that well when using go_stats_func_wm
        # if year == years[0]:
        #     cld_occ_hold = cld_occ
        #     continue
        # if year == years[1]:
        #     hold_test = go_append(cld_occ_hold, cld_occ, 1)
        #     continue
        # hold_test = go_append(hold_test, cld_occ, 0)

    return cld_occ_arr, lat, lon
            
if __name__ == '__main__':
    months_arr = np.arange(3,6)
    # months_arr = np.array([12,1,2])
    months = go_stats.make_np_list(months_arr)
    months = go_stats.fix_days(months)
    inx_months_arr = np.arange(2,5)
    # inx_months_arr = np.array([11,12,1])
    inx_months = go_stats.make_np_list(inx_months_arr)
    inx_months = go_stats.fix_days(inx_months)
    years_arr = np.arange(2007,2017)
    delete = np.argwhere(years_arr==2011)
    years_arr = np.delete(years_arr, int(delete[0]))
    delete_2 = np.argwhere(years_arr==2012)
    years_arr = np.delete(years_arr, int(delete_2[0]))
    years = go_stats.make_np_list(years_arr)
    season = 'mam'

    # index = indexes.enso_full(inx_months, years_arr, 0)
    index = indexes.sam(inx_months, years_arr, 0)
    for year in years:
        data_arr = overview.get_data([year], months, 0)
        lat = data_arr[0,0]
        lon = data_arr[0,1]
        cbtz_bin = data_arr[0,2]
        cbtz_dbin = data_arr[0,3]
        dz_bin = data_arr[0,4]
        dz_dbin = data_arr[0,5]

        restrict_arr = np.array([dz_bin, dz_top, dz_bot, cbtz_bin, z_top, z_bot, dz_dbin, cbtz_dbin], dtype=object)

        cld_cnt_raw = data_arr[:,6]
        cld_tot_raw = data_arr[:,-2]

        #go_stats_arr = jja_global_avg.go_stats_func_wm(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr)
        go_stats_arr = jja_global_avg.go_stats_func(cld_cnt_raw, cld_tot_raw, lat, lon, restrict_arr)
        cld_occ = go_stats_arr[0]
        lat = go_stats_arr[1][0]
        lon = go_stats_arr[2][0]

        if restrict_domain_call == 1:
            bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0], latlon_bounds[1],
                                                        latlon_bounds[2], latlon_bounds[3])
            cld_occ, lat, lon = jja_global_avg.region_arr(cld_occ, lat, lon, bounds_arr)

        if year == years[0]:
            hold_test = np.array([cld_occ])
            continue
        hold_test = np.append(hold_test, np.array([cld_occ]), axis=0)

        # if year == years[0]:
        #     cld_occ_hold = cld_occ
        #     continue
        # if year == years[1]:
        #     hold_test = go_append(cld_occ_hold, cld_occ, 1)
        #     continue
        # hold_test = go_append(hold_test, cld_occ, 0)
        
    #cov_arr = build_cov_arr(hold_test, 2, index)
    #index = index[5:]
    #cld_occ_anom = go_stats.remove_mean(hold_test, 1500,0,1000,0,0)

    # cld_occ_anom_m = np.mean(cld_occ_anom, axis=1)
    # cld_occ_anom_s = np.mean(cld_occ_anom_m, axis=1)
    # cld_occ_m = np.mean(hold_test, axis=1)
    # cld_occ_s = np.mean(cld_occ_m, axis=1)

    # mast_plot.seasonal_cycle_test(cld_occ_anom_s, cld_occ_s)
    # exit()
    # reg_inc = np.mean(cld_occ_anom, axis=0)
    # cov_arr = bam_analysis.build_cov_arr_1st_dem(cld_occ_anom, 2, index, lat, lon, 80, 80)
    # #reg_inc = bam_analysis.build_cov_arr_1st_dem(cld_occ_anom, 'reg_itc', index, lat, lon, 80, 80)
    # reg_slp = bam_analysis.build_cov_arr_1st_dem(cld_occ_anom, 'reg_slp', index, lat, lon, 80, 80)
    corr_arr, reg_slp, reg_inc = bam_analysis.build_cov_arr_1st_dem(hold_test, 'all', index, lat, lon, -80,-20, 0.06, 'noname')
    
    #rval, pval, std_err = bam_analysis.build_cov_arr_1st_dem(hold_test, 'stats', index, lat, lon, 50, 150, 0.05)

    # mast_plot.oned_all_stats(lat, lon, corr_arr, reg_inc, reg_slp, 'ENSO', 'SON')
    #mast_plot.bam_stats_package(pval, rval, std_err)
    mast_plot.bam_cld_corr(corr_arr, reg_inc, reg_slp, season, lat, lon)
    exit()
    # mast_plot.oned_corr(lat, lon, cov_arr)
    # exit()
    # plt.imshow(cov_arr, origin='lower', aspect = 'equal', extent=[-180,180,-90,90], cmap='RdBu')
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')
    # exit()
    cov_arr = jja_global_avg.test_data(cov_arr)
    std_test = np.std(cov_arr.flatten(), ddof=1)
    std_2 = std_test*2
    stat_arr = spatial_sig(cov_arr, std_2)

    #**Relationship significance testing.
    #l_data, h_data = go_stats.merge_data(hold_test, index)
    #atmos_6040_final.ttest(l_data.flatten(), h_data.flatten())
    #atmos_6040_final.permutation(l_data.flatten(), h_data.flatten())

    sig = np.argwhere(stat_arr == 1)
    nsig = np.argwhere(stat_arr == 0)
    lat_arr, lon_arr = get_lat_lon(cov_arr, sig, lat, lon)

    #**This section of code is for merra cloud occurrence data, specifically looking at the south pole
    # c_lat, c_lon, cld_he = overview.call_main(high_enso_years, months, 0)
    # c_lat, c_lon, cld_le = overview.call_main(low_enso_years, months, 0)
    # c_lat, c_lon, cld_mean = overview.call_main(mean, months, 0)

    # mast_plot.cld_occ_sp(lat, lon, cov_arr, cld_le, cld_he, cld_mean)
    # exit()
    if restrict_domain_call == 1:
        mast_plot.mast_corr(lat, lon, hold_test, cov_arr, index, 'SA', lat_arr, lon_arr, 'JJA', 'JJA', z_top, dz_top)
    else:
        mast_plot.oned_corr(lat, lon, cov_arr, lat_arr, lon_arr, 'SON', 'SON', z_top, dz_top, 'sam')



