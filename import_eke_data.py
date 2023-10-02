"""
*This code is used to make animations of the raw data to see what the values look like. Trying to
diagnose the issues with EKE that keep popping up - This code was a copy of make_raw_bam_anim.py, 
because I was running some code on that. I just talked to Jay about a cosine latitudinal weighting 
which is what I'm going to try to implement.

To do:

Leave notes:
#* Don't run right now. Really the only thing that's that useful is the 3d_import data call.


"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
import matplotlib
import datetime as dt
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
import indexes

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

def save_u_anom(u_anom_send, u_d_send, v_anom_send, eke_send, level_send, lat_send):
    time_send = np.arange(0,u_anom_send.shape[0])

    os.chdir(bam_smal_data)
    f = cdf.Dataset('u_anom_test_1986.nc4', 'w', format='NETCDF4')
    time = f.createDimension('time')
    level = f.createDimension('level')
    lat = f.createDimension('lat')

    time = f.createVariable('time', 'f4', ('time',))
    level = f.createVariable('level', 'f4', ('level',))
    lat = f.createVariable('lat', 'f4', ('lat',))
    u_anom = f.createVariable('u_anom', 'f4', ('time','level','lat'))
    u_data = f.createVariable('u_data', 'f4', ('time','level','lat'))
    v_anom = f.createVariable('v_anom', 'f4', ('time','level','lat'))
    eke_data = f.createVariable('eke_data', 'f4', ('time','level','lat'))

    time[:] = time_send
    level[:] = level_send
    lat[:] = lat_send
    u_anom[:] = u_anom_send
    u_data[:] = u_d_send
    v_anom[:] = v_anom_send
    eke_data[:] = eke_send

    f.close()

def import_data():
    os.chdir(bam_smal_data)
    f = cdf.Dataset('u_anom_test.nc4')

    level_get = f.variables['level']
    lat_get = f.variables['lat']
    u_anom_get = f.variables['u_anom']
    u_data_get = f.variables['u_data']
    v_anom_get = f.variables['v_anom']
    eke_data_get = f.variables['eke_data']

    level = level_get[:]
    lat = lat_get[:]
    u_anom = u_anom_get[:]
    u_data = u_data_get[:]
    v_anom = v_anom_get[:]
    eke_data = eke_data_get[:]
    
    return level, lat, u_anom, u_data, v_anom, eke_data

def import_data3d(filename, time_flag):
    """
    Imports the 4d data from make_raw_bam_anim.py. Has the eke_data that we use to create the 
    bam EOF. IN this initial stage, there's one set of the data that has a time label dimension, 
    the time_flag == 1 if time dimension is desired. 
    Returns level, lat, lon, time (if time_flag ==1), u_anom, u3d, v_anom, v3d, eke_data, eke_3d
    """
    os.chdir(bam_smal_data)
    #f = cdf.Dataset('u_anom_test3d_2006_2009.nc4')
    f = cdf.Dataset(filename)

    level_get = f.variables['level']
    lat_get = f.variables['lat']
    lon_get = f.variables['lon']
    if time_flag == 1:
        time_get = f.variables['time']
        time = time_get[:]
    u_anom_get = f.variables['u_anom']
    u3d_get = f.variables['u_anom3d']
    v_anom_get = f.variables['v_anom']
    v3d_get = f.variables['v_anom3d']
    eke_data_get = f.variables['eke_data']
    eke3d_get = f.variables['eke_data3d']

    level = level_get[:]
    lat = lat_get[:]
    lon = lon_get[:]
    u_anom = u_anom_get[:]
    u3d = u3d_get[:]
    v_anom = v_anom_get[:]
    v3d = v3d_get[:]
    eke_data = eke_data_get[:]
    eke3d = eke3d_get[:]
    
    if time_flag ==1:
        return level, lat, lon, time, u_anom, u3d, v_anom, v3d, eke_data, eke3d

    else:
        return level, lat, lon, u_anom, u3d, v_anom, v3d, eke_data, eke3d

def make_dt(time_s):
    """
    This function will make time values from the 2nd dataset into datetime objects. 
    """
    tstep = np.arange(0,time_s.shape[0])
    for t in tstep:
        t_hold = dt.datetime.strptime(str(time_s[t]), '%m%d%Y.0')
        if t == 0:
            t_arr = np.array([t_hold])
            continue
        t_arr = np.append(t_arr, t_hold)

    return t_arr

def fix_date_array(dt_s):
    """
    This function fixes the array of times that was saved with the EKE data. The issue is that the 
    months do not have a leading zero and no seperation between month, day, and year values. Adding 
    a leading zero should solve this problem. 
    """
    tstep = np.arange(0,dt_s.shape[0])
    for t in tstep:
        if len(str(dt_s[t])) == 9:
            print('fixing')
            dt_hold = '0'+str(dt_s[t])
        else:
            dt_hold = dt_s[t]
        if t == 0:
            dt_arr = np.array([dt_hold])
            continue
        dt_arr = np.append(dt_arr, dt_hold)
    print(dt_arr.shape)

    return dt_arr

def merge_data_3d():
    """
    This function is required because eke_data & anomalous u,v,&t values are likely going to 
    come in in multiple netcdfs to improve the efficiency of data generation. This code will merge
    the values together into one array.

    """
    os.chdir(bam_smal_data)
    for filename in os.listdir(bam_smal_data):
        if filename.split('_')[0] == 'u' and filename.split('_')[3] == '2006':
            level1, lat1, lon1, u_anom1, u3d1, v_anom1, v3d1, eke_data1, eke3d1 = import_data3d(filename, 0)
        if filename.split('_')[0] == 'u' and filename.split('_')[3] == '2010':
            level2, lat2, lon2, time2, u_anom2, u3d2, v_anom2, v3d2, eke_data2, eke3d2 = import_data3d(filename, 1)
    test = fix_date_array(time2)
    time2 = make_dt(test)
    d_arr_send = np.arange(1,32,4)
    time1 = indexes.make_dt_arr(d_arr_send, 4, 12, '01-01-2006', 2011)
    time = np.append(time1, time2)
    u_anom = np.append(u_anom1, u_anom2, axis=0)
    u3d = np.append(u3d1, u3d2, axis=0)
    v_anom = np.append(v_anom1, v_anom2, axis=0)
    v3d = np.append(v3d1, v3d2, axis=0)
    eke_data = np.append(eke_data1, eke_data2, axis=0)
    eke3d = np.append(eke3d1, eke3d2, axis=0)

    return level1, lat1, lon1, time, u_anom, u3d, v_anom, v3d, eke_data, eke3d
    
if __name__ == '__main__':
    merge_data_3d()
    exit()

#* This code is for making the animation.
if __name__ == '__main__':
    level, lat, u_anom, u_data, v_anom, eke_data = import_data()
    level, lat, u_anom, u3d, v_anom, v3d, eke_data, eke3d = import_data3d()
    lat_m, lev_m, u_m, v_m, t_m = make_bam_index_smal.get_seasonal_mean('jja')
    tstep = np.arange(u_anom.shape[0])
    bounds_arr = jja_global_avg.restrict_domain(level, lat, lat_height_arr[-1],lat_height_arr[2],lat_height_arr[1],lat_height_arr[0])

    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = mast_plot.shiftedColorMap(cmap_send, 1, 0.5, 0)

    u_m_r, level_r, lat_r = jja_global_avg.region_arr(u_m, lev_m, lat_m, bounds_arr)
    v_m_r, level_r, lat_r = jja_global_avg.region_arr(v_m, lev_m, lat_m, bounds_arr)
    for t in tstep:
        print(t)
        u_anom_r, level_r, lat_r = jja_global_avg.region_arr(u_anom[t,:,:], level, lat, bounds_arr)
        u_data_r, level_r, lat_r = jja_global_avg.region_arr(u_data[t,:,:], level, lat, bounds_arr)
        eke_data_r, level_r, lat_r = jja_global_avg.region_arr(eke_data[t,:,:], level, lat, bounds_arr)
        v_anom_r, level_r, lat_r = jja_global_avg.region_arr(v_anom[t,:,:], level, lat, bounds_arr)

        u_anom_f = go_stats.latitudinal_weight(lat_r, u_anom_r, 'cos')
        u_data_f = go_stats.latitudinal_weight(lat_r, u_data_r, 'cos')
        eke_data_f = go_stats.latitudinal_weight(lat_r, eke_data_r, 'cos')
        v_anom_f = go_stats.latitudinal_weight(lat_r, v_anom_r, 'cos')

        if t == 37:
            mast_plot.make_weight_plot(u_anom_r, lat_r, level_r, t)
            exit()
        else:
            continue

        my_eke = (u_anom_r**2+v_anom_r**2)/2


        if t == 0:
            u_anom_arr = np.array([u_anom_f])
            eke_arr = np.array([eke_data_])
        elif t!=0:
            u_anom_arr = np.append(u_anom_arr, np.array([u_anom_f]), axis=0)
            eke_arr = np.append(eke_arr, np.array([eke_data_f]), axis=0)

        fig, ax = plt.subplots(4, 1, figsize=(10,10))
        ax1 = ax[0]
        norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=u_anom_f.min(), vmax=u_anom_f.max())
        p = ax1.imshow(u_anom_f, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap=shifted_cmap1, norm=norm)
        plt.colorbar(p, ax=ax1)
        ax1.set_title('Anomalous U. Day= 06/01/1986 + ' + str(t) + ' days.')
        ax1.set_xlabel('Latitude')
        ax1.set_ylabel('Pressure (hPa)')

        # ax2 = ax[1]
        # p = ax2.imshow(my_eke, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
        # plt.colorbar(p, ax=ax2)
        # ax2.set_title('My Calculated EKE')
        # ax2.set_xlabel('Latitude')
        # ax2.set_ylabel('Pressure (hPa)')

        # ax2 = ax[1]
        # p = ax2.imshow(u_data_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
        # plt.colorbar(p, ax=ax2)
        # ax2.set_title('U data. Day= 06/01/1986 + ' + str(t) + ' days.')
        # ax2.set_xlabel('Latitude')
        # ax2.set_ylabel('Pressure (hPa)')

        ax3 = ax[1]
        vnorm = mast_plot.norm_maker(v_anom_f)
        p = ax3.imshow(v_anom_f, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap=shifted_cmap1)
        plt.colorbar(p, ax=ax3)
        ax3.set_title('Anomalous V. Day= 06/01/1986 + ' + str(t) + ' days.')
        ax3.set_xlabel('Latitude')
        ax3.set_ylabel('Pressure (hPa)')

        # ax3 = ax[2]
        # p = ax3.imshow(u_m_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
        # plt.colorbar(p, ax=ax3)
        # ax3.set_title('JJA Seasonal Mean U')
        # ax3.set_xlabel('Latitude')
        # ax3.set_ylabel('Pressure (hPa)')

        ax4 = ax[2]
        mean_anom = np.mean(eke_arr, axis=0)
        #mnorm = mcolors.TwoSlopeNorm(vcenter=0,vmin=mean_anom.min(),vmax=mean_anom.max())
        p = ax4.imshow(mean_anom, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
        plt.colorbar(p, ax=ax4)
        ax4.set_title('Rolling mean of EKE Data')
        ax4.set_xlabel('Latitude')
        ax4.set_ylabel('Pressure (hPa)')

        ax5 = ax[3]
        mean_anom = np.mean(u_anom_arr, axis=0)
        munorm = mcolors.TwoSlopeNorm(vcenter=0,vmin=mean_anom.min(),vmax=mean_anom.max())
        p = ax5.imshow(mean_anom, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap=shifted_cmap1, norm=munorm)
        plt.colorbar(p, ax=ax5)
        ax5.set_title('Rolling mean of anomalous U data')
        ax5.set_xlabel('Latitude')
        ax5.set_ylabel('Pressure (hPa)')

        # ax4 = ax[3]
        # p = ax4.imshow(eke_data_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
        # plt.colorbar(p, ax=ax4)
        # ax4.set_title('EKE Values')
        # ax4.set_xlabel('Latitude')
        # ax4.set_ylabel('Pressure (hPa)')

        plt.tight_layout()
        plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/u_anim_wmean2/'+str(t)+'.png')
        plt.close()


#* This code was for making the dataset.
# if __name__ == '__main__':
#     for year in years:
#         for month in months:
#             if mix_years == 1 and month == '01':
#                 year = int(year)+1
#                 year = str(year)
#             season = make_bam_index_smal.season_select(int(month))
#             for day in days:
#                 print('starting: '+ day+'-'+month+'-'+year)
#                 lat, lon, lev, time, u, v, t = bam_merra_calc_smal.call_main(month, year, day)
#                 if lat[0] == 'n':
#                     print('Ommitting: ' + month+ '/'+year+'/'+day)
#                     continue
#                 lat_m, lev_m, u_m, v_m, t_m = make_bam_index_smal.get_seasonal_mean(season)
#                 d_mean_u = np.mean(u, axis=0)
#                 d_mean_v = np.mean(v, axis=0)
#                 d_mean_t = np.mean(t, axis=0)

#                 bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0], 
#                                                               latlon_bounds[1], latlon_bounds[2], latlon_bounds[3])
#                 d_mean_u_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_u, lat, lon, bounds_arr)
#                 d_mean_v_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_v, lat, lon, bounds_arr)
#                 d_mean_t_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_t, lat, lon, bounds_arr) 

#                 f_mean_u = np.mean(d_mean_u_r, axis=2)
#                 f_mean_v = np.mean(d_mean_v_r, axis=2)
#                 f_mean_t = np.mean(d_mean_t_r, axis=2) 

#                 f_mean_u = np.flipud(f_mean_u)
#                 f_mean_v = np.flipud(f_mean_v)
#                 f_mean_t = np.flipud(f_mean_t)

#                 lev_send = np.flipud(lev)      

#                 u_anom = f_mean_u - u_m
#                 v_anom = f_mean_v - v_m
#                 eke_hold = ((f_mean_u-u_m)**2+(f_mean_v-v_m)**2)/2
#                 if day == days[0] and month == months[0] and year==years[0]:
#                     u_anom_arr = np.array([u_anom])
#                     u_arr = np.array([f_mean_u])
#                     v_anom_arr = np.array([v_anom])
#                     eke_arr = np.array([eke_hold])
#                     continue
#                 u_anom_arr = np.append(u_anom_arr, np.array([u_anom]), axis=0)
#                 u_arr = np.append(u_arr, np.array([f_mean_u]), axis=0)
#                 v_anom_arr = np.append(v_anom_arr, np.array([v_anom]), axis=0)
#                 eke_arr = np.append(eke_arr, np.array([eke_hold]), axis=0)

#     save_u_anom(u_anom_arr, u_arr, v_anom_arr, eke_arr, lev_send, lats_r)
