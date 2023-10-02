"""
*This code is run after 'make_clim_bam'. I am using it to take the output from that script 
(the climatological-seasonal mean data) and comparing the daily data to it. From this we should 
be able to build the BAM index which is what this code will do. 

*Between runs, do the following:
    1. Change years/season
    2. Change mixyears

        * I don't think filename needs to be changed.

To do:
    I'm a little worried about the dimensions of the climological mean data. It's being saved
    as (lev,lat) which is the opposite that it should be. It didn't error out when saving so maybe it's 
    okay. Definitely should look into it though.

Leave notes:

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
import jja_global_avg
from pydap.client import open_url
from pydap.cas.urs import setup_session
import requests
import bam_merra_calc
import go_stats

years_arr = np.arange(1986,1991)
years = go_stats.make_np_list(years_arr)
#months = ['12','01','02']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
days_arr = np.arange(1, 32)
days_send = go_stats.make_np_list(days_arr)
days = go_stats.fix_days(days_send)

latlon_bounds = [-20,-70,180,-180]
mix_years = 1

bam_data = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bam'

def get_seasonal_mean_yrs(season, years):
    os.chdir(bam_data)
    fname = season+'_climatological_mean_'+years+'.nc4'
    f = cdf.Dataset(fname)
    
    lat_get = f.variables['lat']
    lev_get = f.variables['level']
    u_get = f.variables['u_tot']
    v_get = f.variables['v_tot']
    t_get = f.variables['t_tot']

    lat = lat_get[:]
    lev = lev_get[:]
    u = u_get[:]
    v = v_get[:]
    t = t_get[:]

    return lat, lev, u, v, t

def get_seasonal_mean(season):
    os.chdir(bam_data)
    fname = season + '_full_climatological_mean.nc4'

    f = cdf.Dataset(fname)
    
    lat_get = f.variables['lat']
    lev_get = f.variables['level']
    u_get = f.variables['u_tot']
    v_get = f.variables['v_tot']
    t_get = f.variables['t_tot']

    lat = lat_get[:]
    lev = lev_get[:]
    u = u_get[:]
    v = v_get[:]
    t = t_get[:]

    return lat, lev, u, v, t

def season_select(month):
    #** Month must be in int form!!
    if month < 3 or month == 12:
        return 'djf'
    elif month > 2 and month < 6:
        return 'mam'
    elif month > 5 and month < 9:
        return 'jja'
    else:
        return 'son'
    
def save_eke(data_send, vt_send, u_send, uv_send, levels, lats):
    time_send = np.arange(0,data_send.shape[0])

    os.chdir(bam_data)
    f = cdf.Dataset('eke_metrics_1986_1990.nc4', 'w', format='NETCDF4')
    time = f.createDimension('time')
    level = f.createDimension('level')
    lat = f.createDimension('lat')

    time = f.createVariable('time', 'f4', ('time',))
    level = f.createVariable('level', 'f4', ('level',))
    lat = f.createVariable('lat', 'f4', ('lat',))
    eke_data = f.createVariable('eke_data', 'f4', ('time', 'level', 'lat',))
    vt_data = f.createVariable('vt_data', 'f4', ('time', 'level', 'lat',))
    u_data = f.createVariable('u_data', 'f4', ('time', 'level', 'lat'))
    uv_data = f.createVariable('uv_data', 'f4', ('time', 'level', 'lat'))

    time[:] = time_send
    level[:] = levels
    lat[:] = lats
    eke_data[:] = data_send
    vt_data[:] = vt_send
    u_data[:] = u_send
    uv_data[:] = uv_send

    f.close()

#First, make array of eddy mean kinetic energy (EKE): ((u'-u^2)+(v'-v^2))/2
if __name__ == '__main__':
    for year in years:
        for month in months:
            if mix_years == 1 and month == '01':
                year = int(year)+1
                year = str(year)
            season = season_select(int(month))
            for day in days:
                print('starting: '+ day+'-'+month+'-'+year)
                lat, lon, lev, time, u, v, t = bam_merra_calc.call_main(month, year, day)

                if lat[0] == 'n':
                    print('Ommitting: ' + month+ '/'+year+'/'+day)
                    continue
                lat_m, lev_m, u_m, v_m, t_m = get_seasonal_mean(season)

                #*Make daily mean:
                d_mean_u = np.mean(u, axis=0)
                d_mean_v = np.mean(v, axis=0)
                d_mean_t = np.mean(t, axis=0)

                bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0], 
                                                              latlon_bounds[1], latlon_bounds[2], latlon_bounds[3])
                d_mean_u_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_u, lat, lon, bounds_arr)
                d_mean_v_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_v, lat, lon, bounds_arr)
                d_mean_t_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_t, lat, lon, bounds_arr) 

                f_mean_u = np.mean(d_mean_u_r, axis=2)
                f_mean_v = np.mean(d_mean_v_r, axis=2)
                f_mean_t = np.mean(d_mean_t_r, axis=2)        

                f_mean_u = np.flipud(f_mean_u)
                f_mean_v = np.flipud(f_mean_v)
                f_mean_t = np.flipud(f_mean_t)

                lev_send = np.flipud(lev)      

                eke_hold = ((f_mean_u-u_m)**2+(f_mean_v-v_m)**2)/2
                vt_hold = ((f_mean_v-v_m)*(f_mean_t-t_m))
                u_hold = f_mean_u - u_m
                vu_hold = (f_mean_u-u_m)*(f_mean_v-v_m)

                print(eke_hold.shape)

                if day == days[0] and month == months[0] and year == years[0]:
                    eke_arr = np.array([eke_hold])
                    vt_arr = np.array([vt_hold])
                    u_arr = np.array([u_hold])
                    vu_arr = np.array([vu_hold])
                    continue
                eke_arr = np.append(eke_arr, np.array([eke_hold]), axis=0) #*Need to test, maybe could just do one bracket on array. Not sure what the correct axis is.
                vt_arr = np.append(vt_arr, np.array([vt_hold]), axis=0)
                u_arr = np.append(u_arr, np.array([u_hold]), axis=0)
                vu_arr = np.append(vu_arr, np.array([vu_hold]), axis=0)

    save_eke(eke_arr, vt_arr, u_arr, vu_arr, lev_send, lats_r)
    