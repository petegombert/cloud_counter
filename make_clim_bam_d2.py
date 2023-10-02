"""
* This script makes the climatological mean ncfile. It will be a 1d file with doy as the 
dimension. This will make generating the monthly BAM numbers less resource intensive. This is the same exact 
script as 'make_clim_bam.py', it just calls 'bam_merra_calc_d2.py' instead of 'bam_merra_calc.py'. The reason
for this is that the merra_hold folder can only have one file in it, so we have multiple folders so I can run
multiple sets of the climatological mean code at 1 time. Furthermore, other variations of 'make_clim_bam.py' 
will have the similar d<#> added to the end of file. This will allow me to submit multiple files as jobs.

#!! Things that need to be changed between runs: 1. Time period
                                                 2. netcdf name
                                                 3. Check year select in loop (for making array)

Leave notes:

To do:
    Could possibly speed up the program by clearing the month and day mast arrays once the mean is 
    made. I think they might not be clearing. - Update they seem to be clearing. Still could probably speed
    it up.


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
import bam_merra_calc_d2
import go_stats

years_arr = np.arange(2003,2004)
years = go_stats.make_np_list(years_arr)
#months = ['01','02','03','04','05','06','07','08','09','10','11','12']
months = ['09','10','11']
days_send = go_stats.make_np_list(np.arange(1,3))
days = go_stats.fix_days(days_send)
# years = ['1980','1981']
# months = ['01','02','03']
mix_years = 0

latlon_bounds = [-20,-70,180,-180]

def get_mean(data):
    tot_u = data[:,4]
    tot_v = data[:,5]

    mean_u = np.mean(tot_u)
    mean_v = np.mean(tot_v)

    return mean_u, mean_v

def save_nc(mast_u, mast_v, mast_t, lat_send, lev):
    os.chdir('/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bam')
    f = cdf.Dataset('son_climatological_mean_2003_2013.nc4', 'w', format='NETCDF4')
    level = f.createDimension('level')
    lat = f.createDimension('lat')

    level = f.createVariable('level', 'f4', ('level',))
    lat = f.createVariable('lat', 'f4', ('lat',))
    u_tot = f.createVariable('u_tot', 'f4', ('level','lat',))
    v_tot = f.createVariable('v_tot', 'f4', ('level','lat',))
    t_tot = f.createVariable('t_tot', 'f4', ('level','lat',))

    level[:] = lev
    lat[:] = lat_send
    u_tot[:] = mast_u
    v_tot[:] = mast_v
    t_tot[:] = mast_t
    f.close()

if __name__ == '__main__':
#* This is my first guess on daily analysis into climatlogical mean, going to save as a seasonal average.
    for year in years:
        for month in months:
            if month == '01' and mix_years == 1:
                year_hold = int(year)+1
                year = str(year_hold)
                print(year)
            for day in days:
                print('starting: '+str(month)+'-'+day+'-'+year)
                lat_hold, lon_hold, lev_hold, time_hold, tot_u_hold, tot_v_hold, t_hold = bam_merra_calc_d2.call_main(month, year, day)
                if lat_hold[0] == 'n':
                    print('Ommitting: ' + month+ '/'+year+'/'+day)
                    continue

                d_mean_u = np.mean(tot_u_hold, axis=0)
                d_mean_v = np.mean(tot_v_hold, axis=0)
                d_mean_t = np.mean(t_hold, axis=0)

                bounds_arr = jja_global_avg.restrict_domain(lat_hold, lon_hold, latlon_bounds[0], 
                                                              latlon_bounds[1], latlon_bounds[2], latlon_bounds[3])
                d_mean_u_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_u, lat_hold, lon_hold, bounds_arr)
                d_mean_v_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_v, lat_hold, lon_hold, bounds_arr)
                d_mean_t_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_t, lat_hold, lon_hold, bounds_arr)

                f_mean_u_r = np.mean(d_mean_u_r, axis=2)
                f_mean_v_r = np.mean(d_mean_v_r, axis=2)
                f_mean_t_r = np.mean(d_mean_t_r, axis=2)

                f_mean_u_r = np.flipud(f_mean_u_r)
                f_mean_v_r = np.flipud(f_mean_v_r)
                f_mean_t_r = np.flipud(f_mean_t_r)

                lev_r = np.flip(lev_hold)

                #Make array of daily zonal mean values
                if day == days[0]:
                    d_mast_arr = np.array([[day, lats_r, lev_r, f_mean_u_r, f_mean_v_r, f_mean_t_r]], dtype=object)
                    continue
                d_mast_arr = np.append(d_mast_arr, np.array([[day, lats_r, lev_r, f_mean_u_r, f_mean_v_r, f_mean_t_r]], dtype=object), axis=0)
                print(d_mast_arr.shape)

            #*Take mean of daily zonal mean array for each value
            print('Making monthly mean: ' + month)
            m_z_mean_u = np.mean(d_mast_arr[:,3])
            m_z_mean_v = np.mean(d_mast_arr[:,4])
            m_z_mean_t = np.mean(d_mast_arr[:,5])

            #*Make array of monthly zonal mean values
            if month == months[0]:
                m_mast_arr = np.array([[month, d_mast_arr[0,1], d_mast_arr[0,2], m_z_mean_u, m_z_mean_v, m_z_mean_t]], dtype=object)
                continue
            m_mast_arr = np.append(m_mast_arr, np.array([[month, d_mast_arr[0,1], d_mast_arr[0,2], m_z_mean_u, m_z_mean_v, m_z_mean_t]], dtype=object), axis=0)
            print(m_mast_arr.shape)

        #*Take mean of monthly zonal mean values -> makes seasonal mean
        print('Making seasonal mean: ' + year)
        y_z_mean_u = np.mean(m_mast_arr[:,3])
        y_z_mean_v = np.mean(m_mast_arr[:,4])
        y_z_mean_t = np.mean(m_mast_arr[:,5])

        #*Make array of yearly values -> the mean of this will be the decade mean for the season.
        #*This needs to be year == years[1] if we want to mix years
        if year == years[0]:
            y_mast_arr = np.array([[year, m_mast_arr[0,1], m_mast_arr[0,2], y_z_mean_u, y_z_mean_v, y_z_mean_t]], dtype=object)
            continue
        y_mast_arr = np.append(y_mast_arr, np.array([[year, m_mast_arr[0,1], m_mast_arr[0,2], y_z_mean_u, y_z_mean_v, y_z_mean_t]], dtype=object), axis=0)
        print(y_mast_arr.shape)

    #*Take mean of yearly values
    z_mean_u = np.mean(y_mast_arr[:,3])
    z_mean_v = np.mean(y_mast_arr[:,4])
    z_mean_t = np.mean(y_mast_arr[:,5])

    save_nc(z_mean_u, z_mean_v, z_mean_t, y_mast_arr[0,1], y_mast_arr[0,2])

    



    #*This was for my intitial guess into monthly climalogical mean. Want to do daily.
# going to step through years to make the monthly climatological data.
    # for month in months:
    #     for year in years:
    #         lat_hold, lon_hold, lev_hold, tot_u_hold, tot_v_hold = bam_merra_calc.call_main(month, year)
    #         if lat_hold[0] == 'n':
    #             print('Ommitting: ' + month +'/'+year)
    #             continue 
    #         if year == years[0]:
    #             mast_arr = np.array([[year, lat_hold, lon_hold, lev_hold, tot_u_hold, tot_v_hold]], dtype=object)
    #             continue
    #         mast_arr = np.append(mast_arr, np.array([[year, lat_hold, lon_hold, lev_hold, tot_u_hold, tot_v_hold]], dtype=object), axis =0)

    #     mean_u, mean_v = get_mean(mast_arr)
    #     if month == months[0]:
    #         mast_u = mean_u
    #         mast_v = mean_v
    #         continue
    #     mast_u = np.append(mast_u, mean_u, axis=0)
    #     mast_v = np.append(mast_v, mean_v, axis=0)

    # #* This is just to send to the ncfile.
    # month_arr = np.arange(0,13)
    # save_nc(mast_u, mast_v, lat_hold, lon_hold, lev_hold, month_arr)





