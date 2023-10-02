"""
*This code is used to make animations of the raw data to see what the values look like. Trying to
diagnose the issues with EKE that keep popping up 

To do:
    The plots after the mass weighting look right but are just very low. Likely a units thing.

Leave notes:

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
import matplotlib
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
import bam_merra_calc
import bam_merra_calc_d2
import make_bam_index_smal
import go_stats
import mast_plot

years_arr = np.arange(2010,2011)
#*I did this to skip year =2011
years_arr = np.append(years_arr, np.arange(2012,2017))
years = go_stats.make_np_list(years_arr)
#months = ['06','07','08']
#months = ['06']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
days_arr = np.arange(1, 32, 4)
#days_arr = np.arange(1, 4)
days_send = go_stats.make_np_list(days_arr)
days = go_stats.fix_days(days_send)

latlon_bounds = [-20,-70,180,-180]
lat_height_arr = [-70,-20,900,200]
mix_years = 0

bam_smal_data = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bamtest'
bam_data = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bam/'

def save_u_anom(u_anom_send, u3d_send, v_anom_send, v3d_send, eke_send, 
                eke3d_send, t_anom_send, level_send, lat_send, lon_send, date_s):
    #time_send = np.arange(0,u_anom_send.shape[0])

    os.chdir(bam_smal_data)
    print(level_send.shape)
    f = cdf.Dataset('u_anom_test3d_2010_2016.nc4', 'w', format='NETCDF4')
    time = f.createDimension('time')
    level = f.createDimension('level')
    lat = f.createDimension('lat')
    lon = f.createDimension('lon')

    time = f.createVariable('time', 'f4', ('time',))
    level = f.createVariable('level', 'f4', ('level',))
    lat = f.createVariable('lat', 'f4', ('lat',))
    lon = f.createVariable('lon', 'f4', ('lon',))
    u_anom = f.createVariable('u_anom', 'f4', ('time','level','lat'))
    u_anom3d = f.createVariable('u_anom3d', 'f4', ('time','lon','level','lat'))
    v_anom = f.createVariable('v_anom', 'f4', ('time','level','lat'))
    v_anom3d = f.createVariable('v_anom3d', 'f4', ('time','lon','level','lat'))
    eke_data = f.createVariable('eke_data', 'f4', ('time','level','lat'))
    eke_data3d = f.createVariable('eke_data3d', 'f4', ('time','lon','level','lat'))
    t_anom = f.createVariable('t_anom', 'f4', ('time','level','lat'))

    time[:] = date_s
    level[:] = level_send
    lat[:] = lat_send
    lon[:] = lon_send
    u_anom[:] = u_anom_send
    u_anom3d[:] = u3d_send
    v_anom[:] = v_anom_send
    v_anom3d[:] = v3d_send
    eke_data[:] = eke_send
    eke_data3d[:] = eke3d_send
    t_anom[:] = t_anom_send

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

def import_3ddata():
    os.chdir(bam_smal_data)
    f = cdf.Dataset('u_anom_test3d_1986.nc4')

    u3d_get = f.variables['u_anom3d']

    u3d = u3d_get[:]

    print(u3d.shape)
    return u3d


def make_meridional_anom(data_send, mean_send, t_send, levels, lats):
    #!!!!Levels needs to be the flipped levels or the domain restriction won't work.
    xstep = np.arange(0,data_send.shape[2])
    bounds_arr = jja_global_avg.restrict_domain(levels, lats, lat_height_arr[-1], lat_height_arr[2],
                                                lat_height_arr[1],lat_height_arr[0])
    mean_r, levels_r, lat_r = jja_global_avg.region_arr(mean_send, levels, lats, bounds_arr)
    
    for x in xstep:
        #Get each height*lat grid
        t = t_send[:,:,x]
        data = data_send[:,:,x]

        #Flip it (why levels need to flipped)
        t = np.flipud(t)
        data = np.flipud(data)

        #restrict data for height*lat grid
        data_r, levels_r, lat_r = jja_global_avg.region_arr(data, levels, lats, bounds_arr)
        t_r, levels_r, lat_r = jja_global_avg.region_arr(t, levels, lats, bounds_arr)

        anom_s = data_r-mean_r
    
        #apply mass weighted * lat weighted filters
        #print('starting lon: '+ str(x))
        #mast_plot.make_weight_plot(anom, lat_r, levels_r, 1)
        #anom = go_stats.mass_weight(levels_r, anom, t_r)
        #* Only need to make the weights the first time.
        # if x == 0:
        #     mass_weight = go_stats.mass_weight2(levels_r, anom, t_r)
        #     lat_weight = go_stats.latitudinal_weight2(lats_r, anom, 'cos')
        # anom_s = anom*mass_weight*lat_weight
        #anom = go_stats.latitudinal_weight(lat_r, anom, 'cos')
        #mast_plot.make_weight_plot(anom, lat_r, levels, 2)

        if x == 0:
            anom_arr = np.array([anom_s])
            continue
        anom_arr = np.append(anom_arr, np.array([anom_s]), axis=0)
    print(anom_arr.shape)
        
    return anom_arr

def make_eke(u_send, v_send, lat_weight_s):
    xstep = np.arange(0,u_send.shape[0])
    for x in xstep:
        eke_hold = (u_send[x,:,:]**2+v_send[x,:,:]**2)/2
        eke_hold = eke_hold*lat_weight_s
        if x == 0:
            eke_arr = np.array([eke_hold])
            continue
        eke_arr = np.append(eke_arr, np.array([eke_hold]), axis=0)

    return eke_arr

def apply_weights(u_send, v_send, lat_weight_s):
    xstep = np.arange(0,u_send.shape[0])
    for x in xstep:
        u_hold = u_send[x,:,:]*lat_weight_s
        v_hold = v_send[x,:,:]*lat_weight_s
        if x == 0:
            u_arr = np.array([u_hold])
            v_arr = np.array([v_hold])
            continue
        u_arr = np.append(u_arr, np.array([u_hold]), axis=0)
        v_arr = np.append(v_arr, np.array([v_hold]), axis=0)

    return u_arr, v_arr

#* This code is for making the animation.
# if __name__ == '__main__':
#     level, lat, u_anom, u_data, v_anom, eke_data = import_data()
#     lat_m, lev_m, u_m, v_m, t_m = make_bam_index_smal.get_seasonal_mean('jja')
#     tstep = np.arange(u_anom.shape[0])
#     bounds_arr = jja_global_avg.restrict_domain(level, lat, lat_height_arr[-1],lat_height_arr[2],lat_height_arr[1],lat_height_arr[0])
    
#     cmap_send = matplotlib.cm.get_cmap('RdBu')
#     shifted_cmap1 = mast_plot.shiftedColorMap(cmap_send, 1, 0.5, 0)

#     u_m_r, level_r, lat_r = jja_global_avg.region_arr(u_m, lev_m, lat_m, bounds_arr)
#     v_m_r, level_r, lat_r = jja_global_avg.region_arr(v_m, lev_m, lat_m, bounds_arr)
#     for t in tstep:
#         print(t)
#         u_anom_r, level_r, lat_r = jja_global_avg.region_arr(u_anom[t,:,:], level, lat, bounds_arr)
#         u_data_r, level_r, lat_r = jja_global_avg.region_arr(u_data[t,:,:], level, lat, bounds_arr)
#         eke_data_r, level_r, lat_r = jja_global_avg.region_arr(eke_data[t,:,:], level, lat, bounds_arr)
#         v_anom_r, level_r, lat_r = jja_global_avg.region_arr(v_anom[t,:,:], level, lat, bounds_arr)

#         my_eke = (u_anom_r**2+v_anom_r**2)/2

#         if t == 0:
#             u_anom_arr = np.array([u_anom_r])
#             eke_arr = np.array([eke_data_r])
#         elif t!=0:
#             u_anom_arr = np.append(u_anom_arr, np.array([u_anom_r]), axis=0)
#             eke_arr = np.append(eke_arr, np.array([eke_data_r]), axis=0)

#         fig, ax = plt.subplots(4, 1, figsize=(10,10))
#         ax1 = ax[0]
#         norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=u_anom_r.min(), vmax=u_anom_r.max())
#         p = ax1.imshow(u_anom_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap=shifted_cmap1, norm=norm)
#         plt.colorbar(p, ax=ax1)
#         ax1.set_title('Anomalous U. Day= 06/01/1986 + ' + str(t) + ' days.')
#         ax1.set_xlabel('Latitude')
#         ax1.set_ylabel('Pressure (hPa)')

#         # ax2 = ax[1]
#         # p = ax2.imshow(my_eke, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
#         # plt.colorbar(p, ax=ax2)
#         # ax2.set_title('My Calculated EKE')
#         # ax2.set_xlabel('Latitude')
#         # ax2.set_ylabel('Pressure (hPa)')

#         # ax2 = ax[1]
#         # p = ax2.imshow(u_data_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
#         # plt.colorbar(p, ax=ax2)
#         # ax2.set_title('U data. Day= 06/01/1986 + ' + str(t) + ' days.')
#         # ax2.set_xlabel('Latitude')
#         # ax2.set_ylabel('Pressure (hPa)')

#         ax3 = ax[2]
#         vnorm = mcolors.TwoSlopeNorm(vcenter=0,vmin=v_anom.min(),vmax=v_anom.max())
#         p = ax3.imshow(v_anom_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap=shifted_cmap1)
#         plt.colorbar(p, ax=ax3)
#         ax3.set_title('Anomalous V. Day= 06/01/1986 + ' + str(t) + ' days.')
#         ax3.set_xlabel('Latitude')
#         ax3.set_ylabel('Pressure (hPa)')

#         # ax3 = ax[2]
#         # p = ax3.imshow(u_m_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
#         # plt.colorbar(p, ax=ax3)
#         # ax3.set_title('JJA Seasonal Mean U')
#         # ax3.set_xlabel('Latitude')
#         # ax3.set_ylabel('Pressure (hPa)')

#         ax4 = ax[1]
#         mean_anom = np.mean(eke_arr, axis=0)
#         #mnorm = mcolors.TwoSlopeNorm(vcenter=0,vmin=mean_anom.min(),vmax=mean_anom.max())
#         p = ax4.imshow(mean_anom, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
#         plt.colorbar(p, ax=ax4)
#         ax4.set_title('Rolling mean of EKE Data')
#         ax4.set_xlabel('Latitude')
#         ax4.set_ylabel('Pressure (hPa)')

#         ax5 = ax[3]
#         mean_anom = np.mean(u_anom_arr, axis=0)
#         munorm = mcolors.TwoSlopeNorm(vcenter=0,vmin=mean_anom.min(),vmax=mean_anom.max())
#         p = ax5.imshow(mean_anom, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
#         plt.colorbar(p, ax=ax5)
#         ax5.set_title('Rolling mean of anomalous U data')
#         ax5.set_xlabel('Latitude')
#         ax5.set_ylabel('Pressure (hPa)')

#         # ax4 = ax[3]
#         # p = ax4.imshow(eke_data_r, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
#         # plt.colorbar(p, ax=ax4)
#         # ax4.set_title('EKE Values')
#         # ax4.set_xlabel('Latitude')
#         # ax4.set_ylabel('Pressure (hPa)')

#         plt.tight_layout()
#         plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/u_anim_wmean2/'+str(t)+'.png')
#         plt.close()


#* This code was for making the dataset.

if __name__ == '__main__':
    for year in years:
        for month in months:
            if mix_years == 1 and month == '01':
                year = int(year)+1
                year = str(year)
            season = make_bam_index_smal.season_select(int(month))
            for day in days:
                print('starting: '+ day+'-'+month+'-'+year)
                lat, lon, lev, time, u, v, t = bam_merra_calc_smal.call_main(month, year, day)
                #lat, lon, lev, time, u, v, t, fname = bam_merra_calc_smal.import_data()
                if lat[0] == 'n':
                    print('Ommitting: ' + month+ '/'+year+'/'+day)
                    continue
                print(season)
                lat_m, lev_m, u_m, v_m, t_m = make_bam_index_smal.get_seasonal_mean(season)
                d_mean_u = np.mean(u, axis=0)
                d_mean_v = np.mean(v, axis=0)
                d_mean_t = np.mean(t, axis=0)

                bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0], 
                                                              latlon_bounds[1], latlon_bounds[2], latlon_bounds[3])
                d_mean_u_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_u, lat, lon, bounds_arr)
                d_mean_v_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_v, lat, lon, bounds_arr)
                d_mean_t_r, lats_r, lons_r = jja_global_avg.region_arr3d(d_mean_t, lat, lon, bounds_arr) 

                zonal_mean_u = np.mean(d_mean_u_r, axis=2)
                zonal_mean_v = np.mean(d_mean_v_r, axis=2)
                zonal_mean_t = np.mean(d_mean_t_r, axis=2) 

                # d_mean_u = np.flipud(d_mean_u_r)
                # d_mean_v = np.flipud(d_mean_v_r)
                # d_mean_t = np.flipud(d_mean_t_r)

                lev_send = np.flipud(lev)   

                #Needed to do this to make the lat weight right. Very bandaid like, should definitely fix.
                if day == days[0] and month == months[0] and year == years[0]:
                    bounds_arr_lat_height = jja_global_avg.restrict_domain(lev_send, lats_r, lat_height_arr[-1], lat_height_arr[2],
                                                                       lat_height_arr[1], lat_height_arr[0])
                    u_rest, lev_r, lats_r = jja_global_avg.region_arr(d_mean_u_r[:,:,10], lev_send, lats_r, bounds_arr_lat_height)
                    lat_weight = go_stats.latitudinal_weight2(lats_r, u_rest, 'cos')

            #* Need to send temperature & levels to do the mass weighting. 
            # Both the mass weight and lat weight all are in this function.
            # !!!Levels needs to be the flipped levels or the domain restriction won't work.
            # I flip the u&t values in the function.
                u_anom3d = make_meridional_anom(d_mean_u_r, u_m, d_mean_t_r, lev_send, lats_r)
                v_anom3d = make_meridional_anom(d_mean_v_r, v_m, d_mean_t_r, lev_send, lats_r)

                u_anom = zonal_mean_u - u_m
                v_anom = zonal_mean_v - v_m
                t_anom = zonal_mean_t - t_m

                u_anom_r, lev_r, lat_r = jja_global_avg.region_arr(u_anom, lev_send, lats_r, bounds_arr_lat_height) 
                v_anom_r, lev_r, lat_r = jja_global_avg.region_arr(v_anom, lev_send, lats_r, bounds_arr_lat_height)
                t_anom_r, lev_r, lats_r = jja_global_avg.region_arr(t_anom, lev_send, lat_r, bounds_arr_lat_height)

                u_anom_r = u_anom_r*lat_weight
                v_anom_r = v_anom_r*lat_weight
                t_anom_r = t_anom_r*lat_weight

                #eke_hold = ((d_mean_u-u_m)**2+(d_mean_v-v_m)**2)/2
                eke_3d = make_eke(u_anom3d, v_anom3d, lat_weight)
                #*Need to weight the u&v data because the weight is applied in the make_eke function
                u_anom3d_w, v_anom3d_w = apply_weights(u_anom3d, v_anom3d, lat_weight)
                print(eke_3d.shape)
                eke_hold = np.mean(eke_3d, axis=0)
                if day == days[0]:
                    mast_plot.show_eke(eke_3d[10,:,:], eke_hold, lev_send, lons_r[10])

                if day == days[0] and month == months[0] and year==years[0]:
                    u_anom_arr = np.array([u_anom_r])
                    u_anom3d_arr = np.array([u_anom3d_w])
                    v_anom_arr = np.array([v_anom_r])
                    v_anom3d_arr = np.array([v_anom3d_w])
                    eke_arr = np.array([eke_hold])
                    eke_arr3d = np.array([eke_3d])
                    t_anom_arr = np.array([t_anom_r])
                    date_arr = np.array([str(month+day+year)])
                    continue

                u_anom_arr = np.append(u_anom_arr, np.array([u_anom_r]), axis=0)
                u_anom3d_arr = np.append(u_anom3d_arr, np.array([u_anom3d]), axis=0)
                v_anom_arr = np.append(v_anom_arr, np.array([v_anom_r]), axis=0)
                v_anom3d_arr = np.append(v_anom3d_arr, np.array([v_anom3d]), axis=0)
                eke_arr = np.append(eke_arr, np.array([eke_hold]), axis=0)
                eke_arr3d = np.append(eke_arr3d, np.array([eke_3d]), axis=0)
                t_anom_arr = np.append(t_anom_arr, np.array([t_anom_r]), axis=0)
                date_arr = np.append(date_arr, str(month+day+year))


                print(u_anom_arr.shape)
                print(u_anom3d_arr.shape)
                print(v_anom_arr.shape)
                print(v_anom3d_arr.shape)
                print(eke_arr.shape)
                print(eke_arr3d.shape)

    save_u_anom(u_anom_arr, u_anom3d_arr, v_anom_arr, v_anom3d_arr, eke_arr,
                 eke_arr3d, t_anom_arr, lev_r, lats_r, lons_r, date_arr)
