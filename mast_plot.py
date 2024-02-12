#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

*This is the first script of the new climate project. This will grab a year of data and plot cloudyness fraction averages for the year.
leave notes: Not sure what cbz_dz or ctz_dz or cbtz is. Might have to ask Jay or look how he uses in IDL.
             	--It appears that there are 2 data arrays. One that the statistics are generated from and another that show global mean statistics. Not sure where the first is going to
			come into play.
	     Data array is going to be very confusing. It's a 5-6 dimensional array. Will need to really think through this one and create a really nice key for myself.
To do: 
        Should add in the if __name__ = '__main__'
         Should maybe put sum_months into get_data.py. I think this would be a useful function to be able to call from numerous scripts 
         Need to transition all the plotting functions into imshow. Need to figure out how to get the colorbar to work with the axis argument. Maybe? 
          Also am going to need 
        For the mast_corr function, I'm calculating the index from indexes code and passing it from oned_corr. Seems inefficient, 
        should probably get rid of one.
        

"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.colors as mcolors
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import os
import numpy as np
import pandas as pd
from pyhdf.HDF import *
from pyhdf.VS import *
import netCDF4 as cdf
import cartopy.crs as ccrs
import cartopy.feature as feature
import julian
from datetime import datetime
import sys
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
import subprocess
from warnings import filterwarnings
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes

np.set_printoptions(threshold=sys.maxsize)
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

try:
    year = sys.argv[1]
except IndexError:
    year = 2006

daytime_filter = 0
min_dbz = -0.
make_plot = 1

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def norm_maker(data):
    #Or could do if np.max == 0.0
    #print(np.min(data), np.max(data))
    try:
        norm = mcolors.TwoSlopeNorm(vcenter=0,vmin=np.min(data),vmax=np.max(data))
    except ValueError:
        norm = mcolors.TwoSlopeNorm(vcenter=0,vmin=np.min(data),vmax=0.0000001)
    return norm

def sum_months():
    data = get_files_dict
    print(data.shape)
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

def make_plot_main():
    data = data_send
    data = np.sum(data, axis=2)
    dbz = get_files_dict[0,4]
    data = sum_along_dbz(data, dbz, min_dbz)

    z = get_files_dict[0][3]-0.5
    lat = get_files_dict[0][1]

#    fig, ax = plt.subplots(1, 1, figsize = (10, 10))
    plt.imshow(data, origin='lower', aspect='auto')
    plt.colorbar()
#    ax.contourf(lat, z, data)
#    ax.set_ylim([0, 15])
#    ax.set_xlim([-75, 75])
#    ax.set_ylabel('Height (km)')
#    ax.set_xlabel('Latitude')
#    ax.set_title('200801-200812 Mean Cloudiness. No min dBZ')
#    ax=plt.gca() #get the current axes
#    PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
#    plt.colorbar(PCM, ax=ax)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/climate_height_imshow.png')
    plt.clf()

def make_zonal_mean_plot(data_send, x_axis, y_axis, title, month):
#pretty basic plotting on the zonal mean
    fig,ax = plt.subplots(figsize=(10,10))

    p=ax.imshow(data_send, origin='lower', aspect='auto', extent=[-90,90,0,20])
    cb = plt.colorbar(p,shrink=0.5)
    ax.set_title(title)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Height (km)')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/monthly_avg_plots/zonal_mean_'+str(month)+'.png')

def make_zonal_anom_plot(data_send, x_axis, y_axis, title, month, year):
#pretty basic plotting on the zonal mean
    fig,ax = plt.subplots(figsize=(10,10))

    p=ax.imshow(data_send, origin='lower', aspect='auto', extent=[-90,90,0,20])
    cb = plt.colorbar(p,shrink=0.5)
    ax.set_title(title)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Height')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/monthly_annom_plots/monthly_annom_'+year+'_'+month+'.png')

def make_mast_anom_plot(data_avg, data_anom, data_mon, x_axis, y_axis, month, year):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    fig,ax = plt.subplots(2, 2, figsize=(10,10))

    ax1 = ax[0,0]
    p1 = ax1.imshow(data_avg, origin='lower', aspect='auto', extent=[-90,90,0,20])
    cb1 = plt.colorbar(p1,shrink=0.5, ax=ax1)
    ax1.set_title('Averaged value for month: '+str(month))
    ax1.set_xlabel('Latitude')
    ax1.set_ylabel('Height (km)')

    ax2 = ax[1,0]
    p2 = ax2.imshow(data_anom, origin='lower', aspect='auto', extent=[-90,90,0,20], cmap='RdBu')
    cb2 = plt.colorbar(p2,shrink=0.5, ax=ax2)
    ax2.set_title('Anomaly values for month: '+str(month)+' year: '+str(year))
    ax2.set_xlabel('Latitude')
    ax2.set_ylabel('Height (km)')

    ax3 = ax[0,1]
    p3 = ax3.imshow(data_mon, origin='lower', aspect='auto', extent=[-90,90,0,20])
    cb3 = plt.colorbar(p3,shrink=0.5,ax=ax3)
    ax3.set_title('Monthly zonal average, month: '+str(month)+', year: '+str(year))
    ax3.set_xlabel('Latitude')
    ax3.set_ylabel('Height (km)')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/monthly_annom_plots_full/'+str(year)+'/monthly_anom_'+str(year)+'_'+month+'.png')
    plt.close()
def make_seas_avg_plot(data_avg, seas):
    fig, ax = plt.subplots(figsize=(10,10))

    p = ax.imshow(data_avg, origin='lower', aspect='auto', extent=[-90,90,0,21])
    cb = plt.colorbar(p,shrink=0.5)
    ax.set_title('Averaged seasonal cloud occurrence for season: '+str(seas))
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Height')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/seas_avg_'+seas+'.png')

def low_cloud(lat, lon, data):
#    data = np.flip(data, axis=0)
    print(data.shape)
    print(lat.shape)
    print(lon.shape)
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1, 1, 1, projection=proj)
    g1 = ax.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.top_labels = False
    g1.right_labels = False
    ax.coastlines(color='white')

#    levels = np.arange(0,700,50)
    p = ax.contourf(lon[0:-1], lat[0:-1], data, cmap='magma')
#    p = ax.imshow(data, origin='lower', aspect='equal', extent=[-180,180,-90,90], cmap='magma')
    c = plt.colorbar(p, shrink=0.5, ax=ax)

    ax.set_title('Global Cloud Occurrence (06/2007) Cloud Base 0-1000m')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/testies.png')
    plt.close()

def overview_plt(lat, lon, data, dz_bot, dz_top, z_bot, z_top, cwd):
    lat = lat[0]
    lon = lon[0]
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1, projection=proj)
    g1 = ax.gridlines(draw_labels=True, linestyle='dotted')
    g1.top_labels = False
    g1.right_labels = False
    ax.coastlines(color='white')

    #p = ax.contourf(lon[0:-1], lat[0:-1], data, cmap='magma')
    p = ax.imshow(data, origin='lower', aspect='equal', extent=[-180,180,-90,90], cmap='magma')
    c = plt.colorbar(p, shrink=0.5, ax=ax)
    fig.suptitle('Global Cloud Occurrence. Cloud base: '+str(z_bot)+'-'+str(z_top)+' Thickness: '+str(dz_bot)+'-'+str(dz_top))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.savefig(cwd+'/'+ str(dz_bot)+'_'+str(dz_top)
                +'__'+str(z_bot)+'_'+str(z_top)+'.png')
    plt.close()

def overview_plt_rd(lat, lon, data, dz_bot, dz_top, z_bot, z_top, bound_arr, r_nickname):
    lat = lat[0]
    lon = lon[0]

    t_bnd = bound_arr[0][0][0]
    b_bnd = bound_arr[1][0][0]
    r_bnd = bound_arr[2][0][0]
    l_bnd = bound_arr[3][0][0]

    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1, projection=proj)
    g1 = ax.gridlines(draw_labels=True, linestyle='dotted')
    g1.top_labels = False
    g1.right_labels = False
    ax.coastlines(color='white')
    p = ax.contourf(lon[0:-1], lat[0:-1], data, cmap='magma')
    c = plt.colorbar(p, shrink=0.5, ax=ax)
    rect = mpatches.Rectangle((lon[l_bnd], lat[b_bnd]), (lon[r_bnd]-lon[l_bnd]), abs(lat[t_bnd]-lat[b_bnd]), fill=False, color='white')
    ax.add_patch(rect)  

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/overview_plots/'+ str(dz_bot)+'_'+str(dz_top)
                +'__'+str(z_bot)+'_'+str(z_top)+'_id_region_'+ r_nickname+'.png')
    plt.close()

def overview_compare(data1, data2):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    
    ax1 = plt.subplot(1,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels=True, linestyle='dotted')
    g1.top_labels = False
    g1.right_labels = False
    ax1.coastlines()
    p = ax1.imshow((data1 - data2), origin='lower', aspect='equal', extent=[-180,180,-90,90], cmap=shifted_cmap, vmin=-0.15, vmax=0.15)
    c = plt.colorbar(p, ax=ax1)

    fig.suptitle('Mean - Low ENSO')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/overview_plots/comps/mean-low_enso.png')

def low_cloud_wregion(lat, lon, data, bound_arr):
#    data = np.flip(data[0], axis=0)
    lon = lon[0]
    lat = lat[0]

    t_bnd = bound_arr[0]
    b_bnd = bound_arr[1]
    r_bnd = bound_arr[2]
    l_bnd = bound_arr[3]

    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1, 1, 1, projection=proj)
    g1 = ax.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.top_labels = False
    g1.right_labels = False
    ax.coastlines(color='white')

#    levels = np.arange(0,700,50)
    p = ax.contourf(lon[0:-1], lat[0:-1], data, cmap='magma')
#    p = ax.imshow(data, origin='lower', aspect='equal', extent=[-180,180,-90,90], cmap='magma')
    c = plt.colorbar(p, shrink=0.5, ax=ax)
    rect = mpatches.Rectangle((lon[l_bnd], lat[b_bnd]), (lon[r_bnd]-lon[l_bnd]), abs(lat[t_bnd]-lat[b_bnd]), fill=False, color='white')
    ax.add_patch(rect)

    ax.set_title('Global Cloud Occurrence (06/2007) Cloud Base 0-1000m, dz 0-1500')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/atmos_6600/jja_avg_showregion_cb1000dz1500.png')
    plt.close()

def time_series(years, index, mean):
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    ax.scatter(years, index, color='red', label='Index')
#    ax.grid(axis='x', which='minor', linestyle='dotted')
    ax.grid(axis='x', which='major', linestyle='dotted')
    ax.legend()
#    ax.set_ylabel('Index Mean', loc='left')
#    ax.set_xlabel('Year')
    ax2 = ax.twinx()
    ax2.scatter(years, mean, color='blue', label='Mean Cloud Occurrence')
    ax2.legend()
#    ax2.set_ylabel('Average JJA Cloud Occurrence', loc='right')
    plt.show()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/test3.png')

def corr(years, index, mean):
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    ax.scatter(index, mean)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/corr_test.png')

def oned_corr(lat, lon, data, lat_sig, lon_sig, ind_months, cld_months, z, dz, osc):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap =shiftedColorMap(cmap_send, 1, 0.48, 0, 'shifted')
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1, 1, 1, projection=proj)
    g1 = ax.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.top_labels = False
    g1.right_labels = False
    ax.coastlines()

    p = ax.imshow(data, origin='lower', aspect = 'equal', extent=[-180,180,-90,90], cmap=shifted_cmap)
    #p = ax.contourf(lon, lat, data, cmap='RdBu')
    
    c = plt.colorbar(p, shrink=0.5, ax=ax)
    tstep = np.arange(0, lat_sig.shape[0])
    for t in tstep:
        ax.scatter(lon_sig[t]+0.75, lat_sig[t]+1.25, color='black')
    fig.suptitle(str(cld_months)+ ' Cloud Occurrence, base z: <'+ str(z) + 'm dz: <'+ str(dz) + 'm correlated with ' + str(ind_months) + ' '+ str(osc)+' index \n grid avg: 8')
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/over_plots/'+str(osc)+'/corr1d_' + str(ind_months) +'lrggrid.png')

def mast_corr(lat, lon, data, corr_arr, index, reg, lat_sig, lon_sig, ind_months, cld_months, z, dz):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap = shiftedColorMap(cmap_send, 0, 0.6, 1, name='shifted')
    fig, ax = plt.subplots(2, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    ax1 = plt.subplot(2, 1, 1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.top_labels = False
    g1.right_labels = False
    ax1.coastlines(color='white')

    #p = ax1.contourf(lon, lat, corr_arr, cmap='RdBu', ylabel='Latitude')
    if reg == 'NEP':
        p = ax1.imshow(corr_arr, origin='lower', aspect='equal', extent=[-140,-116,10,37], cmap=shifted_cmap)
    if reg == 'SEP':
        p = ax1.imshow(corr_arr, origin='lower', aspect='equal', extent=[-100,-70,-40,-5], cmap=shifted_cmap)
    if reg == 'NEA':
        p = ax1.imshow(corr_arr, origin='lower', aspect='equal', extent=[-30,-10,20,40], cmap=shifted_cmap)
    if reg == 'GLB':
        p = ax1.imshow(corr_arr, origin='lower', aspect='equal', extent=[-180,180,-90,90], cmap='RdBu')
    if reg == 'TST':
        p = ax1.imshow(corr_arr, origin='lower', aspect='equal', extent=[-65,-55,-20,0], cmap='RdBu')
    if reg == 'SOP':
        p = ax1.imshow(corr_arr, origin='lower', aspect='equal', extent=[-140,-80,-63,-43], cmap=shifted_cmap)
    if reg == 'SA':
        p = ax1.imshow(corr_arr, origin='lower', aspect='equal', extent=[-50,2,-63,-43], cmap=shifted_cmap)
    c = plt.colorbar(p, shrink=0.5, ax=ax1)
    tstep = np.arange(0, lat_sig.shape[0])
    for t in tstep:
        ax1.scatter(lon_sig[t]+0.75, lat_sig[t]+1.25, color='black')
    #ax1.text(-0.10, 0.5, 'Latitude')
    #ax1.set_xlabel('Longitude')
    #ax1.set_ylabel('Latitude')

    ax2 = ax[1]

    #box_arr = np.transpose(box_arr)
    #**Needed to redo the labels and order of the data to have it be increasing order. Maybe should color based off of years.
    #Probably should ask Jay, also should make see if the contour or the imshow is better.
    #labels = ['-1.2', '-1.1', '1.2', '-2.1', '-0.2', '-0.2', '0.3', '1.9']
    # labels = ["0.4", "-0.7", "-0.95", "0.55", "-2.5", "0.2", "-0.65", "0.3", "1.8", "-0.35"]
    # labels_sorted = [labels[4],labels[2],labels[1],labels[6],labels[9],labels[5],labels[7],labels[0],labels[3],labels[8]]
    # box_arr_sorted = [box_arr[4],box_arr[2],box_arr[1],box_arr[6],box_arr[9],box_arr[5],box_arr[7],
    #                     box_arr[0],box_arr[3],box_arr[8]]
    # print(labels_sorted)
    test_index, test_arr = indexes.sort_index(index, data)
    #box_arr_sorted = test_arr.tolist()
    labels_sorted = test_index.tolist()
    #labels_sorted = str(labels_sorted)
    labels_sorted = [str(x) for x in labels_sorted]

    tstep = np.arange(0, data.shape[-1])
    for t in tstep:
        hold_arr = test_arr[t,:,:].flatten()
        if t == 0:
            box_arr = [hold_arr]
            continue
        box_arr.append(hold_arr)
    
    print(len(box_arr))
    print(len(labels_sorted))

    ax2.boxplot(box_arr, vert=True, patch_artist=True, labels=labels_sorted)
    ax2.set_ylabel('Cloud Occurrence Frequency')
    ax2.set_xlabel('MEI.v2 Index')
    ax2.set_title('Box and whisker plot of cloud occurrence in the study region for different ENSO indices')
    #plt.title('Covariance between cloud occurrence and MEI.v2 ENSO index')
    plt.xlabel('MEIv2 Index (sorted)')
    plt.ylabel('Cloud Occurrence Frequency')
    fig.suptitle(str(cld_months)+ ' Cloud Occurrence base z: <' + str(z) + 'm dz: <'+ str(dz) + 'm \n correlate with the '+ str(ind_months)+ ' ENSO index')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/mast_plots/mast_corr'+ str(reg)+ '_' + str(ind_months)+ '.png')

def oned_zonal(data):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap =shiftedColorMap(cmap_send, 1, 0.54, 0, 'shifted')
    fig, ax = plt.subplots(1, 1, figsize=(10,10))

    p = ax.imshow(data, origin='lower', aspect='auto', extent=[-90,90,0,20], cmap=shifted_cmap)
    c = plt.colorbar(p, shrink=0.5, ax=ax)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Height (km)')
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/zonal_plts/zonal_sam_JJA_corr.png')

def merra_1st_look(lat, lon, data_high, data_low):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap =shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    fig, ax = plt.subplots(3, 1, figsize =(10,20))
    proj = ccrs.Orthographic(180, -90)

    ax1 = plt.subplot(3, 1, 1, projection=proj)
    ax1.coastlines()
    p = ax1.contourf(lon, lat, data_high, transform=ccrs.PlateCarree(), cmap='jet')
    c = plt.colorbar(p, ax=ax1)
    ax1.set_title('High ENSO Years')

    ax2 = plt.subplot(3, 1, 2, projection=proj)
    ax2.coastlines()
    p = ax2.contourf(lon, lat, data_low, transform=ccrs.PlateCarree(), cmap='jet')
    c = plt.colorbar(p, ax=ax2)
    ax2.set_title('Low ENSO Years')

    ax3 = plt.subplot(3,1,3, projection=proj)
    ax3.coastlines()
    p = ax3.contourf(lon, lat, (data_high-data_low), transform=ccrs.PlateCarree(), cmap=shifted_cmap)
    c = plt.colorbar(p, ax=ax3)
    ax3.set_title('High ENSO Years - Low ENSO Years')
    
    fig.suptitle('JJ Mean 2 Meter Temperature')
    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/merra_plots/mb_500_polar.png')

def merra_mean(lat, lon, data_high, data_low, data_mean):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 =shiftedColorMap(cmap_send, 1, 0.6, 0, 'shifted')
    shifted_cmap2 =shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')
    fig, ax = plt.subplots(3, 1, figsize =(10,20))
    #proj = ccrs.SouthPolarStereo()
    proj=ccrs.Orthographic(180, -90)

    ax1 = plt.subplot(3, 1, 1, projection=proj)
    ax1.coastlines()
    p = ax1.contourf(lon, lat, data_mean, transform = ccrs.PlateCarree(), cmap='jet')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7%', pad=0.4, axes_class=maxes.Axes)
    cb = fig.colorbar(p, cax=cax, orientation='vertical')
    cax.set_ylabel('Height (m)')
    ax1.set_title('Climatological mean')

    ax2 = plt.subplot(3, 1, 2, projection=proj)
    ax2.coastlines()
    p = ax2.contourf(lon, lat, (data_low - data_mean), transform = ccrs.PlateCarree(), cmap=shifted_cmap1)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='7%', pad=0.4, axes_class=maxes.Axes)
    cb = fig.colorbar(p, cax=cax, orientation='vertical')
    cax.set_ylabel('height anomalies (m)')
    ax2.set_title('Low ENSO Years - Mean')

    ax3 = plt.subplot(3,1,3, projection=proj)
    ax3.coastlines()
    p = ax3.contourf(lon, lat, (data_high-data_mean), transform = ccrs.PlateCarree(), cmap=shifted_cmap2)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='7%', pad=0.4, axes_class=maxes.Axes)
    cb = fig.colorbar(p, cax=cax, orientation='vertical')
    cax.set_ylabel('height anomalies (m)')
    ax3.set_title('High ENSO Years - Mean')
    
    fig.suptitle('JJA Mean 500 mb height (m)')
    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/merra_plots/jjamean_FD_mb500_polar.png')


def cld_occ_sp(lat, lon, c_data, l_data, h_data, m_data):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    fig, ax = plt.subplots(3, 1, figsize=(10,20))
    proj = ccrs.Orthographic(180,-90)

    ax1 = plt.subplot(3,1,1, projection = proj)
    ax1.coastlines()
    ax1.set_extent([-179,179,-90,-45], crs=ccrs.PlateCarree())
    p = ax1.pcolormesh(lon, lat, c_data, transform=ccrs.PlateCarree(), cmap=shifted_cmap)
    c = plt.colorbar(p, ax=ax1)
    ax1.set_title('Correlation coefficient between ENSO and low cloud occurrence')

    ax2 = plt.subplot(3,1,2, projection =proj)
    shifted_cmap1 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    ax2.coastlines()
    ax2.set_extent([-179,179,-90,-45], crs=ccrs.PlateCarree())
    p = ax2.pcolormesh(lon, lat, (l_data-m_data)/m_data, transform=ccrs.PlateCarree(), cmap=shifted_cmap1, vmin=-0.2, vmax=0.2)
    c = plt.colorbar(p, ax=ax2)
    ax2.set_title('low ENSO years cloud occurrence - 2006-2016 mean cloud occurrence')

    ax3 = plt.subplot(3,1,3, projection =proj)
    shifted_cmap2 = shiftedColorMap(cmap_send, 0, 0.48, 1, 'shifted')
    ax3.coastlines()
    ax3.set_extent([-179,179,-90,-45], crs=ccrs.PlateCarree())
    p = ax3.pcolormesh(lon, lat, (h_data-m_data)/m_data, transform=ccrs.PlateCarree(), cmap=shifted_cmap2)
    c = plt.colorbar(p, ax=ax3)
    ax3.set_title('high ENSO years cloud occurrence - 2006-2016 mean cloud occurrence')

    fig.suptitle('Low cloud: 0-1000m base / 0-1500m thickness')
    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/merra_plots/cloud_overview2.png')

def merra_cld_occ_mast(m_lat, m_lon, mh_data, ml_data, mm_data, c_lat, c_lon, ch_data, cl_data, cm_data):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    fig, ax = plt.subplots(3,2, figsize=(20,20))
    proj = ccrs.Orthographic(180,-90)

    ax1 = plt.subplot(3,2,1, projection=proj)
    shifted_cmap1 = shiftedColorMap(cmap_send, 0, 0.58, 1, 'shifted')
    ax1.coastlines()
    p = ax1.contourf(m_lon, m_lat, (mm_data-mh_data), transform=ccrs.PlateCarree(), cmap=shifted_cmap1)
    c = plt.colorbar(p, ax=ax1)
    ax1.set_title('Climatological mean 500 mb height - High ENSO years')

    ax2 = plt.subplot(3,2,2, projection=proj)
    shifted_cmap2 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    ax2.coastlines()
    p = ax2.contourf(m_lon, m_lat, (mm_data-ml_data), transform=ccrs.PlateCarree(), cmap=shifted_cmap2)
    c = plt.colorbar(p, ax=ax2)
    ax2.set_title('Climatological mean 500 mb height - Low ENSO years')

    ax3 = plt.subplot(3,2,3, projection=proj)
    shifted_cmap3 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    ax3.coastlines()
    p = ax3.contourf(m_lon, m_lat, (mh_data-ml_data), transform=ccrs.PlateCarree(), cmap=shifted_cmap3)
    c = plt.colorbar(p, ax=ax3)
    ax3.set_title('High ENSO years 500 mb height - low ENSO years')

    ax4 = plt.subplot(3,2,4, projection=proj)
    shifted_cmap4 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    ax4.coastlines()
    ax4.set_extent([-179,179,-90,-45], crs=ccrs.PlateCarree())
    p = ax4.pcolormesh(c_lon, c_lat, (cm_data-ch_data), transform=ccrs.PlateCarree(), cmap=shifted_cmap4)
    c = plt.colorbar(p, ax=ax4)
    ax4.set_title('Mean cloud occurrence - high ENSO years')

    ax5 = plt.subplot(3,2,5, projection = proj)
    shifted_cmap5 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    ax5.coastlines()
    ax5.set_extent([-179,179,-90,-45], crs=ccrs.PlateCarree())
    p = ax5.pcolormesh(c_lon, c_lat, (cm_data-cl_data), transform=ccrs.PlateCarree(), cmap=shifted_cmap5, vmin=-0.3, vmax=0.3)
    c = plt.colorbar(p, ax=ax5)
    ax5.set_title('Mean cloud occurrence - low ENSO years')

    ax6 = plt.subplot(3,2,6, projection=proj)
    shifted_cmap6 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    ax6.coastlines()
    ax6.set_extent([-179,179,-90,-45], crs=ccrs.PlateCarree())
    p = ax6.pcolormesh(c_lon, c_lat, (ch_data-cl_data), transform=ccrs.PlateCarree(), cmap=shifted_cmap6, vmin=-0.3, vmax=0.3)
    c = plt.colorbar(p, ax=ax6)
    ax6.set_title('High ENSO years cloud occurrence - low ENSO years')
    fig.suptitle('500 mb height & 0-1000m base / 0-1500dz thickness cloud occurrence overview')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/merra_plots/mast_mb500.png')

def zonal_enso_comps(m_data, h_data, l_data, z, lat):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 0, 0.36, 1, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 0, 0.58, 1, 'shifted')
    shifted_cmap3 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    fig, ax = plt.subplots(4, 1, figsize=(10,20))
    
    ax1 = ax[0]
    p =ax1.imshow(m_data, origin='lower', extent=[lat[0], lat[-1], z[0], z[-1]], aspect='auto', cmap='magma')
    ax1.set_title('Mean zonal cloud occurrence')
    ax1.set_xlabel('Latitude')
    ax1.set_ylabel('Height (km)')
    c=plt.colorbar(p, ax=ax1)

    ax2 = ax[1]
    p=ax2.imshow((l_data-m_data)/m_data, origin='lower', extent=[lat[0], lat[-1], z[0], z[-1]], aspect='auto', cmap=shifted_cmap1)
    ax2.set_title('(Low ENSO years zonal mean cloud occurrence - mean zonal cloud occurrence)/mean zonal cloud occurrence')
    ax2.set_xlabel('Latitude')
    ax2.set_ylabel('Height (km)')
    c=plt.colorbar(p, ax=ax2)

    ax3 = ax[2]
    p=ax3.imshow((h_data-m_data)/m_data, origin='lower', extent=[lat[0], lat[-1], z[0], z[-1]], aspect='auto', cmap=shifted_cmap2)
    ax3.set_title('(High ENSO years zonal mean cloud occurrence - mean zonal cloud occurrence)/mean zonal cloud occurrence')
    ax3.set_xlabel('Latitude')
    ax3.set_ylabel('Height (km)')
    c=plt.colorbar(p, ax=ax3)

    ax4 = ax[3]
    p=ax4.imshow((h_data-l_data), origin='lower', extent=[lat[0], lat[-1], z[0], z[-1]], aspect='auto', cmap=shifted_cmap3)
    ax4.set_title('High ENSO years zonal mean cloud occurrence - low ENSO years zonal mean cloud occurrence')
    ax4.set_xlabel('Latitude')
    ax4.set_ylabel('Height (km)')
    c=plt.colorbar(p, ax=ax4)

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/enso_zonal_comp_rel.png')

def test_trend(years, data):
    plt.scatter(years, data)
    plt.xlabel('Years')
    plt.ylabel('PBL top pressure (Pa)')
    plt.title('PBL top pressure at lat: -60 lon: 100')
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/tester.png')

def eof_inital(lat, lon, data):
    proj = ccrs.Orthographic(180,-90)
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    ax1 = plt.subplot(1,1,1, projection=proj)
    p = ax1.contourf(lon, lat, data, transform=ccrs.PlateCarree(), cmap='magma')
    plt.colorbar(p, ax=ax1)

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/tester.png')

def sam_enso_rel(p_data, n_data, mean):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.68, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.6, 0, 'shifted')
    shifted_cmap3 = shiftedColorMap(cmap_send, 1, 0.38, 0, 'shifted')
    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(3,1, figsize=(10,14))
    ax1 = plt.subplot(3,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines()
    #p = ax1.imshow(p_data, origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap='magma')
    p = ax1.imshow((p_data-mean), origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap=shifted_cmap2)
    c = plt.colorbar(p, ax=ax1)
    ax1.set_title('Positive relationship years (2013-2015) - mean low cloud occurrence')

    ax2 = plt.subplot(3,1,2, projection=proj)
    g2 = ax2.gridlines(draw_labels = True, linestyle = 'dotted')
    g2.xlabels_top = False
    g2.ylabels_right = False
    ax2.coastlines()
    #p = ax2.imshow(n_data, origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap='magma')
    p = ax2.imshow((n_data-mean), origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap=shifted_cmap3)
    c = plt.colorbar(p, ax=ax2)
    ax2.set_title('Negative relationship years (2008-2010) - mean low cloud occurrence')

    ax3 = plt.subplot(3,1,3, projection=proj)
    g3 = ax3.gridlines(draw_labels = True, linestyle = 'dotted')
    g3.xlabels_top = False
    g3.ylabels_right = False
    ax3.coastlines()
    p = ax3.imshow((p_data - n_data), origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap=shifted_cmap1)
    c = plt.colorbar(p, ax=ax3)
    ax3.set_title('Positive relationship years - negative relationship years')

    #plt.suptitle('low cloud: z:<1000m dz:<1500m')
    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/index_plts/low_cld_diff_pn_2008_2008_2013_2014.png')

def sam_enso_timep(cld_mean, cld_ear, cld_mid, cld_late):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.43, 0, 'shifted')
    shifted_cmap3 = shiftedColorMap(cmap_send, 1, 0.54, 0, 'shifted')

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(3,1, figsize=(12,14))
    ax1 = plt.subplot(3,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines()
    p = ax1.imshow((cld_ear-cld_mean), origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap=shifted_cmap1)
    plt.colorbar(p, ax=ax1)
    ax1.set_title('Early years (2006-2010) - mean')
    
    ax2 = plt.subplot(3,1,2, projection=proj)
    g2 = ax2.gridlines(draw_labels = True, linestyle='dotted')
    g2.xlabels_top = False
    g2.ylabels_right = False
    ax2.coastlines()
    p = ax2.imshow((cld_mid-cld_mean), origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap=shifted_cmap2)
    plt.colorbar(p, ax=ax2)
    ax2.set_title('Mid years (2010-2014) - mean')

    ax3 = plt.subplot(3,1,3, projection=proj)
    g3 = ax3.gridlines(draw_labels = True, linestyle='dotted')
    g3.xlabels_top = False
    g3.ylabels_right = False
    ax3.coastlines()
    p = ax3.imshow((cld_late-cld_mean), origin='lower', extent=[-180,180,-90,90], aspect='equal', cmap=shifted_cmap3)
    plt.colorbar(p, ax=ax3)
    ax3.set_title('Late years (2014-2017) - mean')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/index_plts/time_period_thrds_SON.png')

def sam_enso_timep_so(lat, lon, cld_mean, cld_ear, cld_mid, cld_late):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.38, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.48, 0, 'shifted')
    shifted_cmap3 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')

    fig, ax = plt.subplots(3, 1, figsize=(10,20))
    proj = ccrs.Orthographic(180,-90)

    ax1 = plt.subplot(3,1,1, projection = proj)
    ax1.coastlines()
    ax1.set_extent([-179,179,-90,-30], crs=ccrs.PlateCarree())
    p = ax1.pcolormesh(lon, lat, (cld_ear-cld_mean), transform=ccrs.PlateCarree(), cmap=shifted_cmap1)
    c = plt.colorbar(p, ax=ax1)
    ax1.set_title('Early years (2006-2010) - mean *Months:all')

    ax2 = plt.subplot(3,1,2, projection =proj)
    shifted_cmap1 = shiftedColorMap(cmap_send, 0, 0.5, 1, 'shifted')
    ax2.coastlines()
    ax2.set_extent([-179,179,-90,-30], crs=ccrs.PlateCarree())
    p = ax2.pcolormesh(lon, lat, (cld_mid-cld_mean), transform=ccrs.PlateCarree(), cmap=shifted_cmap2)
    c = plt.colorbar(p, ax=ax2)
    ax2.set_title('Mid years (2010-2014) - mean *Months:all')

    ax3 = plt.subplot(3,1,3, projection =proj)
    shifted_cmap2 = shiftedColorMap(cmap_send, 0, 0.48, 1, 'shifted')
    ax3.coastlines()
    ax3.set_extent([-179,179,-90,-30], crs=ccrs.PlateCarree())
    p = ax3.pcolormesh(lon, lat, (cld_late-cld_mean), transform=ccrs.PlateCarree(), cmap=shifted_cmap3)
    c = plt.colorbar(p, ax=ax3)
    ax3.set_title('Late years (2014-2017) - mean *Months:all')

    fig.suptitle('Low cloud: 0-1000m base / 0-1500m thickness')
    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/index_plts/time_period_thrds_SO')

    

def merra_sam_enso(lat, lon, pos_data, neg_data, mean_data):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.53, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.4, 0, 'shifted')

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(2,1, figsize=(10,10))
    ax1 = plt.subplot(2,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines()
    p = ax1.contourf(lon, lat, (pos_data - mean_data), cmap=shifted_cmap1)
    plt.colorbar(p, ax=ax1)
    ax1.set_title('Positive correlation years SLP (Pa) - mean SLP (Pa)')

    ax2 = plt.subplot(2,1,2, projection=proj)
    g2 = ax2.gridlines(draw_labels = True, linestyle = 'dotted')
    g2.xlabels_top = False
    g2.ylabels_right = False
    ax2.coastlines()
    p = ax2.contourf(lon, lat, (neg_data - mean_data), cmap=shifted_cmap2)
    plt.colorbar(p, ax=ax2)
    ax2.set_title('Negative correlation years SLP (Pa) - mean SLP (Pa)')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/merra_plots/sam_enso_corr_plts/all_periods_avg.png')

def index_histo(data, osc, months_name):
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    ax.hist(data)
    ax.set_title('Histogram of ' + osc +' index values for months: ' + months_name)
    
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/index_plts/index_histograms/'+osc+'_'+months_name+'.png')

def merra_corr(lat, lon, mean, p1, p2, p3, p4, p5, p6):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.53, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.53, 0, 'shifted')
    shifted_cmap3 = shiftedColorMap(cmap_send, 1, 0.53, 0, 'shifted')
    shifted_cmap4 = shiftedColorMap(cmap_send, 1, 0.53, 0, 'shifted')
    shifted_cmap5 = shiftedColorMap(cmap_send, 1, 0.63, 0, 'shifted')
    shifted_cmap6 = shiftedColorMap(cmap_send, 1, 0.53, 0, 'shifted')

    fig, ax = plt.subplots(3,2, figsize=(10,10))
    proj = ccrs.PlateCarree()

    ax1 = plt.subplot(3,2,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines()
    p=ax1.contourf(lon, lat, (p1-mean), cmap=shifted_cmap1)
    plt.colorbar(p, ax=ax1)
    ax1.set_title('(1980-1984) mean - Climatological mean')

    ax2 = plt.subplot(3,2,2, projection=proj)
    g2 = ax2.gridlines(draw_labels = True, linestyle = 'dotted')
    g2.xlabels_top = False
    g2.ylabels_right = False
    ax2.coastlines()
    p=ax2.contourf(lon, lat, (p2-mean), cmap=shifted_cmap2)
    plt.colorbar(p, ax=ax2)
    ax2.set_title('(1984-1993) mean - climatological mean')

    ax3 = plt.subplot(3,2,3, projection=proj)
    g3 = ax3.gridlines(draw_labels = True, linestyle = 'dotted')
    g3.xlabels_top = False
    g3.ylabels_right = False
    ax3.coastlines()
    p=ax3.contourf(lon, lat, (p3-mean), cmap=shifted_cmap3)
    plt.colorbar(p, ax=ax3)
    ax3.set_title('(1993-1997) mean - climatological mean')

    ax4 = plt.subplot(3,2,4, projection=proj)
    g4 = ax4.gridlines(draw_labels = True, linestyle = 'dotted')
    g4.xlabels_top = False
    g4.ylabels_right = False
    ax4.coastlines()
    p=ax4.contourf(lon, lat, (p4-mean), cmap=shifted_cmap4)
    plt.colorbar(p, ax=ax4)
    ax4.set_title('(1997-2012) mean - climatological mean')

    ax5 = plt.subplot(3,2,5, projection=proj)
    g5 = ax5.gridlines(draw_labels = True, linestyle = 'dotted')
    g5.xlabels_top = False
    g5.ylabels_right = False
    ax5.coastlines()
    p=ax5.contourf(lon, lat, (p5-mean), cmap=shifted_cmap5)
    plt.colorbar(p, ax=ax5)
    ax5.set_title('(2012-2017) mean - climatological mean')

    ax6 = plt.subplot(3,2,6, projection=proj)
    g6 = ax6.gridlines(draw_labels = True, linestyle = 'dotted')
    g6.xlabels_top = False
    g6.ylabels_right = False
    ax6.coastlines()
    p=ax6.contourf(lon, lat, (p6-mean), cmap=shifted_cmap6)
    plt.colorbar(p, ax=ax6)
    ax6.set_title('(2017-2020) mean - climatological mean')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/merra_plots/sam_enso_corr_plts/allmonths_1980_2020.png')


def make_weight_plot(data, lat_r, level_r, day):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')

    fig, ax = plt.subplots(1,1,figsize=(10,10))

    norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=data.min(), vmax=data.max())
    #p = ax.contourf(lat_r, level_r, data, cmap=shifted_cmap1, norm=norm)
    p = ax.imshow(data, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap=shifted_cmap1, norm=norm)
    plt.colorbar(p,ax=ax)
    #ax.set_title('Weighted Anomalous U data. Day=06/01/1985 + ' + str(day))
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Pressure (hPa)')

    if day == 1:
        ax.set_title('Anomalous U data. Day=06/01/1985 + ' + str(day))
        plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg/pre_weight.png')
    elif day == 2: 
        ax.set_title('Weighted Anomalous U data. Day=06/01/1985 + ' + str(day))
        plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg/post_weight.png')
    plt.close()

def show_eke(data1, data2, level_r, lon_s):
    fig, ax = plt.subplots(2,1,figsize=(10,10))

    ax1 = ax[0]
    p = ax1.imshow(data1, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
    plt.colorbar(p, ax=ax1)
    ax1.set_xlabel('Latitude')
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_title('Cross Section of EKE at Longitude '+ str(lon_s))

    ax2 = ax[1]
    p = ax2.imshow(data2, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
    plt.colorbar(p, ax=ax2)
    ax2.set_xlabel('Latitude')
    ax2.set_ylabel('Pressure (hPa)')
    ax2.set_title('Zonal Mean EKE')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg/eke_test_w_lon.png')
    plt.close()

def eke_plot_loop(data_eke, data_eke_m, td_eke, level_r, lat_r, lon_r, day):
    fig, ax = plt.subplots(3,1, figsize=(10,10))
    proj_so = ccrs.Orthographic(180,-90)
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')

    ax1 = ax[0]
    #norm1 = norm_maker(data_eke)
    p = ax1.imshow(data_eke, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
    plt.colorbar(p, ax=ax1)
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_xlabel('Latitude')
    ax1.set_title('Zonal Mean EKE for date=01/01/2008 + ' + str(day*4) + ' days')

    ax2 = ax[1]
    #norm2 = norm_maker(data_eke_m)
    p = ax2.imshow(data_eke_m, origin='upper',extent=[-70,-20,level_r[-1],level_r[0]],aspect='auto',cmap='jet')
    plt.colorbar(p, ax=ax2)
    ax2.set_ylabel('Pressure (hPa)')
    ax2.set_xlabel('Latitude')
    ax2.set_title('Rolling mean of zonal mean EKE data')

    td_eke_r = np.transpose(td_eke)
    td_eke_r = np.flipud(td_eke_r)
    lat_f = np.flipud(lat_r)
    ax3 = plt.subplot(3,1,3,projection=proj_so)
    ax3.coastlines(color='white')
    ax3.set_extent([-179,179,-90,-15], crs=ccrs.PlateCarree())
    #norm3 = norm_maker(td_eke_r)
    p = ax3.pcolormesh(lon_r, lat_f, td_eke_r, transform=ccrs.PlateCarree(), cmap='jet')
    plt.colorbar(p, ax=ax3)
    # ax3.set_xlabel('Longitude')
    # ax3.set_ylabel('Latitude')
    ax3.set_title('500 mb EKE for date=01/01/2008 + ' + str(day*4) + ' days')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/eke_3d_1986_mw/'+str(day)+'.png')
    plt.close()

def zonal_eof(data_send, lat_bins, z_bins, pc1, enso_indx):
    fig, ax = plt.subplots(2,1, figsize=(10,10))
    ax1 = ax[0]
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')

    norm = norm_maker(data_send)
    #p = plt.contourf(np.arange(0,45),np.arange(0,21),corr_arr, cmap='RdBu')
    p = ax1.imshow(data_send, cmap=shifted_cmap1, norm=norm, extent=[lat_bins[-1],lat_bins[0],z_bins[0],z_bins[-1]], origin='upper', aspect='auto')
    plt.colorbar(p, ax=ax1)
    ax1.set_xlabel('Latitude')
    ax1.set_ylabel('Height (km)')
    ax1.set_title('Regression of EOF2 with zonal mean cloud occurrence')

    ax2 = ax[1]
    ax2.plot(np.arange(0,pc1.shape[0]), pc1, color='Red', label='PC2')
    ax2.plot(np.arange(0,pc1.shape[0]), enso_indx, color='Blue', label='ENSO Index')
    ax2.set_xlabel('Time (months of 2008)')
    ax2.set_ylabel('Standardized index (PC2)')
    ax2.set_title('PC1 of zonal mean cloud occurrence')
    ax2.legend()
    

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/cld_occ_zonal_eof2.png')

def bam_eke_corr(corr_s, reg_ins_s, reg_slp_s):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.45, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.42, 0, 'shifted')
    fig, ax = plt.subplots(3,1, figsize=(10,10))

    ax1 = ax[0]
    p= ax1.imshow(corr_s, extent=[-70,-20,1000,200], origin='upper', aspect='auto', cmap='jet')
    plt.colorbar(p, ax=ax1)
    ax1.set_xlabel('Latitude')
    ax1.set_ylabel('Height (hPa)')
    ax1.set_title('BAM PC1 correlated with EKE (2006-2009)')
    
    ax2 = ax[1]
    p = ax2.imshow(reg_ins_s, extent=[-70,-20,1000,200], origin='upper', aspect='auto', cmap='jet')
    plt.colorbar(p, ax=ax2)
    ax2.set_xlabel('Latitude')
    ax2.set_ylabel('Height (hPa)')
    ax2.set_title('BAM PC1 / EKE (2006-2009), regression intercept')
    
    ax3 = ax[2]
    p = ax3.imshow(reg_slp_s, extent=[-70,-20,1000,200], origin='upper', aspect='auto', cmap='jet')
    plt.colorbar(p, ax=ax3)
    ax3.set_xlabel('Latitude')
    ax3.set_ylabel('Height (hPa)')
    ax3.set_title('BAM PC1 / EKE (2006-2009), regression slope')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg/bam_eof1_eke_mast_t2.png')

def bam_cld_corr(corr_s, reg_ins_s, reg_slp_s):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    cmap_send1 = matplotlib.cm.get_cmap('PiYG')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.28, 0, 'shifted')
    shifted_cmap3 = shiftedColorMap(cmap_send, 1, 0.49, 0, 'shifted')

    fig, ax = plt.subplots(3,1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    
    ax1 = plt.subplot(3,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines()
    p = ax1.imshow(corr_s, extent=[-180,180,-90,90], origin='lower', aspect='auto', cmap=shifted_cmap1)
    plt.colorbar(p, ax = ax1)
    ax1.set_title('Monthly averaged BAM index lag-1, correlated with low cloud occurrence MAM (2007-2009)')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = plt.subplot(3,1,2, projection=proj)
    g2 = ax2.gridlines(draw_labels = True, linestyle = 'dotted')
    g2.xlabels_top = False
    g2.ylabels_right = False
    ax2.coastlines()
    p = ax2.imshow(reg_ins_s, extent=[-180,180,-90,90], origin='lower',aspect='auto',cmap=shifted_cmap2)
    plt.colorbar(p, ax=ax2)
    ax2.set_title('Monthly averaged BAM index lag-1, regressed with low cloud occurrence, intercept MAM (2007-2009)')
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')

    ax3 = plt.subplot(3,1,3, projection=proj)
    g3 = ax3.gridlines(draw_labels = True, linestyle = 'dotted')
    g3.xlabels_top = False
    g3.ylabels_right = False
    ax3.coastlines()
    p = ax3.imshow(reg_slp_s, extent=[-180,180,-90,90], origin='lower',aspect='auto',cmap=shifted_cmap3)
    plt.colorbar(p, ax=ax3)
    ax3.set_title('Monthly averaged BAM index lag-1, regressed with low cloud occurrence, slope MAM (2007-2009)')
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg_cld_occ/lng_ts_lag-1_mam_sm.png')

def bam_stats_package(p_vals, r_vals):
    fig, ax = plt.subplots(2, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()

    ax1 = plt.subplot(2,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines()
    p = ax1.imshow(r_vals, extent=[-180,180,-90,90], origin='lower', aspect='auto', cmap='jet')
    plt.colorbar(p, ax = ax1)
    ax1.set_title('R^2 values of BAM index correlated with low cloud occurrence (2006-2009)')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = plt.subplot(2,1,2, projection=proj)
    g2 = ax2.gridlines(draw_labels = True, linestyle = 'dotted')
    g2.xlabels_top = False
    g2.ylabels_right = False
    ax2.coastlines()
    p = ax2.imshow(p_vals, extent=[-180,180,-90,90], origin='lower',aspect='auto',cmap='jet', vmin=0.01, vmax=0.1)
    plt.colorbar(p, ax=ax2)
    ax2.set_title('P values of BAM index regressed with low cloud occurrence (2006-2009)')
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg/stats_pkg_storms.png')
    

def bam_cld_corr_zonal(corr_s, reg_ins_s, reg_slp_s):
    cmap_send = matplotlib.cm.get_cmap('RdBu')
    shifted_cmap1 = shiftedColorMap(cmap_send, 1, 0.48, 0, 'shifted')
    shifted_cmap2 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')
    shifted_cmap3 = shiftedColorMap(cmap_send, 1, 0.5, 0, 'shifted')

    fig, ax = plt.subplots(3,1, figsize=(10,10))

    ax1 = ax[0]
    p = ax1.imshow(corr_s, origin='upper', extent=[-90,90,0,20], aspect='auto', cmap=shifted_cmap1)
    plt.colorbar(p, ax=ax1)
    ax1.set_xlabel('Latitude')
    ax1.set_ylabel('Height (km)')
    ax1.set_title('Monthly averaged BAM index correlated with anomalous zonal mean cloud occurrence (2006-2009)')

    ax2 = ax[1]
    p = ax2.imshow(reg_ins_s, origin='upper', extent=[-90,90,0,20], aspect='auto', cmap=shifted_cmap2)
    plt.colorbar(p, ax=ax2)
    ax2.set_xlabel('Latitude')
    ax2.set_ylabel('Height (km)')
    ax2.set_title('Monthly averaged BAM index regressed with anomalous zonal mean cloud occurrence (2006-2009), intercept')

    ax3 = ax[2]
    p = ax3.imshow(reg_slp_s, origin='upper', extent=[-90,90,0,20], aspect='auto', cmap=shifted_cmap3)
    plt.colorbar(p, ax=ax3)
    ax3.set_xlabel('Latitude')
    ax3.set_ylabel('Height (km)')
    ax3.set_title('Monthly averaged BAM index regressed with anomalous zonal mean cloud occurrence (2006-2009), slope')

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/zonal_bam_corr.png')

def seasonal_cycle_test(data_anom_mean, data_mean):
    fig, ax = plt.subplots(2,1, figsize=(10,10))

    ax1 = ax[0]
    ax1.plot(np.arange(0,data_anom_mean.shape[0]), data_anom_mean, label='Anomalous cloud occurrence', color='red')
    ax1.set_ylabel('Cloud occurrence')
    ax1.set_xlabel('Time, months since 6/2006')
    ax1.set_title('Anomalous (Seasonal mean removed)')

    ax2 = ax[1]
    ax2.plot(np.arange(0,data_mean.shape[0]), data_mean, label='Mean cloud occurrence', color='blue')
    ax2.set_ylabel('Cloud occurrence')
    ax2.set_xlabel('Time, months since 6/2006')
    ax2.set_title('Raw data')

    plt.legend()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/cloud_occ_ts_zscore.png')
    
def bam_ts_plot(raw_s, mm_s, anom_s, mm_anom_s):
    fig, ax = plt.subplots(2,1,figsize=(10,10))
    ax1 = ax[0]
    p1, =ax1.plot(np.arange(0,raw_s.shape[0]), raw_s, label='Raw BAM data')
    ax1t = ax1.twiny()
    p2, =ax1t.plot(np.arange(0,mm_s.shape[0]), mm_s, color='orange', label='Monthly mean raw data')
    ax1.set_xlabel('Time, days + 1/1/2006')
    ax1t.set_xlabel('Time, months + 1/2006')
    ax1.set_ylabel('Standardized index')
    ax1.set_title('PC1 of EKE')

    ax2 = ax[1]
    p3, =ax2.plot(np.arange(0,anom_s.shape[0]), anom_s, label='Anomalous Raw BAM data')
    ax2t = ax2.twiny()
    p4, =ax2t.plot(np.arange(0,mm_anom_s.shape[0]), mm_anom_s, color='orange', label='monthly mean anomalous data')
    ax2.set_xlabel('Time, days + 1/1/2006')
    ax2t.set_xlabel('Time, months + 1/2006')
    ax2.set_ylabel('Standardized index')

    ax1.legend(handles=[p1, p2])
    ax2.legend(handles=[p3, p4])

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/timeseries_mast_sr1.png')

def lin_reg_func(slp_s, ints_s, x_shape):
    """
    This function will return an array of y values based on the linear regression of shape x_shape.
    !x_shape must be int shape.
    """
    tstep = np.arange(-10, x_shape)
    for t in tstep:
        y_h = slp_s*t+ints_s
        if t == tstep[0]:
            y_arr = np.array([y_h])
            continue
        y_arr = np.append(y_arr, y_h)

    return y_arr

def make_ts_bam_cld_occ(x_send, y_send, slp, ints, corr, p_val, lat_s, lon_s):
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    
    p1, =ax.plot(np.arange(0,x_send.shape[0]), x_send, color='blue', label='BAM index')
    axt = ax.twinx()
    p2, =axt.plot(np.arange(0,y_send.shape[0]), y_send, color='red', label='Cloud Occurrence')
    ax.legend(handles=[p1,p2])
    ax.set_ylabel('Standardized index')
    axt.set_ylabel('Anomalous cloud occurrence')
    ax.set_xlabel('Time, months + 1/2006')
    ax.set_title('Time series of cloud occurrence and bam index for point: '+ str(lat_s)+ ', '+ str(lon_s)+
                 '\nRegression slope: '+ str(slp)+ ' intercept: '+ str(ints)+ ' Corr coef.: '+ str(corr))

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/onep_ts.png')

def std_err_bam(std_err):
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    proj = ccrs.PlateCarree()

    ax1 = plt.subplot(1,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines(color='white')
    p = ax1.imshow(std_err, extent=[-180,180,-90,90], origin='lower', aspect='equal', cmap='jet')
    plt.colorbar(p, ax = ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.set_title('Standard Error Values for DJF (2006/2007-2008/2009)')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/std_err_djf.png')
    plt.close()

def one_point_anal_bsc(x_send, y_send, slp, ints, corr, p_val, lat_s, lon_s):
    reg_ln = lin_reg_func(slp, ints, x_send.shape[0])
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    
    ax.scatter(x_send, y_send)
    ax.plot(np.arange(-10,x_send.shape[0]), reg_ln, color='red', label='Linear Regression')
    ax.vlines(0,np.min(y_send),np.max(y_send), color='black')
    ax.set_xlabel('PC1 of EKE (Monthly mean BAM index)')
    ax.set_ylabel('Anomalous cloud occurrence')
    ax.set_title('Regression slope: '+ str(slp)+ ' intercept: '+ str(ints)
                 +'\n at point Lat: '+ str(lat_s)+ ' Lon: '+ str(lon_s)+ '\n Correlation coefficient: '+ str(corr)+ ' P-value: '+ str(p_val))

    ax.set_xlim([np.min(x_send),np.max(x_send)])
    ax.set_ylim([np.min(y_send),np.max(y_send)])

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/onep_test.png')

def one_point_anal(x_send, y_send, slp, ints, lat_s, lon_s, p_val, corr):
    reg_ln = lin_reg_func(slp, ints, x_send.shape[0])
    proj = ccrs.PlateCarree()

    fig, ax = plt.subplots(2,1, figsize=(10,10))

    ax1 = ax[0]
    ax1.scatter(x_send, y_send)
    ax1.plot(np.arange(-10,x_send.shape[0]), reg_ln, color='red', label='Linear Regression')
    ax1.vlines(0,np.min(y_send),np.max(y_send), color='black')
    ax1.set_xlabel('PC1 of EKE (Monthly mean BAM index)')
    ax1.set_ylabel('Anomalous cloud occurrence')
    ax1.set_title('Regression slope: '+ str(slp)+ ' intercept: '+ str(ints)
                 +'\n at point Lat: '+ str(lat_s)+ ' Lon: '+ str(lon_s)+ '\n Correlation coefficient: '+ str(corr)+ ' P-value: '+ str(p_val))

    ax1.set_xlim([np.min(x_send),np.max(x_send)])
    ax1.set_ylim([np.min(y_send),np.max(y_send)])

    ax2 = plt.subplot(2,1,2, projection = proj)
    g2 = ax2.gridlines(draw_labels = True, linestyle = 'dotted')
    g2.xlabels_top = False
    g2.ylabels_right = False
    ax2.coastlines(color='black')
    ax2.scatter(lon_s, lat_s)
    ax2.set_xlim([-180,180])
    ax2.set_ylim([-90,90])

    ax2.set_title('Location of grid select')
    # plt.xlim([-1,1])
    # plt.ylim([0,50])
    plt.tight_layout()
    #plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/oned_bam_gridstep_t0/'+ str(lat_s)+ '_'+ str(lon_s)+ '.png')
    plt.close()

def test_bam_nao_rel(bam_idx, nao_idx, slp, intc, pval, corr_coef):
    reg_lin = lin_reg_func(slp, intc, nao_idx.shape[0])
    fig, ax = plt.subplots(2,1,figsize=(10,10))

    ax1 = ax[0]
    l1, = ax1.plot(np.arange(0,bam_idx.shape[0]), bam_idx, color='red', label = 'Monthly Mean BAM index')
    l2, = ax1.plot(np.arange(0,nao_idx.shape[0]), nao_idx, color='blue', label = 'Monthly Mean NAO index')
    ax1.set_xlabel('Time, months + 1/2006')
    ax1.set_ylabel('Standardized index values')
    ax1.set_title('Time series of anomalous NAO and BAM index')
    ax1.legend(handles=[l1,l2])

    ax2 = ax[1]
    l3, = ax2.plot(np.arange(-10,nao_idx.shape[0]), reg_lin, color='red', label='Regression Line')
    ax2.scatter(nao_idx, bam_idx)
    ax2.vlines(0,np.min(bam_idx),np.max(bam_idx), color='black')
    ax2.set_xlim([np.min(nao_idx),np.max(nao_idx)])
    ax2.set_ylim([np.min(bam_idx),np.max(bam_idx)])

    ax2.set_xlabel('Anomalous NAO Standardized Index Values')
    ax2.set_ylabel('Anomalous BAM Standardized Index Values')
    ax2.set_title('Linear Regression\nRegression slope: '+ str(slp)+ ' intercept: '+ str(intc)+ 
                  '\n Correlation coefficient: '+ str(corr_coef))
    ax2.legend(handles=[l3])

    plt.tight_layout()
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_nao_ts_anom.png')
    
def test_svd(data, vh, u):
    fig, ax = plt.subplots(3, 1, figsize=(10,10))

    ax1 = ax[0]
    p = ax1.imshow(data[0,:,:], cmap='RdBu')
    plt.colorbar(p, ax=ax1)

    ax2 = ax[1]
    p = ax2.imshow(vh[0,:,:], cmap='RdBu')
    plt.colorbar(p, ax=ax2)

    ax3 = ax[2]
    ax3.plot(np.arange(0,u.shape[0]), u[0,:])

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/tester.png')

def raw_cld_occ_loop(data, year_s):
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    proj = ccrs.PlateCarree()

    ax1 = plt.subplot(1,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines(color='black')
    norm = norm_maker(data)
    p = ax1.imshow(data, origin='lower', aspect='equal', extent=[-180,180,-90,90], cmap='RdBu', norm=norm)
    plt.colorbar(p, ax = ax1)
    ax1.set_title('Year: '+ str(year_s))

    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/raw_cld_bam_loop/'+str(year_s)+'.png')
    plt.close()

def bam_mean_cld_occ(data_s):
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    proj = ccrs.PlateCarree()

    ax1 = plt.subplot(1,1,1, projection=proj)
    g1 = ax1.gridlines(draw_labels = True, linestyle = 'dotted')
    g1.xlabels_top = False
    g1.ylabels_right = False
    ax1.coastlines(color='black')
    #norm = norm_maker(data_s)
    p = ax1.imshow(data_s, origin='lower', aspect='equal', extent=[-180,180,-90,90], cmap='jet')
    plt.colorbar(p, ax=ax1)
    ax1.set_title('Mean cloud occurrence (2006-2016) bad data removed')
    
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/mean_cld_occ_bam.png')

if __name__ == '__main__':
    data_id = '08'
    print('main code running')
    #*Uses get_data.py to get/test data.
    get_files_dict = get_data.find_files(data_id, 0)
    get_data.test_monthly_data(get_files_dict)
    month_sum_data = sum_months()
    data_send = month_sum_data[1,:,:,:,:]
    if make_plot == 1:
        make_plot_main()
