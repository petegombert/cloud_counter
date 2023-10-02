"""
*This code builds the long-term mean reanalysis values.

To do:

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
import merra_first
import merra_analize
import indexes
import go_stats
import mast_plot
import overview

high_enso_years = ['2006','2009','2012','2015']
low_enso_years = ['2013','2014','2016']
mean = ['2006','2007','2008','2009','2010','2012','2013','2014','2015','2016']
months = ['06','07','08']


def save_year(lat_s, lon_s, mb_500_s, slp_s, t2m_s, pbltop_s, q500_s, year, cat):
    #os.chdir('/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/FD_ENSO/'+str(cat))
    os.chdir('/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/SAM_GLB/'+str(cat))
    f = cdf.Dataset('keep'+str(year)+'.nc4', 'w', format='NETCDF4')
    lat = f.createDimension('lat')
    lon = f.createDimension('lon')

    lat = f.createVariable('lat', 'f4', ('lat',))
    lon = f.createVariable('lon', 'f4', ('lon',))
    mb500 = f.createVariable('mb500', 'f4', ('lat','lon',))
    slp = f.createVariable('slp', 'f4', ('lat','lon',))
    t2m = f.createVariable('t2m', 'f4', ('lat','lon',))
    pbltop = f.createVariable('pbltop', 'f4', ('lat','lon',))
    q500 = f.createVariable('q500', 'f4', ('lat','lon',))

    print(mb_500_s.shape)

    lat[:] = lat_s
    lon[:] = lon_s
    slp[:] = slp_s
    mb500[:] = mb_500_s
    t2m[:] = t2m_s
    pbltop[:] = pbltop_s
    q500[:] = q500_s
    f.close()

def import_year(year):
    tstep = np.arange(6,9)
    for t in tstep:
        if t < 10:
            month = go_stats.add_zeros(str(t), 1)
        else:
            month = str(t)
        merra_first.merra_call_get(year, month, 'jja_mean', [-40,-90,-180,180])
    mast_arr = merra_analize.build_ds('jja_mean')
    print(mast_arr[0,2].shape)
    lat = mast_arr[0,0]
    lon = mast_arr[0,1]
    m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mast_arr)
    save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, year, 'jja_mean')

#Run this to build the mean dataset. Need to add the variables in all the other codes. 
# if __name__ == '__main__':
#    tstep = np.arange(2006,2017,1)
#    for t in tstep:
#        import_year(str(t))

#** This code builds the time series data for detrending
# if __name__ == '__main__':
#     years = np.arange(1984, 2020)
#     slp_ts = merra_analize.build_ts(years, 'jja_mean', -55, 0)
#     mast_plot.test_trend(years, slp_ts)

#*This will run the initial code for ENSO connection (6-23)
# if __name__ == '__main__':
#     m_mast_arr = merra_analize.build_ds('jja_mean')
#     h_mast_arr = merra_analize.build_ds('high')
#     l_mast_arr = merra_analize.build_ds('low')
#     m_lat = m_mast_arr[0,0]
#     m_lon = m_mast_arr[0,1]
#     m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(m_mast_arr)
#     h_mb500, h_slp, h_t2m, h_pbltop, h_q500 = merra_analize.get_mean_vals(h_mast_arr)
#     l_mb500, l_slp, l_t2m, l_pbltop, l_q500 = merra_analize.get_mean_vals(l_mast_arr)

#     mast_plot.merra_mean(m_lat, m_lon, h_mb500, l_mb500, m_mb500)
    #mast_plot.merra_mean(m_lat, m_lon, (h_slp-h_pbltop), (l_slp-l_pbltop), (m_slp-m_pbltop))
#    exit()
#     c_lat, c_lon, cld_he = overview.call_main(high_enso_years, months, 0)
#     c_lat, c_lon, cld_le = overview.call_main(low_enso_years, months, 0)
#     c_lat, c_lon, cld_mean = overview.call_main(mean, months, 0)

#     mast_plot.merra_cld_occ_mast(m_lat, m_lon, h_mb500, l_mb500, m_mb500, c_lat[0][1:], c_lon[0][1:], 
#                                  cld_he, cld_le, cld_mean)

#*This will run the main code for SAM/ENSO correlation question (7-23)
if __name__ == '__main__':
    m_mast_arr = merra_analize.build_ds('mean')
    p_mast_arr = merra_analize.build_ds('pos_years')
    n_mast_arr = merra_analize.build_ds('neg_years')
    m_lat = m_mast_arr[0,0]
    m_lon = m_mast_arr[0,1]
    m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(m_mast_arr)
    p_mb500, p_slp, p_t2m, p_pbltop, p_q500 = merra_analize.get_mean_vals(p_mast_arr)
    n_mb500, n_slp, n_t2m, n_pbltop, n_q500 = merra_analize.get_mean_vals(n_mast_arr)

    mast_plot.merra_sam_enso(m_lat, m_lon, p_mb500, n_mb500, m_mb500)
    


# if __name__ == '__main__':
#     m_mast_arr = merra_analize.build_ds('test_m')
#     t_mast_arr = merra_analize.build_ds('test')
#     lat = m_mast_arr[0,0]
#     lon = m_mast_arr[0,1]
#     m_mb500, m_slp, m_t2m = merra_analize.get_mean_vals(m_mast_arr)
#     t_mb500, t_slp, t_t2m = merra_analize.get_mean_vals(t_mast_arr)

#     p = plt.contourf(lon, lat, (m_mb500 - t_mb500), cmap ='jet')
#     c = plt.colorbar(p)
#     plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/tester.png')
