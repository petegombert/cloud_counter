"""
* This code will gather the merra data for high and low ENSO years.

To do:

Leave notes:
    `Need to make sure it's being saved as 'keep'

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
import indexes
import merra_analize
import merra_mean

years = ['2006','2007','2008','2009','2010','2012','2013','2014','2015','2016']
#years = np.arange(1980,2023)
#months = ['06','07','08']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
#indices = indexes.enso(['07','08'])
indices = indexes.enso_full(['07','08'])
sam_pos_years = ['2015','2016']
sam_neg_years = ['2010','2012']

# step = 0
# #Before running this - merra_analize.build_ds() needs to be ready!!
# for index in indices:
#     merra_mean.import_year(str(years[step]))
#     if float(index) >= 0.35:
#         print('high: ' + str(years[step]))
#         merra_first.merra_call_get(str(years[step]), '06', 'high', [-40,-90,-180,180])
#         merra_first.merra_call_get(str(years[step]), '07', 'high', [-40,-90,-180,180])
#         merra_first.merra_call_get(str(years[step]), '08', 'high', [-40,-90,-180,180])
#         mast_arr = merra_analize.build_ds('high')
#         lat = mast_arr[0,0]
#         lon = mast_arr[0,1]
#         m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mast_arr)
#         merra_mean.save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, str(years[step]), 'high')

#     if float(index) <= -0.35:
#         print('low: '+ str(years[step]))
#         merra_first.merra_call_get(str(years[step]), '06', 'low', [-40,-90,-180,180])
#         merra_first.merra_call_get(str(years[step]), '07', 'low', [-40,-90,-180,180])
#         merra_first.merra_call_get(str(years[step]), '08', 'low', [-40,-90,-180,180])
#         mast_arr = merra_analize.build_ds('low')
#         lat = mast_arr[0,0]
#         lon = mast_arr[0,1]
#         m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mast_arr)
#         merra_mean.save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, str(years[step]), 'low')

#     step = step+1

#*This is the SAM correlation with ENSO section 
for year in years:
    for month in months:
        merra_first.merra_call_get(year, month, 'mean', [90,-90,-180,180])
    mast_arr = merra_analize.build_ds('mean')
    lat = mast_arr[0,0]
    lon = mast_arr[0,1]
    m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mast_arr)
    merra_mean.save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, str(year), 'mean') 

exit()

for year in sam_pos_years:
    for month in months:
        merra_first.merra_call_get(year, month, 'pos_years', [90,-90,-180,180])
    mast_arr = merra_analize.build_ds('pos_years')
    lat = mast_arr[0,0]
    lon = mast_arr[0,1]
    m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mast_arr)
    merra_mean.save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, str(year), 'pos_years')

for year in sam_neg_years:
    for month in months:
        merra_first.merra_call_get(year, month, 'neg_years', [90,-90,-180,180])
    mast_arr = merra_analize.build_ds('neg_years')
    lat = mast_arr[0,0]
    lon = mast_arr[0,1]
    m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mast_arr)
    merra_mean.save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, str(year), 'neg_years')
