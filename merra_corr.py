"""
* This code initiates the merra analysis for the sam enso correlation analysis. 

To do:

Leave notes:
    Need to figure out the mean dataset, since there's no monthly data in there
    Just need to append the data for each year.
    Not sure how to do the DJF month select...
    Annoying because it is the most important one.

    * We've decided to abandon this analysis for right now. I was trying to build the 3 month averaged comparison but 
    got a little stuck and the next day Jay told me he doesn't think this is worth our time anymore.
    Maybe should keep investigating on my own though. I think that the data isn't super indicative of anything right 
    now though. 

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
import merra_mean
import indexes
import go_stats
import mast_plot
import overview

earliest_year = ['1980','1981','1982','1983','1984','1985','1986','1987','1988','1989']
early_year = ['1990','1991','1992','1993','1994','1995','1996','1997','1998','1999']
mid_year = ['2000','2001','2002','2003','2004','2005','2006','2006','2008','2009']
late_year = ['2010','2012','2013','2014','2015','2016','2017','2018','2019']

months = ['01','02','03','04','05','06','07','08','09','10','11','12']

#!!Need to make sure path is right in all merra scripts before running!!

def call_main1(years, months, cat):
    for year in years:
        for month in months:
            merra_first.merra_call_get(year, month, cat, [90,-90,-180,180])
        mast_arr = merra_analize.build_ds(year, cat)
        lat = mast_arr[0,0]
        lon = mast_arr[0,1]
        m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mast_arr)
        merra_mean.save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, year, cat)
        merra_mean.save_year(lat, lon, m_mb500, m_slp, m_t2m, m_pbltop, m_q500, year, 'mean')

#* This code builds the merra dataset
# if __name__ == '__main__':
#     ana_direct = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/SAM_GLB'
#     os.chdir(ana_direct)
#     t =4
#     for filename in os.listdir(ana_direct):
#         if filename.split('_')[0][-1] == str(t):
#             print('starting: ' + filename)
#             year_start = filename.split('_')[1]
#             year_end = filename.split('_')[2]
#             #*This removes 2011 from the analysis.
#             if year_end == '2012':
#                 year_end = '2011'
#             years = np.arange(int(year_start), int(year_end))
#             years = years.tolist()
#             years_str = map(str, years)
#             years = list(years_str)
#             call_main(years, months, filename)
#             t = t+1
#     exit()

#* This code analizes the dataset & plots results. 
# if __name__ == '__main__':
#     ana_direct = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/SAM_GLB'
#     os.chdir(ana_direct)
    # p1_mast = merra_analize.build_ds('all', 'p1_1980_1984')
    # p2_mast = merra_analize.build_ds('all', 'p2_1984_1993')
    # p3_mast = merra_analize.build_ds('all', 'p3_1993_1997')
    # p4_mast = merra_analize.build_ds('all', 'p4_1997_2012')
    # p5_mast = merra_analize.build_ds('all', 'p5_2012_2017')
    # p6_mast = merra_analize.build_ds('all', 'p6_2017_2020')
    # mean_mast = merra_analize.build_ds('all', 'mean')

    # slp_ts = mean_mast[:,5]
    # tstep = np.arange(0, slp_ts.shape[0])
    # for t in tstep:
    #     mean_hold = np.mean(slp_ts[t])
    #     if t == 0:
    #         mean_arr = np.array([mean_hold])
    #         continue
    #     mean_arr = np.append(mean_arr, mean_hold)

    # plt.plot(np.arange(1980,2019), mean_arr)
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')

    # exit()
    # lat = p1_mast[0,0]
    # lon = p1_mast[0,1]

    # p1_500mb, p1_slp, p1_t2m, p1_pbltop, p1_q500 = merra_analize.get_mean_vals(p1_mast)
    # p2_500mb, p2_slp, p2_t2m, p2_pbltop, p2_q500 = merra_analize.get_mean_vals(p2_mast)
    # p3_500mb, p3_slp, p3_t2m, p3_pbltop, p3_q500 = merra_analize.get_mean_vals(p3_mast)
    # p4_500mb, p4_slp, p4_t2m, p4_pbltop, p4_q500 = merra_analize.get_mean_vals(p4_mast)
    # p5_500mb, p5_slp, p5_t2m, p5_pbltop, p5_q500 = merra_analize.get_mean_vals(p5_mast)
    # p6_500mb, p6_slp, p6_t2m, p6_pbltop, p6_q500 = merra_analize.get_mean_vals(p6_mast)
    # m_500mb, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mean_mast)

    # mast_plot.merra_corr(lat, lon, m_slp, p1_slp, p2_slp, p3_slp, p4_slp, p5_slp, p6_slp)

#* This code makes just positive and negative correlation years. Maybe removes some time stuff
# if __name__ == '__main__':

#     p1_mast = merra_analize.build_ds('all', 'p1_1980_1984')
#     p3_mast = merra_analize.build_ds('all', 'p3_1993_1997')
#     p5_mast = merra_analize.build_ds('all', 'p5_2012_2017')
#     pos_mast = np.append(p1_mast, p3_mast, axis=0)
#     pos_mast = np.append(pos_mast, p5_mast, axis=0)

#     p2_mast = merra_analize.build_ds('all', 'p2_1984_1993')
#     p4_mast = merra_analize.build_ds('all', 'p4_1997_2012')
#     p6_mast = merra_analize.build_ds('all', 'p6_2017_2020')
#     neg_mast = np.append(p2_mast, p4_mast, axis=0)
#     neg_mast = np.append(neg_mast, p6_mast, axis=0)

#     mean_mast = merra_analize.build_ds('all', 'mean')

#     lat = mean_mast[0,0]
#     lon = mean_mast[0,1]

#     p_500mb, p_slp, p_t2m, p_pbltop, p_q500 = merra_analize.get_mean_vals(pos_mast)
#     n_500mb, n_slp, n_t2m, n_pbltop, n_q500 = merra_analize.get_mean_vals(neg_mast)
#     m_500mb, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mean_mast)

#     mast_plot.merra_sam_enso(lat, lon, p_slp, n_slp, m_slp)

if __name__ == '__main__':
    jja_p1 = merra_analize.build_ds(['1980','1981','1982','1983'], 'p1_1980_1984', ['06','07','08'])
    jja_p2 = merra_analize.build_ds(['1984','1985','1986','1987','1988','1989','1990','1991','1992'], 'p2_1984_1993', ['06','07','08'])
    jja_p3 = merra_analize.build_ds(['1993','1994','1995','1996'], 'p3_1993_1997', ['06','07','08'])
    jja_p4 = merra_analize.build_ds(['1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010'], 'p4_1997_2012', ['06','07','08'])
    jja_p5 = merra_analize.build_ds(['2012','2013','2014','2015','2016'], 'p5_2012_2017', ['06','07','08'])
    jja_p6 = merra_analize.build_ds(['2017','2018','2019'], 'p6_2017_2020', ['06','07','08'])


    p1_500mb, p1_slp, p1_t2m, p1_pbltop, p1_q500 = merra_analize.get_mean_vals(jja_p1)
    p2_500mb, p2_slp, p2_t2m, p2_pbltop, p2_q500 = merra_analize.get_mean_vals(jja_p2)
    p3_500mb, p3_slp, p3_t2m, p3_pbltop, p3_q500 = merra_analize.get_mean_vals(jja_p3)
    p4_500mb, p4_slp, p4_t2m, p4_pbltop, p4_q500 = merra_analize.get_mean_vals(jja_p4)
    p5_500mb, p5_slp, p5_t2m, p5_pbltop, p5_q500 = merra_analize.get_mean_vals(jja_p5)
    p6_500mb, p6_slp, p6_t2m, p6_pbltop, p6_q500 = merra_analize.get_mean_vals(jja_p6)
  #  m_500mb, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(mean_mast)
        

