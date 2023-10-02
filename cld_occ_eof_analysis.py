"""
*This script will be used for the EOF analysis of the cloud occurrence statistics.

To do:

Leave notes:
Here's where I am. Take the first PC and multiply each value in the EOF by the PC, then do the correlation.
The PC is the time series of changes and the EOF is the shape. Then you correlate it with the data.

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

def proj_eof(eof1, pc1):
    ystep = np.arange(0,eof1.shape[0])
    xstep = np.arange(0,eof1.shape[1])
    for y in ystep:
        for x in xstep:
            eof1_p = eof1[y,x]*pc1

            if x == 0:
                eof1_p_xarr = np.array([eof1_p])
                continue
            eof1_p_xarr = np.append(eof1_p_xarr, np.array([eof1_p]), axis=0)
        if y == 0:
            eof1_p_arr = np.array([eof1_p_xarr])
            continue
        eof1_p_arr = np.append(eof1_p_arr, np.array([eof1_p_xarr]), axis=0)

    return eof1_p_arr

years = go_stats.make_np_list(np.arange(2006,2017))
twen011 = years.pop(5)

months = ['12','01','02']

if __name__ == '__main__':
    cld_occ, lat, lon = oned_corr.eof_call_main(years, months, 0)
    m_cld_occ = np.mean(cld_occ, axis=2)
    #Subtract the mean from each step
    tstep = np.arange(cld_occ.shape[-1])
    for t in tstep:
        cld_occ[:,:,t] = cld_occ[:,:,t]-m_cld_occ
    cld_occ = np.flipud(cld_occ)
    print(cld_occ.shape)
    test = np.transpose(cld_occ)

    eof1, pcs = eof_analysis.go_eof_cld(test)
    pc1 = pcs[:,0]
    pc1 = sp.zscore(pc1, ddof=1)

    eof1 = np.rot90(eof1, k=-1)
    eof1 = np.fliplr(eof1)
    print(eof1.shape)

    proj_arr = proj_eof(eof1, pc1)

    corr_arr = eof_analysis.build_cov_arr_1st_dem(cld_occ, proj_arr, 'reg', pc1)

    norm = mast_plot.norm_maker(corr_arr)
    p = plt.imshow(corr_arr, cmap='RdBu', norm=norm)
    plt.colorbar(p)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')


    