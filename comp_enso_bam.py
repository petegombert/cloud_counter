"""
*This function is designed to compare the correlation between ENSO and BAM. The reasoning behind this
is because there 
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
from scipy.optimize import curve_fit
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
import new_make_raw_bam_anim
import make_raw_bam_anim
import bam_analysis

months_arr = np.arange(1,13)
months = go_stats.make_np_list(months_arr)
months = go_stats.fix_days(months)


if __name__ == '__main__':
    enso_index = indexes.enso_full(months, np.arange(2006,2009), 0)
    eof1, bam_index = bam_analysis.run_main_bam()
    enso_index = enso_index[:-1]

    corr_mtx = np.corrcoef(enso_index, bam_index)
    corr = corr_mtx[0,1]
    print(corr)

    plt.plot(np.arange(0,enso_index.shape[0]), enso_index, color='Red', label='ENSO')
    plt.plot(np.arange(0,bam_index.shape[0]), bam_index, color='Blue', label='BAM')
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testinger.png')
      