"""
*This code looks into the relationship between BAM and NAO. After looking at the cloud occurrence 
we wondered whether the BAM and NAO were correlated. This code looks into this question.
    To do:

    Leave notes:

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

years = np.arange(2006,2010)
months_arr = np.arange(1,13)
months = go_stats.make_np_list(months_arr)
months = go_stats.fix_days(months)

if __name__ == '__main__':
    raw_nao_indx = indexes.nao(months, years, 0)
    eof1, raw_bam_indx = bam_analysis.run_main_bam()

    daily_vals_nao = go_stats.make_daily_mean_vals(raw_nao_indx, 32, 3)
    nao_indx = go_stats.remove_mean_bam(raw_nao_indx, daily_vals_nao)

    daily_vals_bam = go_stats.make_daily_mean_vals(raw_bam_indx, 4, 3)
    daily_bam_indx = go_stats.remove_mean_bam(raw_bam_indx, daily_vals_bam)
    bam_indx = go_stats.make_monthly_mean(daily_bam_indx, 4, 4, 12)

    corrtest = np.corrcoef(nao_indx, bam_indx)
    test_corr = corrtest[0,1]

    model, intercept, pval, rval, stderr = sp.linregress(nao_indx, bam_indx)
    print(model)
    print(test_corr)

    mast_plot.test_bam_nao_rel(bam_indx, nao_indx, model, intercept, pval, test_corr)
   


