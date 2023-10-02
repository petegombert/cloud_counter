"""
* This code explores the teleconnections between ENSO and the SAM. Specifically the differences
in cloud occurrence between positive relationship years (2006-2008), and negative relationship years
(2010-2012). Literature reference (not from code but theory): https://doi.org/10.5194/cp-16-743-2020

Leave notes:

To do:
    *Data only starts at 062006, should ask Jay if I should limit the correlation coefficient.

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
import get_data
import mast_plot
import go_stats
import build_ann_anom
import coarsify
import diveinnewdata
import jja_global_avg
import overview
import atmos_6040_final
import indexes
import overview

#pos_rel_years = ['2015','2016']
#neg_rel_years = ['2010','2012']
pos_rel_years = ['2013','2014']
neg_rel_years = ['2008','2009']
early_years = ['2006','2007','2008','2009']
middle_years = ['2010','2012','2013']
late_years = ['2014','2015','2016']

mean = ['2006','2007','2008','2009','2010','2012','2013','2014','2015','2016']

#months = ['12','01','02']
months = 'all'
#months = ['01','02','03','04','05','06','07','08','09','10','11','12']

# if __name__ == '__main__':
#     lat, lon, cld_mean = overview.call_main(mean, months, 0)
#     lat, lon, cld_pos = overview.call_main(pos_rel_years, months, 0)
#     lat, lon, cld_neg = overview.call_main(neg_rel_years, months, 0)

#     mast_plot.sam_enso_rel(cld_pos, cld_neg, cld_mean)

# *For the timeperiod plots (split into early, mid, late anomalies)
# if __name__ == '__main__':
#     lat, lon, cld_mean = overview.call_main(mean, months, 0, 'yes', [60,-90,-180,180])
#     lat, lon, cld_ear = overview.call_main(early_years, months, 0, 'yes', [60,-90,-180,180])
#     lat, lon, cld_mid = overview.call_main(middle_years, months, 0, 'yes', [60,-90,-180,180])
#     lat, lon, cld_lat = overview.call_main(late_years, months, 0, 'yes', [60,-90,-180,180])

#     mast_plot.sam_enso_timep(cld_mean, cld_ear, cld_mid, cld_lat)

# *For the timeperiod plots, just in orthoginal view focused on the southern ocean    
# if __name__ == '__main__':
#     lat, lon, cld_mean = overview.call_main(mean, months, 0, 'yes', [-30,-90,-179,179])
#     lat, lon, cld_ear = overview.call_main(early_years, months, 0, 'yes', [-30,-90,-179,179])
#     lat, lon, cld_mid = overview.call_main(middle_years, months, 0, 'yes', [-30,-90,-179,179])
#     lat, lon, cld_lat = overview.call_main(late_years, months, 0, 'yes', [-30,-90,-179,179])

#     print(cld_mean.shape)

#     mast_plot.sam_enso_timep_so(lat, lon, cld_mean, cld_ear, cld_mid, cld_lat)



