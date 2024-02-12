"""
Likely the first contact with the repository. This code currently (2/2024) only will gather the raw occurrence data and plot it.
If you wish to do something with it, import it like any other package (script must be in the same directory as other scripts it imports and the 
ones that import it) and call 'call_main.give_me_data1() or call_main.give_me_data2()'
"""
import oned_corr
import go_stats
import os
import mast_plot
import numpy as np

#*****This is where the dataset is selected. id=1 is the zonal dataset & id=2 is the global cloud occurrence dataset.

def give_me_data1(years, months, dbz_bot, dbz_top):
    """
  This function is the callable version of this script for dataset 1.
    """

def give_me_data2(years, months, z_bot, z_top, dz_bot, dz_top):
    """
  This function is the callable version of this script for dataset 2.
    """

if __name__ == '__main__':
    cwd = os.getcwd()
    dataset_id=1
    z_bot = 0
    z_top = 1000
    dz_bot = 0
    dz_top = 1500
    
    dbz_bot = 0
    dbz_top = 22
    
    years_arr = np.arange(2006,2016)
    years = go_stats.make_np_list(years_arr)
    months_arr = np.arange(1,13)
    months = go_stats.make_np_list(months_arr)
    months = go_stats.fix_days(months)
    
    latlon_bounds = [45,15,-120,-60]
    restrict_domain_call = 0
    if dataset_id == 1:
        print('Dataset 1')

    elif dataset_id == 2:
        print('Dataset 2')
        cld_arr, lat, lon = oned_corr.eof_call_main_wrestrict(years, months, restrict_domain_call, 0, latlon_bounds, z_bot, z_top, dz_bot, dz_top)

        mast_plot.overview_plt(lat, lon, cld_arr, dz_bot, dz_top, z_bot, z_top, time, cwd)
        
  
