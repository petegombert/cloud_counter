# cloud_counter
Graduate Research project

This repository was created to organize my cloud occurrence master's thesis project. An additional benefit is that it will make it easier for future students to use my research and code to continue investigating this topic.

Sections:
  1. Raw Data Analysis (zonal&global)
  2. Initial Cloud Occurrence Calculations

****

1. RAW DATA ANALYSIS
1.1 Data access
These scripts are designed to efficiently get the data required. They are designed to be dynamic, where you input the 'dataid' and a selection. This selection specifies what the dataid looks like. The options are as follows:
    Yearly flag = 1: Get all months of a specific year
    Yearly flag = 0: Get all years of a specific month
    Yearly flag = 2: Get specific month & specific year

This selection is the same for both datasets. The scripts that get data from each dataset are:
    Zonal dataset: get_data.py
    Global/thickness dataset: diveinnewdata.py

The function that initiates these codes is called find_files().

The output of this function will be an array with multiple elements/dimensions. For each time step given the input 'dataid', there will be an array of: 
    get_data.py: (year/month, lat, lon, z, dbz, cld_occ)
    diveinnewdata.py: (lat, lon, cbtz_bin, cbtz_dbin, dz_bin, dz_dbin, single_layer_occ, tot_count, year)

1.2 Interpretation/Visualization
After gathering the correct data, visualizing the raw data is important. However, some manipulations are required to get the data in 
the correct units. The data comes in a multidimensional numpy array, for each dataset the 'cloud occurrence' values are technically cloud counts. They must be divided by the total number of cloud counts (overpasses) to return the cloud occurrence. Also, we are looking at specific clouds which is where the dbz, cbtz_bin, and dz_bin values come from. These are bins of different cloud type. Dbz is radar reflectivity & should be pretty self-explanatory. cbtz_bin & dz_bin are the middle values of the cloud base/top & thickness bin values. The ctbz_dbin & dz_dbin values are the shape of each bin (ie distance from the middle value to the next bin.)

To functionally do these operations the script go_stats.py is used. This script houses most statistical operations/array manipulations that are required throughout the repository. In this case, the functions that are called to preform the data manipulation are:
    Zonal dataset: go_stats.build_stats()
    Global/thickness dataset: jja_global_avg.go_stats_func()
**These functions will return the mean cloud occurrence value of the input array. Thus, if the input array consists of months (1, 2, 3) are being looped by years, the result will be the mean values of Jan., Feb., and March for each year. If the individual values are desired, call the functions:
    Zonal dataset: go_stats.build_stats_wm()
    Global/thickness dataset: jja_global_avg.go_stats_func_wm()

The global/thickness dataset has a script that is well equipt & function for building any subset of raw data that one could want.
This script is called overview.py, 

