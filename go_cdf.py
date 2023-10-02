#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:19:51 2022

@author: petergombert

* This script builds a ncdf from data provided
leave notes: getting data from the import run. Need to figure out how to return the data back here.

To do:
	Need to figure out how to aggrigate multiple months. Should just be able to add them up, not that clean though...
"""

import os
import numpy as np
import pandas as pd
import netCDF4 as cdf
import julian
from datetime import datetime
import sys
import subprocess
from warnings import filterwarnings
import get_data
import mast_plot
import go_stats

filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

"""
must provide build_cdf func with the following:
1: number of dimensions
2: number of variables
3: np.array of dimension names
4: np.array of dimension data
5: np.array of dimensions that a variable uses
5: np.array of variable names
6: np.array of data
7: path for ncdf
8: cdf name
"""
def build_cdf(num_d, num_v, d_name, d_data, v_dim, v_name, v_data, path, cdf_name):
    os.chdir(path)
    ncfile = cdf.Dataset(cdf_name + '.cdf', mode='w', format='NETCDF4')
    for d in d_name:
        d = ncfile.createDimension(d)
    vstep = np.arange(0, num_v.shape[0])
    for v in vstep:
        v = ncfile.createVariable(v[v], np.float32, (v_dim[v],)) #this might need to be reworked for variables that have multiple dimensions
    for v in num_v:

