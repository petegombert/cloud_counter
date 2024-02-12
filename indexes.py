"""
This script is where the general circulation indeces are imported and handled.

Leave notes:
    I'm making a full merra/enso analysis so i need the indices for all merra years.
    Want to see if enso patterns are persistant. 
    Also going to need to detrend the data...
    Trying to see if leading pc of U data matches up with the SAM index. But im not sure the 
    sam index code is really working that well. Need to make sure it's getting the right values.
    
    
To do:
    Truncate the index values. Since they're imported they're longer than when I just typed them.
    Need to make a year shifter - I did, it just needs to be tested.
    Need to fix the year shifter, it doesn't work right. The SAM function works right.
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
import datetime

def enso(months, years, mix_years):
    #** months must be passed as list.
    #data = np.genfromtxt('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/indices/enso_index.txt', delimiter='    ')
    data = pd.read_csv('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/indices/enso_index.txt', sep='     ', header=None, engine='python')
    data = data.to_numpy()
    hit2011 = np.argwhere(data[:,0] == 2011)
    data = np.delete(data, hit2011, 0)
    year_min = years[0]
    year_max = years[-1]+2
    inx_min = np.argwhere(year_min == data[:,0])
    inx_max = np.argwhere(year_max == data[:,0])
    
    tstep = np.arange(inx_min, inx_max)
    for t in tstep:
        #print(data[t,0])
        for month in months:
            if mix_years == 1 and month == '01':
                t_send = t+1
            else:
                t_send = t
            print(data[t_send, 0])
            month_int = int(month)
            index_hold = data[t_send,month_int]
            if month == months[0]:
                index_arr = np.array([index_hold])
                continue
            index_arr = np.append(index_arr, np.array([index_hold]), axis=0)
        index_send = np.mean(index_arr)
        if t == inx_min:
            mast_arr = np.array([index_send])
            continue
        mast_arr = np.append(mast_arr, np.array([index_send]), axis=0)
        #**This code gets values without mean.
        # index_send = index_arr
        # if t == inx_min:
        #     mast_arr = index_arr
        #     continue
        # mast_arr = np.append(mast_arr, index_send)

    return mast_arr

def sam(months, years, mix_years):
    #Years must be in order
    data = pd.read_csv('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/indices/sam_index_data.txt', sep='  ', header=None, engine='python')
    data = data.to_numpy()
    hit2011 = np.argwhere(data[:,0] == 2011)
    data = np.delete(data, hit2011, 0)
    year_min = years[0]
    year_max = years[-1]+2
    inx_min = np.argwhere(year_min == data[:,0])
    inx_max = np.argwhere(year_max == data[:,0])

    tstep = np.arange(inx_min, inx_max)
    for t in tstep:
        for month in months:
            if mix_years == 1 and month == '01':
                t_send = t+1
            elif month == months[0]:
                t_send = t
            month_int = int(month)
            index_hold = data[t_send,month_int]
            if month == months[0]:
                index_arr = np.array([index_hold])
                continue
            index_arr = np.append(index_arr, np.array([index_hold]), axis=0)

        index_send = np.mean(index_arr)
        if t == inx_min:
            mast_arr = np.array([index_send])
            continue
        mast_arr = np.append(mast_arr, np.array([index_send]), axis=0)
        #**This code gets values without mean.
        # index_send = index_arr
        # if t == inx_min:
        #     mast_arr = index_arr
        #     continue
        # mast_arr = np.append(mast_arr, index_send)

    return mast_arr    

def enso_full(months, years, mix_years):
    #** months must be passed as list.
    data = pd.read_csv('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/indices/enso_full.txt', sep='     ', header=None, engine='python')
    data = data.to_numpy()
    data = data[:-1,:]
    hit2011 = np.argwhere(data[:,0] == 2011)
    data = np.delete(data, hit2011, 0)
    year_min = years[0]
    if years.shape[0] == 1:
        year_max = year_min+1
    else:
        year_max = years[-1]+2
    print(year_max)
    inx_min = np.argwhere(year_min == data[:,0])
    inx_max = np.argwhere(year_max == data[:,0])
    
    tstep = np.arange(inx_min, inx_max)
    print(tstep)
    for t in tstep:
        print(data[t,0])
        for month in months:
            if mix_years == 1 and month == '01':
                t_send = t+1
            else:
                t_send = t
            print(data[t_send,0])
            month_int = int(month)
            index_hold = data[t_send,month_int]
            index_hold = float(index_hold)
            if month == months[0]:
                index_arr = np.array([index_hold])
                continue
            index_arr = np.append(index_arr, np.array([index_hold]), axis=0)
        # index_send = np.mean(index_arr)
        # if t == inx_min:
        #     mast_arr = np.array([index_send])
        #     continue
        # mast_arr = np.append(mast_arr, np.array([index_send]), axis=0)
        #**This code gets values without mean.
        index_send = index_arr
        if t == inx_min:
            mast_arr = index_arr
            continue
        mast_arr = np.append(mast_arr, index_send)

    return mast_arr

def nao(months, years, mix_years):
    """
    This function gets the nao index for the months & years selected. Mix years will skip to the next year when ==1 & month ==1
    Notes:
        Not sure why but the indexing removes the first year (1950)
    """
    data = pd.read_csv('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/indices/nao_data.txt', sep='   ', header=1, engine='python')
    data = data.to_numpy()
    #data = data[1:]
    hit2011 = np.argwhere(data[:,0] == 2011)
    data = np.delete(data, hit2011, 0)
    year_min = years[0]
    if years.shape[0] == 1:
        year_max = year_min+1
    else:
        year_max = years[-1]+1
    print(year_max)
    inx_min = np.argwhere(year_min == data[:,0])
    inx_max = np.argwhere(year_max == data[:,0])
    
    tstep = np.arange(inx_min, inx_max)
    print(tstep)
    for t in tstep:
        for month in months:
            if mix_years == 1 and month == '01':
                t_send = t+1
            else:
                t_send = t
            print(data[t_send,0])
            month_int = int(month)
            index_hold = data[t_send,month_int]
            index_hold = float(index_hold)
            if month == months[0]:
                index_arr = np.array([index_hold])
                continue
            index_arr = np.append(index_arr, np.array([index_hold]), axis=0)
        # index_send = np.mean(index_arr)
        # if t == inx_min:
        #     mast_arr = np.array([index_send])
        #     continue
        # mast_arr = np.append(mast_arr, np.array([index_send]), axis=0)
        #**This code gets values without mean.
        index_send = index_arr
        if t == inx_min:
            mast_arr = index_arr
            continue
        mast_arr = np.append(mast_arr, index_send)

    return mast_arr

def sort_index(index, data):
    #Sorts the sent index and the data in the same way for making the box&whisker plots.
    def apply_to_data():
        print('DATA:')
        print(data.shape)
        sstep = np.arange(0, data_arr.shape[0])
        for s in sstep:
            t = int(data_arr[s])
            print(t)
            print(type(t))
            try:
                data_arr_f = np.append(data_arr_f, np.array([data[:,:,t]]),axis=0)
            except UnboundLocalError:
                data_arr_f = np.array([data[:,:,t]])
        print(data_arr_f.shape)
        return(data_arr_f)
    
    tstep = np.arange(0, index.shape[0])
    print(index.shape[0])
    index_arr = np.array([])
    data_arr = np.array([[]])
    for t in tstep:
        if t == 0:
            index_val = index[t]
            index_arr = np.insert(index_arr, 0, index_val)
            #data_arr = np.array([data[:,:,t].flatten()])
            print(data_arr.shape)
            data_arr = np.insert(data_arr, 0, t)
            continue
        if t == 1:
            index_test = index[t]
            if index_test >= index_val:
                index_arr = np.insert(index_arr, 1, index_test)
                #data_arr = np.insert(data_arr, 1, data[:,:,t])
                data_arr = np.insert(data_arr, 1, t)
            else:
                index_arr = np.insert(index_arr, 0, index_test)
                #data_arr = np.insert(data_arr, 0, np.array([data[:,:,t].flatten()]))
                data_arr = np.insert(data_arr, 0, t)
                print(data_arr.shape)
            continue
        index_test = index[t]
        diff = abs(index_arr - index_test)
        hit = np.argwhere(diff == np.min(diff))
        factor = hit.shape[0]
        if index_arr[int(hit[0])] < index_test:
            index_arr = np.insert(index_arr, int(hit[0])+factor, index_test)
            #data_arr = np.insert(data_arr, int(hit[0])+factor, data[:,:,t])
            data_arr = np.insert(data_arr, int(hit[0])+factor, t)
        else:
            index_arr = np.insert(index_arr, int(hit[0]), index_test)
            #data_arr = np.insert(data_arr, int(hit[0]), data[:,:,t])
            data_arr = np.insert(data_arr, int(hit[0]), t)
    
    print(data_arr)
    data_arr = apply_to_data()
    #print(data_arr.shape)

    return index_arr, data_arr   

def enso_sam_corr():
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    enso_ind = enso_full(months, np.arange(1980,2020), 0)
    #enso_ind = enso(months, np.arange(2006, 2016), 0)
    sam_ind = sam(months, np.arange(1980, 2020), 0)
    corr_hold = np.corrcoef(enso_ind, sam_ind)
    corr = corr_hold[0,1]
    corr_arr = go_stats.make_rolling_corr(enso_ind, sam_ind)
    years = np.arange(1980, 2021)
    years = np.delete(years, -10)
#* if there is a 0.08333333333333333333, it is getting monthly values and it must be changed to monthly mode in the sam and enso functions
    fig, ax = plt.subplots(2,1, figsize=(10,10))
    ax1 = ax[0]
    ax1.plot(years, enso_ind, color='blue', label='enso')
    ax1.plot(years, sam_ind, color='red', label='sam')
    ax1.axvline(x=1980, color='black')
    ax1.axvline(x=1984, color='black')
    ax1.axvline(x=1993, color='black')
    ax1.axvline(x=1997, color='black')
    ax1.axvline(x=2012, color='black')
    ax1.axvline(x=2017, color='black')
    ax1.axvline(x=2020, color='black')
    ax1.set_title('correlation coefficient: '+ str(corr))
    ax1.legend()

    ax2 = ax[1]
    ax2.plot(np.arange(1980, 2020),  corr_arr)
    ax2.set_title('Rolling Correlation Coefficient')
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/data/index_plts/enso_sam_1980_2020_ann_corr.png')

def season_histo(months, osc):
    if osc == 'sam':
        vals = sam(months, np.arange(1980,2020), 0)
    else:
        vals = enso_full(months, np.arange(1980,2020), 0)
    
    mast_plot.index_histo(vals, osc, 'SON')

def make_dt_arr(day_arr_s, num_years, num_months, starting_date, years_skip):
    """
    This function is to deal with the time averaging in a daily index.
    It will be called by a function in go_stats.py and will ultimately help return the
    BAM index as monthly mean values rather than every 4 days as it is now.
    How it works is it takes the dates each month, the number of years, and number of months
    and it loops through the years and months, and for each month it loops through the days
    to get an array of datetime values. These will be used in go_stats.make_monthly_mean()
    *Right now years_skip only works if it is 1 year. Should make it work for any year.

    Starting day should be string of shape (m-d-y)
    """
    starting_day = datetime.datetime.strptime(starting_date,'%m-%d-%Y')
    print('Starting day: '+ datetime.datetime.strftime(starting_day,'%m-%d-%Y'))
    ystep = np.arange(0,num_years)
    mstep = np.arange(0,num_months)
    dstep = np.arange(0,day_arr_s.shape[0])
    for y in ystep:
        year = starting_day.year+y
        if int(year) == years_skip:
            continue
        print(year)
        for m in mstep:
            month = starting_day.month+m
            for d in dstep:
                try:
                    date_hold = datetime.date(year,month,day_arr_s[d])
                except ValueError:
                    continue

                try:
                    date_arr = np.append(date_arr, date_hold)
                except UnboundLocalError:
                    date_arr = np.array([date_hold])
    
    return date_arr
    
if __name__=='__main__':
    test = nao(['11'],np.array([1998,1999]),0)
    print(test)
    exit()

if __name__ == '__main__':
    make_dt_arr(np.arange(1,32,4),4,12)
    exit()

if __name__ == '__main__':
    sam_years = np.arange(1985,1986)
    months = ['12','01','02']
    index = sam(months, sam_years, 1)
    plt.plot(index)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg/sam_index.png')
    exit()

    # season_histo(['09','10','11'], 'sam')
    # exit()

    #enso(['06','07'], np.arange(2006,2016), 0)
#* This section of code generates averages for years in the sam and enso indexes
    # test1 = enso(['01','02','03','04','05','06','07','08','09','10','11','12'], np.arange(2008,2009), 0)
    # print(np.mean(test1))
    # test1sam = sam(['01','02','03','04','05','06','07','08','09','10','11','12'], np.arange(2008,2009), 0)
    # print(np.mean(test1sam))
    # test2 = enso(['01','02','03','04','05','06','07','08','09','10','11','12'], np.arange(2013,2014), 0)
    # print(np.mean(test2))
    # test2sam = sam(['01','02','03','04','05','06','07','08','09','10','11','12'], np.arange(2013,2014), 0)
    # print(np.mean(test2sam))
    # exit()

    enso_sam_corr()
    exit()
    #index = enso_full(['07','08'])
    #sort_index(index, 0)
    sam_indices_win = sam(['06'], np.arange(2006,2016))
    sam_indices_sum = sam(['12'], np.arange(2006,2016))

    plt.plot(np.arange(2006,2016), sam_indices_win, color='blue', label='winter')
    plt.plot(np.arange(2006,2016), sam_indices_sum, color='orange', label='summa')
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/testing.png')
    #print(sam_indices)
