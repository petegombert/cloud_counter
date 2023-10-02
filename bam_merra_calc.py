"""
*This script will gather the data for the BAM calculation. We will get merra data, it will likely 
look similar to some of the merra analysis. It will also need to be combined with EOF analysis. 

To do:

Leave notes:

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
from time import sleep

merra_data_path = '/uufs/chpc.utah.edu/common/home/u1113223/merra_hold'

def download_file(month, year, day):
    #print('starting '+month+'/'+year)
    os.chdir('/uufs/chpc.utah.edu/common/home/u1113223/merra_hold')
    if int(year) < 1992:
        URL = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I6NPANA.5.12.4/'+ str(year) + '/' + str(month) + '/MERRA2_100.inst6_3d_ana_Np.'+str(year)+str(month)+str(day)+'.nc4'
    elif int(year) < 2001:
        URL = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I6NPANA.5.12.4/'+ str(year) + '/' + str(month) + '/MERRA2_200.inst6_3d_ana_Np.'+str(year)+str(month)+str(day)+'.nc4'
    elif int(year) < 2011:
        URL = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I6NPANA.5.12.4/'+ str(year) + '/' + str(month) + '/MERRA2_300.inst6_3d_ana_Np.'+str(year)+str(month)+str(day)+'.nc4'
    else:
        URL = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I6NPANA.5.12.4/'+ str(year) + '/' + str(month) + '/MERRA2_400.inst6_3d_ana_Np.'+str(year)+str(month)+str(day)+'.nc4'
    print(URL)

    username = 'petegombert'
    password = 'peteG@115'

    FILENAME = str(year)+str(month)+'.nc4'


    with requests.Session() as session:
        session.auth = (username, password)
        r1 = session.request('get', URL)

        r = session.get(r1.url, auth=(username, password))

    #result = requests.get(URL, auth=())
    try:
        r.raise_for_status()
        f = open(FILENAME,'wb')
        f.write(r.content)
        f.close()
        print('contents of URL written to '+FILENAME)
        return 'data'
    
    except requests.exceptions.HTTPError:
        print('requests.get() returned an error code '+str(r.status_code))
        return 'no data'
    except requests.exceptions.ConnectionError:
        #* This handles the time error by sleeping. 
        #Not sure if it is an actual fix but I'm going to test and see if it solves the problem
        sleep(240)

        with requests.Session() as session:
            session.auth = (username, password)
            r1 = session.request('get', URL)

            r = session.get(r1.url, auth=(username, password))
            r.raise_for_status()
            f = open(FILENAME,'wb')
            f.write(r.content)
            f.close()
            print('contents of URL written to '+FILENAME)

def import_data():
    os.chdir(merra_data_path)
    for filename in os.listdir(merra_data_path):
        f = cdf.Dataset(filename)

    lat_get = f.variables['lat']
    lon_get = f.variables['lon']
    lev_get = f.variables['lev']
    time_get = f.variables['time']

    tot_u_get = f.variables['U']
    tot_v_get = f.variables['V']
    t_get = f.variables['T']

    lat = lat_get[:]
    lon = lon_get[:]
    lev = lev_get[:]
    time = time_get[:]
    tot_u = tot_u_get[:]
    tot_v = tot_v_get[:]
    t = t_get[:]

    return lat, lon, lev, time, tot_u, tot_v, t, filename

def import_mean():
    os.chdir('/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bam')
    f = cdf.Dataset('climatological_mean_1980_2000.nc4')

    u_get = f.variables['u_tot']

    u = u_get[:]
    print(u.shape)

def call_main(month, year, day):
    test = download_file(month, year, day)
    if test[0] == 'n':
        #* Have to return 5 no datas because there is 5 values passed out in a normal call.
        return 'no data', 'no data', 'no data', 'no data', 'no data', 'no data', 'no data'
    lat, lon, lev, time, tot_u, tot_v, t, filename = import_data()
    os.chdir(merra_data_path)
    os.remove(filename)

    return lat, lon, lev, time, tot_u, tot_v, t

if __name__ == '__main__':
    #download_file('08','1980')
    #lat, lon, lev, tot_u, tot_v, filename = import_data()
    import_mean()



    
    


