"""
*This code will clean the merra code to only include the specified region and specified variables.

To do:
    Need to change the data call to get the MERRA data that I was using for the ENSO analysis.
    I pulled the current data call from the bam_merra_calc.py script because I lost the old save of 
    merra_first.py in the catastrophy.
Leave notes:
    Was focused on looking at the meterology of the area correlated with the BAM. I've restricted the data but need to 
    re-work how it integrates into all the other code that made it. 

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
import requests
from time import sleep

merra_data_path = '/uufs/chpc.utah.edu/common/home/u1113223/merra_hold_enso'

month = '11'
year = '1999'
day = '12'

latlon_bounds = [75,35,180,-180]

def download_file(month, year, day):
    #print('starting '+month+'/'+year)
    os.chdir('/uufs/chpc.utah.edu/common/home/u1113223/merra_hold_enso')
    if int(year) < 1992:
        URL = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXSLV.5.12.4/'+ str(year) + '/MERRA2_100.tavgM_2d_slv_Nx.'+str(year)+str(month)+'.nc4'
    elif int(year) < 2001:
        URL = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXSLV.5.12.4/'+ str(year) + '/MERRA2_200.tavgM_2d_slv_Nx.'+str(year)+str(month)+'.nc4'
    elif int(year) < 2011:
        URL = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXSLV.5.12.4/'+ str(year) + '/MERRA2_300.tavgM_2d_slv_Nx.'+str(year)+str(month)+'.nc4'
    else:
        URL = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXSLV.5.12.4/'+ str(year) + '/MERRA2_400.tavgM_2d_slv_Nx.'+str(year)+str(month)+'.nc4'
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
    # for data in data_sel:
    #     data_get = f.variables[data]
    #     data = np.array([data_get[:]])
    #     print(data.shape[0])

    #     try:
    #         data_arr = np.append(data_arr, np.array([data]), axis=1)
    #     except UnboundLocalError:
    #         data_arr = np.array([data])
    #         print(data_arr.shape[0])

    lat_get = f.variables['lat']
    lon_get = f.variables['lon']
    mb500_get = f.variables['H500']
    slp_get = f.variables['SLP']
    t2m_get = f.variables['T2M']
    pbltop_get = f.variables['PBLTOP']

    lat = lat_get[:]
    lon = lon_get[:]
    mb500 = mb500_get[:]
    slp = slp_get[:]
    t2m = t2m_get[:]
    pbltop = pbltop_get[:]

    return lat, lon, mb500, slp, t2m, pbltop

def test(data, lat, lon):
    #plt.imshow(data, origin='lower', aspect='equal', extent=[-180,180,-90,90])
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1, projection=proj)
    ax.coastlines()
    p = plt.contourf(lon, lat, data, cmap='jet')
    c = plt.colorbar(p, ax=ax)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/tester.png')

if __name__ == '__main__':
    data_get = download_file(month, year, day)
    #if data_get[0] == 'n':
    lat, lon, mb500, slp, t2m, pbltop = import_data()
    print(pbltop.shape)
    bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0],latlon_bounds[1],latlon_bounds[2],
                                                latlon_bounds[3])
    mb500_r, lat_r, lon_r = jja_global_avg.region_arr(mb500, lat, lon, bounds_arr)
    slp_r, lat_r, lon_r = jja_global_avg.region_arr(slp, lat, lon, bounds_arr)
    t2m_r, lat_r, lon_r = jja_global_avg.region_arr(t2m, lat, lon, bounds_arr)
    pbltop_r, lat_r, lon_r = jja_global_avg.region_arr(pbltop[0,:,:], lat, lon, bounds_arr)

# if __name__ == '__main__':
#     download_file('07','1988','24')
#     lat, lon, mb500, slp, t2m = import_data()
#     bound_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0],latlon_bounds[1],
#                                                latlon_bounds[2],latlon_bounds[3])
#     mb500_r, lat_r, lon_r = jja_global_avg.region_arr(mb500, lat, lon, bound_arr)
#     print(mb500_r.shape)


# if __name__ == '__main__':
#     data_sel = ['lat', 'lon']
#     lat, lon, mb500, slp, t2m = import_data(data_sel)
#     #test(t2m[0,:,:], lat, lon)
#     bounds_arr = jja_global_avg.restrict_domain(lat, lon, latlon_bounds[0],latlon_bounds[1],latlon_bounds[2], 
#                                 latlon_bounds[3])
#     test = jja_global_avg.region_arr(t2m[0,:,:], lat, lon, bounds_arr)
#     print(test.shape)




    