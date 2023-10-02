"""
*This code is the first code to look into making eofs for SH variability in the values
we're looking at.

Leave notes:
    Not sure why the cloudsat_occurrence EOF analysis isn't working. First an issue with the 
    number of samples, second an issue with infinite values in x? not sure what that means...

To do:
    Probably should recreate the plots in the Thompson paper to test my eof analysis.

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
#from concurrent.futures.process import _MAX_WINDOWS_WORKERS
from sklearn.decomposition import PCA
from eofs.standard import Eof
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
import bam_analysis
import make_bam_index_smal
import import_eke_data
import make_raw_bam_anim

eke_data_path = '/uufs/chpc.utah.edu/common/home/mace-group6/merra_monthly_data/BAM/data4bamtest'

years_arr = np.arange(2006, 2017)
bad_year = np.argwhere(years_arr==2011)
years_arr = np.delete(years_arr, bad_year, 0)
years = go_stats.make_np_list(years_arr)
months = ['01', '02', '03']

def make_eof(data):
    pca_model = PCA().fit(data)
    eofs = pca_model.components_
    test = pca_model.transform(data)
    score = pca_model.singular_values_
    
    print(data.shape)
    print(eofs.shape)
    print(test.shape)
    exit()

    return eofs

def build_cov_arr_1st_dem(data1, data2, corsel, index):
#* This is the same code as oned_corr.build_cov_arr(). I just needed to rework it because for the 
#cloudsat data, the time dimension is the last dimension, but for this data it is the first dimension.
    xstep = np.arange(0,data1.shape[1])
    ystep = np.arange(0,data1.shape[0])
    print(index)
    for y in ystep:
        for x in xstep:
            if corsel == 1:
                cov_m_hold = np.correlate(data1[:,y,x], index)
                cov_hold = cov_m_hold
            elif corsel == 2:
                cov_m_hold = np.corrcoef(data1[y,x,:], data2[y,x,:])
                cov_hold = cov_m_hold[0,1]
            elif corsel == 'reg':
                model, intercept, rval, pval, stderr = sp.linregress(index, data1[y,x,:])
                cov_hold = model
            else:
                cov_m_hold = np.cov(data1[:,y,x], index)
                cov_hold = cov_m_hold[0,1]
            if x == 0:
                cov_x_arr = np.array([cov_hold])
                continue
            cov_x_arr = np.append(cov_x_arr, np.array([cov_hold]), axis=0)
        if y == 0:
            cov_arr = np.array([cov_x_arr])
            continue
        cov_arr = np.append(cov_arr, np.array([cov_x_arr]), axis =0)
    
    return cov_arr

def import_eke_data():
    os.chdir(eke_data_path)
    f = cdf.Dataset('u_anom_test3d_1986.nc4')
    level_get = f.variables['level']
    lat_get = f.variables['lat']
    data_get = f.variables['eke_data']
    vt_data_get = f.variables['vt_data']
    u_data_get = f.variables['u_data']
    uv_data_get = f.variables['uv_data']


    level = level_get[:]
    lat = lat_get[:]
    data = data_get[:]
    vt_data = vt_data_get[:]
    u_data = u_data_get[:]
    uv_data = uv_data_get[:]

    return data, level, lat, vt_data, u_data, uv_data

def make_df(data, col1, col2, col3):
    c1step = np.arange(0,col1.shape[0])
    c2step = np.arange(0,col2.shape[0])
    c3step = np.arange(0,col3.shape[0])
    for c1 in c1step:
        for c2 in c2step:
            for c3 in c3step:
                if c1 == 0 and c2 == 0 and c3 == 0:
                    mast_arr = np.array([[col1[c1], col2[c2], col3[c3], data[c1,c2,c3]]], dtype=object)
                    continue
                mast_arr = np.append(mast_arr, np.array([[col1[c1], col2[c2], col3[c3], data[c1,c2,c3]]], dtype=object), axis=0)
                print(mast_arr.shape)

    return mast_arr

def make_cols(df_len, col, col_leng):
    arr_hold = np.ones(df_len)
    col_step = 0
    tstep = np.arange(0, df_len, col_leng)
    for t in tstep:
        if t == 0:
            t_hold = t
            continue
        input = col[col_step]
        arr_hold[t_hold:t] = arr_hold[t_hold:t]*input
        t_hold = t
        col_step=col_step+1
    
    arr_hold[t:] = arr_hold[t:]*col[-1]
    return arr_hold 

def repeater(big_col, col_send):
    tstep = np.arange(0, big_col.shape[0], col_send.shape[0])
    for t in tstep:
        if t == 0:
            col_arr = col_send
            continue
        col_arr = np.append(col_arr, col_send, axis=0)

    return col_arr

def add_df_vals(df, col1, col2, col3, label1, label2, label3):
    #*One of the most confusing pieces of code I've written. 
    #   See the readme file for a more in depth overview.
    c1step = np.arange(0,df.shape[0], col1.shape[0])
    #First dimension is the biggest, doesn't repeat
    col1_send = make_cols(df.shape[0], col1, col2.shape[0]*col3.shape[0])
    #Second dimension will repeat every new 1st dimension
    col2_hold = make_cols(int(df.shape[0]/col1.shape[0]), col2, col3.shape[0])
    #Third dimension will also repeat every new 1st dimension and every new 2nd dimension
    #Think stacking rows vertically on one another to build a picture.
    col3_hold = make_cols(int(df.shape[0]/(col1.shape[0]*col2.shape[0])), col3, 1)

    col2_send = repeater(col1_send, col2_hold)
    col3_send = repeater(col1_send, col3_hold)

    df[str(label1)] = col1_send
    df[str(label2)] = col2_send
    df[str(label3)] = col3_send

    return df

def bam_eofs():
    data, level, lat, vt_data = import_eke_data()
    time = np.arange(0,180)
    df = pd.DataFrame(data.flatten(), columns=['Data'])
    df_send = add_df_vals(df, time, level, lat, 'Time', 'Level', 'Lat')

    eof_send = get_time_space(df_send, time_dim='Time', lumped_space_dims=['Level','Lat'])
    n = 4
    pca = df_eof(eof_send, pca_type='varimax', n_components=n)
    eofs = pca.eofs(2, n)
    eofs_da = eofs.stack(['Level','Lat']).to_numpy()
    pcs = pca.pcs(s=2, n=n)
    evfs = pca.evf(n=n)

    eof1 = eofs.iloc[0]
    print(eof1.head)
    eof1_test = eof1.to_frame()
    eof1_arr = eof1_test.unstack().to_numpy()

    #print(eofs.iat[0, 10])

    print(eof1_arr[0,:])


    pcs = pcs.to_numpy()

    return eof1_arr, pcs, vt_data, data

def bam_eof2():
    time = np.arange(0,180)
    lat_height_arr = [-70,-20.5,1000,200]
    eke_data, level, lat, vt_data, u_data, uv_data = import_eke_data()
    bounds_arr = jja_global_avg.restrict_domain(level,lat,lat_height_arr[-1],lat_height_arr[2],lat_height_arr[1],
                                    lat_height_arr[0])
    eke_data, level_rest, lat_rest = jja_global_avg.region_arr3d(eke_data, level, lat, bounds_arr)
    u_data, level_rest, lat_rest = jja_global_avg.region_arr3d(u_data, level, lat, bounds_arr)
    vt_data, level_rest, lat_rest = jja_global_avg.region_arr3d(vt_data, level, lat, bounds_arr)
    uv_data, level_rest, lat_rest = jja_global_avg.region_arr3d(uv_data, level, lat, bounds_arr)

    time_steps, spatial_points = eke_data.shape[0], np.prod(eke_data.shape[1:])
    reshaped_data = eke_data.reshape(time_steps, spatial_points)

    solver = Eof(reshaped_data)
    eofs = solver.eofs(neofs=4)
    eigenvalues = solver.varianceFraction(neigs=4)
    pcs = solver.pcs(npcs=4)
    print('eigenvalues:')
    print(eigenvalues)

    lev_points, lat_points = eke_data.shape[1], eke_data.shape[-1]

    eof1 = eofs[0]
    eof1_reshaped = eof1.reshape(lev_points,lat_points)
    
    return eof1_reshaped, pcs, vt_data, eke_data, uv_data, level_rest, lat_rest

def myeof(x):
    print(x.shape)
    n = x.shape[0]
    x_t = np.transpose(x)
    print(x_t.shape)
    cov_mtx = (x_t@x)/n-1
    cov_mtx_test = np.cov(x, rowvar=False, ddof=1)
    print(cov_mtx_test.shape)
    print(cov_mtx.shape)
    print(np.mean(cov_mtx_test))
    print(np.mean(cov_mtx))
    
    eig_val, eig_vec = np.linalg.eigh(cov_mtx_test)

    sorted_index = np.argsort(eig_val)[::-1]

    sorted_eig_val = eig_val[sorted_index]
    sorted_eig_vec = eig_vec[sorted_index]

    eof1 = sorted_eig_vec[0,:]
    eof1_r = np.reshape(eof1, (19,100))

    eofs = sorted_eig_vec[:,0:6]
    test = np.dot(eofs.transpose(), x.transpose()).transpose()
    print(test.shape)
    print(sorted_eig_val[0]/np.sum(sorted_eig_val))

    return eof1_r, test[:,0]

def remove_mean(eke_data_s, eke_mean_s):
    tstep = np.arange(0,eke_data_s.shape[0])

    for t in tstep:
        eke_data_anom = eke_data_s[t,:,:]-eke_mean_s
        if t == 0:
            eke_anom_arr = np.array([eke_data_anom])
            continue
        eke_anom_arr = np.append(eke_anom_arr, np.array([eke_data_anom]), axis=0)

    return eke_anom_arr
    
def bam_eof_3():
    #Need to restrict the levels to below 200 hPa
    lat_height_arr = [-70,-20.5,900,400]
    #eke_data, level, lat, vt_data, u_data, uv_data = import_eke_data()
    level, lat, lon, u_anom, u3d, v_anom, v3d, eke_data, eke3d = import_eke_data.import_data3d()
    lat_m, lev_m, u_m, v_m, t_m = make_bam_index_smal.get_seasonal_mean('jja')
    # bounds_arr = jja_global_avg.restrict_domain(level,lat,lat_height_arr[-1],lat_height_arr[2],lat_height_arr[1],
    #                                 lat_height_arr[0])
    # eke_data, level_rest, lat_rest = jja_global_avg.region_arr3d(eke_data, level, lat, bounds_arr)
    # u_data, level_rest, lat_rest = jja_global_avg.region_arr3d(u_data, level, lat, bounds_arr)
    # vt_data, level_rest, lat_rest = jja_global_avg.region_arr3d(vt_data, level, lat, bounds_arr)
    # uv_data, level_rest, lat_rest = jja_global_avg.region_arr3d(uv_data, level, lat, bounds_arr)
    # t_data_r, level_rest, lat_rest= jja_global_avg.region_arr(t_m, level, lat, bounds_arr)

    # print(eke_data.shape)
    # eke_data_f = bam_analysis.make_eke_filter_test(eke_data, lat_rest, level_rest, t_data_r)
    # u_data_f = bam_analysis.make_eke_filter_test(u_data, lat_rest, level_rest, t_data_r)
    # print(eke_data_f.shape)
    eke3d = eke3d*2
    eke3d_m = np.mean(eke3d, axis=1)
    eke3d_m = eke3d_m/2
    u3d_m = np.mean(u3d, axis=1)
    mass_weights = go_stats.mass_weight_3(level, eke3d[10,5,:,:])
    eke3d_w, uanom_w = make_raw_bam_anim.apply_weights(eke3d_m, u3d_m, mass_weights)
    #eke_m = np.mean(eke3d_w, axis=0)

    #eke_m_r = remove_mean(eke3d_w, eke_m)

    time_steps, level_steps, lat_steps = eke3d_w.shape
    eke_data_send = np.reshape(eke3d_w, (time_steps, (level_steps*lat_steps)))
    #u_data_send = u_data.reshape(time_steps, level_steps*lat_steps)
    print(eke_data_send.shape)

    time = np.arange(0,180)

    pca = PCA(n_components=6)
    pca.fit(eke_data_send)
    eofs = pca.components_
    pcs = pca.transform(eke_data_send)

    eofs_reshaped = eofs[0].reshape(level_steps, lat_steps)
    #print(level_rest.shape)
    print(pca.explained_variance_ratio_)

    return eofs_reshaped, pcs, eke3d_w, uanom_w, v_anom, level, lat 
    #return eofs_reshaped, pcs, vt_data, u_data, eke_data, uv_data, level_rest, lat_rest

def go_eof_cld(cld_data):
    time_steps, lon_steps, lat_steps = cld_data.shape
    cld_data_send = cld_data.reshape(time_steps, lat_steps*lon_steps)

    pca = PCA(n_components=6)
    pca.fit(cld_data_send)
    eofs = pca.components_
    pcs = pca.transform(cld_data_send)

    print(pca.explained_variance_ratio_)

    eofs_reshaped = eofs[0].reshape(lon_steps, lat_steps)

    return eofs_reshaped, pcs

def eof_main(data_s, n_components_s):
    time_shape, y_shape, x_shape = data_s.shape
    data_send = np.reshape(data_s, (time_shape,(y_shape*x_shape)))
    
    pca = PCA(n_components = n_components_s)
    pca.fit(data_send)
    eofs = pca.components_
    pcs = pca.transform(data_send)

    return eofs, pcs, pca.explained_variance_, y_shape, x_shape

#*This is the eof analysis for the BAM using a different EOF calculator:
if __name__ == '__main__':
    bam_eof_3()
    exit()

if __name__ == '__main__':
    eke_data, level, lat, vt_data, u_data = import_eke_data()
    time = np.arange(0,180)
    print(eke_data.shape)
    time_steps, spatial_points = eke_data.shape[0], np.prod(eke_data.shape[1:])
    reshaped_data = eke_data.reshape(time_steps, spatial_points)

    solver = Eof(reshaped_data)
    eofs = solver.eofs(neofs=4)
    eigenvalues = solver.varianceFraction(neigs=4)
    pcs = solver.pcs(npcs=4)

    lev_points = 42
    lat_points = 100

    eof1 = eofs[0]
    eof1_reshaped = eof1.reshape(lev_points,lat_points)

    p = plt.imshow(eof1_reshaped, cmap='RdBu', origin='upper')
    plt.colorbar(p)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/plots/bam_reg/tester.png')
    exit()


#*This is the eof analysis for the BAM
# if __name__ == '__main__':
#     eke_data, level, lat = import_eke_data()
#     #1st, convert to dataframe
#     #2nd, use get_time_space function to convert to the eof
#     #3rd, 

#* This is the eof analysis on the bam occurrence dataset.
# if __name__ == '__main__':
#     data, level, lat, vt_data = import_eke_data()
#     time = np.arange(0,180)
#     df = pd.DataFrame(data.flatten(), columns=['Data'])
#     df_send = add_df_vals(df, time, level, lat, 'Time', 'Level', 'Lat')

#     eof_send = get_time_space(df_send, time_dim='Time', lumped_space_dims=['Level','Lat'])
#     n = 4
#     pca = df_eof(eof_send, pca_type='varimax', n_components=n)
#     eofs = pca.eofs(2, n)
#     eofs_da = eofs.stack(['Level','Lat']).to_numpy()
#     pcs = pca.pcs(s=2, n=n)
#     evfs = pca.evf(n=n)

#     print(evfs)
#     print(pcs.shape)

#     eof1 = eofs.iloc[0]

#     eof1_test = eof1.to_frame()

#     eof1_arr = eof1_test.unstack().to_numpy()

#     pcs = pcs.to_numpy()
#     #eofs = eofs.to_numpy()

#     mast_plot.bam_eof(pcs[:,0], eof1_arr)

#*This was used to make the EKE animation plots. I didn't know where to put it.
# if __name__ == '__main__':
#     data, level, lat, vt_data = import_eke_data()
#     tstep = np.arange(0,data.shape[0])
#     for t in tstep:
#         data_plt = data[t,:,:]

#         mast_plot.eke_vis(data_plt, t)

#*This is the code for looking at the EOF of cloud occurrence data.
if __name__ == '__main__':
    cld_occ_data, lat, lon = oned_corr.eof_call_main(years, months, 0)

    lat = lat[:-1]
    lon = lon[:-1]

    df = pd.DataFrame(cld_occ_data.flatten(), columns=['Data'])
    df_send = add_df_vals(df, lat, lon, years_arr, 'Lat', 'Lon', 'Years')

    eof_send = get_time_space(df_send, time_dim='Years', lumped_space_dims=['Lat','Lon'])
    n = 4
    pca = df_eof(eof_send, n_components=n)
    eofs = pca.eofs(2, n)
    eofs_da = eofs.stack(['Lat','Lon']).to_numpy()
    pcs = pca.pcs(s=2, n=n)
    evfs = pca.evf(n=n)

    pcs = pcs.to_numpy()
    eof1 = eofs.iloc[0]
    eof1_test = eof1.to_frame()
    eof1_arr = eof1_test.unstack().to_numpy()

    print(evfs)
    #plt.plot(pcs[:,0])
    #p = plt.contourf(lon, lat, eof1_arr, cmap='RdBu')
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1, projection=proj)
    ax.coastlines()
    p = ax.imshow(eof1_arr, cmap='RdBu', origin='lower', aspect='equal', extent=[-180,180,-90,90])
    plt.colorbar(p, ax=ax)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1113223/grad_research/my_code/climate_proj/tester.png')


#* This was looking at taking the EOF of the merra southern ocean data. Not sure if I ever really
# got anywhere with it.

# if __name__ == '__main__':
#     m_mast_arr = merra_analize.build_ds('jja_mean')
#     h_mast_arr = merra_analize.build_ds('high')
#     l_mast_arr = merra_analize.build_ds('low')
#     m_lat = m_mast_arr[0,0]
#     m_lon = m_mast_arr[0,1]
#     m_mb500, m_slp, m_t2m, m_pbltop, m_q500 = merra_analize.get_mean_vals(m_mast_arr)
#     h_mb500, h_slp, h_t2m, h_pbltop, h_q500 = merra_analize.get_mean_vals(h_mast_arr)
#     l_mb500, l_slp, l_t2m, l_pbltop, l_q500 = merra_analize.get_mean_vals(l_mast_arr)
#     years = np.arange(2006,2017)
#     ts_data = merra_analize.build_ts(years, 'jja_mean', m_lat, m_lon)
#     eofs = make_eof(ts_data)

#     mast_plot.eof_inital(m_lat, m_lon, eofs)
