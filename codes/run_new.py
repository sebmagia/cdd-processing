import TimeSeries as ts
import os
import multiprocessing
from itertools import product
from functools import partial
from scipy.special import jv
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,LinearOperator,aslinearoperator
import time

import itertools
from scipy.signal import detrend,hann,coherence,correlate,periodogram,correlation_lags,tukey,butter,sosfilt,spectrogram,welch,csd
from scipy import interpolate,integrate
from scipy.fft import fft,ifft,fftshift,ifftshift,fftfreq,rfftfreq,rfft,fft2,irfft

import cc_processing

if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plt

    print('cordero')
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_DENSEARRAY/M4/'


    path_prim='/home/doctor/Doctor/Magister/Tesis/databases/SUR2A/'
    path_prim='/home/doctor/Doctor/Magister/Tesis/databases/Vs_Los_Presidentes/Lin_1_20/'

    #path_prim='/home/doctor/Doctor/Magister/Tesis/databases/process_data/Pent_7_30/'
    path2_prim='/home/doctor/Doctor/Magister/Tesis/databases/process_data/SUR2A/Outputs/'
    xd=filter(os.path.isdir,os.listdir(path_prim))
    carpetas=[path_prim+x for x in os.listdir(path_prim)]
    carpetas=sorted(list(filter(os.path.isdir,carpetas)))
    paths = [x + '/EqualizedFile.dat' for x in carpetas]
    # paths=[x+path_ele[0] for x in paths]
    TS = ts.TimeSeries(paths, mp=False)
    pata=False
    if pata:
        #TS.plot_traces(name2, closefig=False)
        for x in carpetas[0:1]:
            path=path_prim+x
            path2=path2_prim+x
            if 'SUR1A' in path:
                name = 'SUR1A/' + path.split('/')[-1] + '.png'
                name2 = 'SUR1A/' + path.split('/')[-1] + '_TS.png'
                print(name,name2)
            if 'SUR1B' in path:
                name = 'SUR1B/' + path.split('/')[-1] + '.png'
                name2 = 'SUR1B/' + path.split('/')[-1] + '_TS.png'

            if 'SUR2A' in path:
                name = 'SUR2A/' + path.split('/')[-1] + '.png'
                name2 = 'SUR2A/' + path.split('/')[-1] + '_TS.png'

            if 'SUR2B' in path:
                name = 'SUR2B/' + path.split('/')[-1] + '.png'
                name2 = 'SUR2B/' + path.split('/')[-1] + '_TS.png'

            #print(name)
            print(path)
            print(path2)
            path_ele = sorted(os.listdir(path))
            ## check if elements are directories or files
            dire = [os.path.isdir(path+'/'+x) for x in path_ele]
            paths = [path +'/'+x +'/EqualizedFile.dat' for x, y in zip(path_ele, dire) if y]
            #paths=[x+path_ele[0] for x in paths]
            TS=ts.TimeSeries(paths,mp=False)
            TS.plot_traces(name2,closefig=False)

            #del TS
            TS2=ts.TimeSeries(path2,mp=False,create_data=False)
            TS2.get_NCF(TS2.C,TS2.w,name,closefig=False)


            idx = np.where((TS2.w < 30) & (TS2.w > (2)))[0]
            wp = 2 * np.pi * TS.w[idx]
            wmin = wp[0]
            wmax = wp[-1]
            sig=5
            data = TS.C[sig][idx]
            N = len(wp)
            tic = time.time()
            r = TS2.md_ncf.Interstation_Distance[sig]
            cbounds = [80, 300]
            tic = time.time()
            clow = cbounds[0]
            chigh = cbounds[1]
            ## numero de nodos
            NN = 4
            ## numero de frecuencias
            NC = 40
            cgrid = np.linspace(clow, chigh, NC)
            # ns = np.arange(0, 30)
            f = partial(product, repeat=NN)
            g = f(cgrid)
            cnodes = np.split(np.array(list(g)), 256)

            tic = time.time()
            A_best, rho_best, err_best, cnodes_best, c_best, wgrid,w= ts.grid_search_cw2(N, wmin, wmax, data, r, NN, NC, cbounds,
                                                                                 cnodes, True)

        # def tanhfit(omega, e, d, f):
        #     #print('xd',e*omega)
        #     c = d * np.arctanh(e * omega) + f / np.sqrt(omega)
        #     #print(c)
        #
        #     return c
        #
        #
        #
        # from scipy.optimize import curve_fit
        #
        # popt, pcov = curve_fit(tanhfit, w, c_best,bounds=([0,80,0],[wgrid[-1]**-1,400,900]))
        #
        # plt.plot(w,tanhfit(w,popt[0],popt[1],popt[2]))
        # plt.plot(w,c_best)
            toc = time.time()
            print('gridsearch linear Done in {:.4f} seconds'.format(toc - tic))
            rho0 = A_best[0] * jv(0, wf * r / c_best)

            tic = time.time()
            rho, cp, ap = ts.newton_search(wf, data, rho0, c_best, A_best, r)
        # fig,ax=plt.subplots(2,1)
        # ax[0].plot(TS.w[idx],cp)
        # plt.xlabel('Freq [Hz]')
        # plt.ylabel('phase velocity [m/s]')
        # plt.savefig('')
        # #wtod, cc, zc = ts.zero_crossings(wp, data, r)
        # #toc = time.time()
        print('newtonsearch linear Done in {:.4f} seconds'.format(toc - tic))


    #TS = ts.TimeSeries(paths, mp=False)
    #del TS
    # path = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/SUR1A/Outputs/M11'

    #TS = ts.TimeSeries(path2, mp=False, create_data=False)
    #TS.get_NCF(TS.C, TS.w, name)
    #path2 = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/SUR2B/Outputs/M3'

    #path='/home/doctor/Doctor/Magister/Tesis/databases/SUR1B/M4/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/SUR2A/M4/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/SUR2B/M3/'


    #name='M4.png'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_SN_P2/M1/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_SN_S1_PILZ/M5/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_SN_S2_PILZ/M8B/'
    if 'SUR1A' in path:
        name='SUR1A/'+path[-3:-1]+'.png'
    if 'SUR1B' in path:
        name = 'SUR1B/' + path[-3:-1] + '.png'
    if 'SUR2A' in path:
        name = 'SUR2A/' + path[-3:-1] + '.png'
    if 'SUR2B' in path:
        name = 'SUR2B/' + path[-3:-1] + '.png'

