import datetime
import itertools
import obspy
import funciones_linda as fl
import os
from datetime import datetime
import multiprocessing
from scipy.special import jv,jn_zeros
from scipy.optimize import curve_fit
from scipy import interpolate
import matplotlib.pyplot as plt
from pandas import DataFrame,read_csv
import numpy as np
from scipy.fft import rfft, fft, ifft, fftshift, fftfreq,rfftfreq
from scipy.signal import detrend, correlate, correlation_lags, tukey, butter, sosfilt, \
    welch,csd,resample,coherence
import datetime
from functools import partial
from itertools import product
import time
from datetime import datetime
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,LinearOperator,aslinearoperator
from collections import Counter
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
import cc_processing as cc
step0=False ## read each mseed and obtain coherence function
step1=False ## create .txt with coherences
step2=False ## test each coherence to find proper freq. limits
step2d5=False
step3=True ## run all coherences with correct freq. limits


Dic21=True
Hs21=False
Mar22=False
Sept22=False
Dic22=False
from scipy import signal
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt



folders=['Dic21','HS21','Mar22','Sept22','Dic22']
#folders=['Dic21']
#folders=['HS21']
#folders=['Mar22']
#folders=['Sept22']
folders_low=[x.lower() for x in folders]
nc=0
path_ori = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'
path_out= '/home/doctor/Doctor/Magister/Tesis/databases/process_data/AJUSTES_NUEVOS/'
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots(2)
coords_x=[]
coords_y=[]
estado=[]
for num,folder in enumerate(folders):

    nf=500 ## number of files to read per folder
    fold_low=folders_low[num]

    path=path_ori+folder+'/STACKS/'
    files = sorted(os.listdir(path))
    coordenadas = np.loadtxt(path[:-7] + 'coords_all.txt')
    mediciones=[]
    with open(path[:-7] + 'values_nodes.txt','r') as fileobject:
        lines=fileobject.readlines()
        for line in lines:
            split=line.strip().split(' ')
            mediciones.append('_'.join(split))

    coordenadas_todas=np.loadtxt(path_ori + 'coords_all_sets.txt')
    coord_split = np.split(coordenadas_todas, 2, axis=1)
    coord_uni = np.unique(np.vstack((coord_split[0], coord_split[1])), axis=0)
    lims = np.loadtxt(path[:-7] + 'lims_cc_new_jan.txt')

    cps = []
    fxs = []
    wnews = []

    idx=np.where(lims[:,3]> 1e-8) ## search for resolved dispersion curves

    cohs_new=['par'+str(x).zfill(2)+'.txt' for x in idx[0] ]
    lims=lims[:nf]
    lims=np.array([x for x in lims if int(x[2]) != 0])
    nfold=0
    #coordenadas=coordenadas[idx]

    for lim, file in zip(lims, cohs_new):
        data=np.loadtxt(path+file)
        f=data[:,0]
        coh=data[:,1]
        Cx=coh[idx]
        nfold=int(lim[0])
        f1=lim[1]
        f2=lim[2]
        f3=lim[3]
        dist=lim[4]
        smooth_factor=lim[5]
        med=mediciones[nfold]
        kind=lim[6]
        sigmad=1 ## varianza de los datos
        sigma=0.75 ## varianza de la suavizacion
        alpha=1 ## varianza de la solucion inicial
        p1=coordenadas[nfold][0:2]
        p2=coordenadas[nfold][2:4]

        cms, zcs, wzero, r, jn, c_est, A_est, rho_pre, fx, CCX, CCX_smooth,fxp,CCXp,CCXp_smooth = cc.get_dispersion_curve2(f, coh, dist,
                                                                                                      int(smooth_factor), f1,
                                                                                                      f2, nc, coord_uni,
                                                                                                      np.vstack(
                                                                                                          (p1, p2)),
                                                                                                      create_figures=False)

        idx2 = np.where(fx <= f3)
        fx2 = fx[idx2]
        CCX2_smooth = CCX_smooth[idx2]
        rho_pre2 = rho_pre[idx2]
        c_est2 = c_est[idx2]
        zcs2 = [x for x in zcs if x <= idx2[0][-1]]
        wzero2 = wzero[:len(zcs2)]

        std = np.abs(CCX2_smooth - rho_pre2)
        std = std * np.eye(len(std))
        std_inv = np.linalg.inv(std)

        rhop, cp, Ap, H2, CM, err, G, G1, G2 = cc.get_dispersion_curve_newton22(fx2 * 2 * np.pi, dist, rho_pre2, CCX2_smooth, fxp,
                                                                               CCXp,
                                                                               CCXp_smooth, c_est2, A_est, sigmad, sigma,
                                                                               alpha, nc,nfold, coord_uni, np.vstack((p1, p2)), wzero2,
                                                                               zcs2,std_inv,fold_low,med)

        lembdas=cp/fx2
        print('dist',dist)
        print('lambdamax',lembdas[0])
        print('lambdamin',lembdas[-1])
        print('lambda/r',lembdas[0]/dist)
        print('lambda/r',lembdas[-1]/dist)
        cps.append(cp)
        wnews.append(fx2 * 2 * np.pi)
        coords_x.append(p1)
        coords_y.append(p2)
        estado.append(int(kind))

        nfold+=1
        nc+=1
        print(nfold,nc)
        print(f1,f2,dist)
        if int(kind) == 1:
            ax1.plot(fx2, cp, 'g-', linewidth=0.5)
            ax2[0].plot((f1,f2),(dist,dist),'g-')
            ax2[1].plot((f1,f2),(cp[0]-cp[-1],cp[0]-cp[-1]),'g-')
        if int(kind) == 2:
            ax1.plot(fx2, cp, 'y-', linewidth=0.5)
            ax2[0].plot((f1,f2),(dist,dist),'y-')
            ax2[1].plot((f1,f2),(cp[0]-cp[-1],cp[0]-cp[-1]),'y-')

        if int(kind) == 3:
            ax1.plot(fx2, cp, 'r-', linewidth=0.5)
            ax2[0].plot((f1,f2),(dist,dist),'r-')
            ax2[1].plot((f1,f2),(cp[0]-cp[-1],cp[0]-cp[-1]),'r-')
np.savetxt(path_out+'estado.txt',estado, fmt='%i')

