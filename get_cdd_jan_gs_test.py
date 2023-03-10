import datetime
import itertools
import obspy
import funciones_linda as fl
import os
from datetime import datetime
import multiprocessing
from scipy.special import jv,jn_zeros,j0
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
import matplotlib
matplotlib.use('Agg')


folders=['Dic21','HS21','Mar22','Sept22','Dic22']
#folders=['Dic21','HS21','Mar22','Dic22']
#folders=['Dic21']
#folders=['HS21']
#folders=['Mar22']
#folders=['Sept22']
folders_low=[x.lower() for x in folders]
nc=0
path_ori = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'
path_out= '/home/doctor/Doctor/Magister/Tesis/databases/process_data/AJUSTES_NUEVOS_GS_FINAL/'
#fig1, ax1 = plt.subplots()
#fig2, ax2 = plt.subplots(2)
coords_x=[]
coords_y=[]
estado=[]
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots(2)
for num,folder in enumerate(folders):

    nf=200 ## number of files to read per folder
    fold_low=folders_low[num]

    path=path_ori+folder+'/STACKS_2/'
    files = sorted(os.listdir(path))
    coordenadas = np.loadtxt(path[:-9] + 'coords_all2.txt')
    mediciones = np.genfromtxt(path[:-9] + 'values_separado.txt', dtype=str)

    coordenadas_todas=np.loadtxt(path_ori + 'coords_all_sets2.txt')
    coord_split = np.split(coordenadas_todas, 2, axis=1)
    coord_uni = np.unique(np.vstack((coord_split[0], coord_split[1])), axis=0)
    lims = np.loadtxt(path[:-9] + 'lims_cc_new_feb.txt')

    cps = []
    fxs = []
    wnews = []

    idx=np.where(lims[:,3]> 1e-8) ## search for resolved dispersion curves

    cohs_new=['par'+str(x).zfill(3)+'.txt' for x in idx[0] ]
    lims = lims[:nf]
    lims = np.array([x for x in lims if int(x[2]) != 0])
    nfold=0
    #sel=0
    #nc=0
    #nfold=0
    #file=cohs_new[sel]
    #lim=lims[sel]
    for lim,file in zip(lims,cohs_new):
        data = np.loadtxt(path + file)
        f = data[:, 0]
        coh = data[:, 1]
        Cx = coh[idx]
        nfold = int(lim[0])
        f1 = lim[1]
        f2 = lim[2]

        f3 = lim[3]
        #f2=20
        #f3=20
        dist = lim[4]
        smooth_factor = lim[5]
        med = mediciones[nfold]
        kind = lim[6]
        sigmad = 1  ## varianza de los datos
        sigma = 0.75  ## varianza de la suavizacion
        alpha = 1  ## varianza de la solucion inicial
        p1 = coordenadas[nfold][0:2]
        p2 = coordenadas[nfold][2:4]


        cmin=90
        cmax=320
        fx,CCX,CCX_smooth,c_combs,ftrial,wtrial,Ab,rho_pre,errb,ccb,cb_int,fxp,CCXp,CCXp_smooth = \
                            cc.get_dispersion_curve2_gs(f, coh, dist,int(smooth_factor),f1,f2,
                            cmax,cmin, nc, coord_uni,np.vstack((p1, p2)),create_figures=False)


        print('err2',np.linalg.norm(CCX_smooth-rho_pre))
        std = np.abs(CCX_smooth - rho_pre)
        std = std * np.eye(len(std))
        std_inv = np.linalg.inv(std)

        rhop, cp, Ap, H2, CM, err, G, G1, G2 = cc.get_dispersion_curve_newton22(fx * 2 * np.pi, dist, rho_pre, CCX_smooth, fxp,
                                                                                    CCXp,
                                                                                    CCXp_smooth, cb_int, Ab, sigmad, sigma,
                                                                                    alpha, nc,nfold, coord_uni, np.vstack((p1, p2)), wtrial,
                                                                                    ccb,std_inv,fold_low,med)

        lembdas = cp / fx
        print('dist', dist)
        print('lambdamax', lembdas[0])
        print('lambdamin', lembdas[-1])
        print('lambda/r', lembdas[0] / dist)
        print('lambda/r', lembdas[-1] / dist)
        cps.append(cp)
        wnews.append(fx * 2 * np.pi)
        coords_x.append(p1)
        coords_y.append(p2)
        estado.append(int(kind))

        nfold += 1
        nc += 1
        print(nfold, nc)
        print(f1, f2, dist)
        if int(kind) == 1:
            ax1.plot(fx, cp, 'g-', linewidth=0.5)
            ax2[0].plot((f1, f2), (dist, dist), 'g-')
            ax2[1].plot((f1, f2), (cp[0] - cp[-1], cp[0] - cp[-1]), 'g-')
        if int(kind) == 2:
            ax1.plot(fx, cp, 'y-', linewidth=0.5)
            ax2[0].plot((f1, f2), (dist, dist), 'y-')
            ax2[1].plot((f1, f2), (cp[0] - cp[-1], cp[0] - cp[-1]), 'y-')

        if int(kind) == 3:
            ax1.plot(fx, cp, 'r-', linewidth=0.5)
            ax2[0].plot((f1, f2), (dist, dist), 'r-')
            ax2[1].plot((f1, f2), (cp[0] - cp[-1], cp[0] - cp[-1]), 'r-')

plt.savefig(path_out+'Todas.png',format='png', dpi=300, bbox_inches='tight',pad_inches=0.1)
    # idx2 = np.where(fx <= f3)
    # fx2 = fx[idx2]
    # CCX2_smooth = CCX_smooth[idx2]
    # rho_pre2 = rho_pre[idx2]
    # c_est2 = c_est[idx2]
    # zcs2 = [x for x in zcs if x <= idx2[0][-1]]
    # wzero2 = wzero[:len(zcs2)]
    #
    # std = np.abs(CCX2_smooth - rho_pre2)
    # std = std * np.eye(len(std))
    # std_inv = np.linalg.inv(std)
    #
    # rhop, cp, Ap, H2, CM, err, G, G1, G2 = cc.get_dispersion_curve_newton22(fx2 * 2 * np.pi, dist, rhob, CCX2_smooth, fxp,
    #                                                                             CCXp,
    #                                                                             CCXp_smooth, c_est2, A_est, sigmad, sigma,
    #                                                                             alpha, nc,nfold, coord_uni, np.vstack((p1, p2)), wzero2,
    #                                                                             zcs2,std_inv,fold_low,med)

    # cms, zcs, wzero, r, jn, c_est, A_est, rho_pre, fx, CCX, CCX_smooth,fxp,CCXp,CCXp_smooth = cc.get_dispersion_curve2(f, coh, dist,
    #                                                                                                   int(smooth_factor), f1,
    #                                                                                                   f2, nc, coord_uni,
    #                                                                                                   np.vstack(
    #                                                                                                       (p1, p2)),
    #                                                                                                   create_figures=False)

    # ip=interpolate.interp1d(wtrial,c_combs,'linear',fill_value='extrapolate')
    # wx=fx*2*np.pi
    # ipw=ip(wx).astype(np.float32)
    # rho_pre2 = j0(0, wx * r / ipw)
    # # rho_pre  = J0[np.floor(DJ * ((wx * r / ipw) - jargmin)).astype(int)]
    # A_est = (rho_pre2 @ CCX_smooth2) / (np.linalg.norm(rho_pre2, axis=1)) ** 2
    # A_est_t = A_est.reshape(len(A_est), 1)
    # rho_pre2 = np.multiply(A_est_t, rho_pre2)
    # Drho = CCX_smooth - rho_pre2
    # err = np.linalg.norm(Drho, axis=1) ** 2
    # agmin = np.argmin(err)
    # Asol = A_est[agmin]
    # rhosol = rho_pre2[agmin]
    # errsol = err[agmin]
    # csol = ipw[agmin]
    # ccombsol = c_combs[agmin]
    # fig,ax=plt.subplots(2,1)
    # ax[0].plot(fx,CCX_smooth,'k')
    # ax[0].plot(fx,rho_pre,'r')
    # ax[0].plot(fx,rhosol,'g')
    #
    # ax[1].plot(fx, c_est, 'k')
    # #ax[1].plot(fx, rho_pre, 'r')
    # ax[1].plot(ftrial, ccombsol, 'g')

    #print('best', best_sol)
    # idx2 = np.where(fx <= f3)
    # fx2 = fx[idx2]
    # CCX2_smooth = CCX_smooth[idx2]
    # rho_pre2 = rho_pre[idx2]
    # c_est2 = c_est[idx2]
    # zcs2 = [x for x in zcs if x <= idx2[0][-1]]
    # wzero2 = wzero[:len(zcs2)]
    #
    # std = np.abs(CCX2_smooth - rho_pre2)
    # std = std * np.eye(len(std))
    # std_inv = np.linalg.inv(std)
    #
    # rhop, cp, Ap, H2, CM, err, G, G1, G2 = cc.get_dispersion_curve_newton22(fx2 * 2 * np.pi, dist, rho_pre2, CCX2_smooth, fxp,
    #                                                                            CCXp,
    #                                                                            CCXp_smooth, c_est2, A_est, sigmad, sigma,
    #                                                                            alpha, nc,nfold, coord_uni, np.vstack((p1, p2)), wzero2,
    #                                                                            zcs2,std_inv,fold_low,med)
    #
    # lembdas=cp/fx2
    # print('dist',dist)
    # print('lambdamax',lembdas[0])
    # print('lambdamin',lembdas[-1])
    # print('lambda/r',lembdas[0]/dist)
    # print('lambda/r',lembdas[-1]/dist)
    # cps.append(cp)
    # wnews.append(fx2 * 2 * np.pi)
    # coords_x.append(p1)
    # coords_y.append(p2)
    # estado.append(int(kind))




