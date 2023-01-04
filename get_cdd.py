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


Dic21=False
Hs21=False
Mar22=False
Sept22=False
Dic22=True
from scipy import signal
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt






path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'

if Dic21:
    path+='Dic21/STACKS/'
    path2='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Dic21/M1B.mseed'
    path3='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Dic21/M2A.mseed'
    nw='SA'
if Hs21:
    path+='HS21/STACKS/'
    nw='SB'
if Mar22:
    path+='Mar22/STACKS/'
    nw='SC'

if Sept22:
    path+='Sept22/STACKS/'
    nw='SD'

if Dic22:
    path += 'Dic22/STACKS/'
    nw = 'SE'

files=sorted(os.listdir(path))
coordenadas=np.loadtxt(path[:-7]+'coords_all.txt')
coord_split=np.split(coordenadas,2,axis=1)
coord_uni=np.unique(np.vstack((coord_split[0],coord_split[1])),axis=0)
cohs=[x for x in files if '.txt' in x]
files=[x for x in files if '.mseed' in x]


if step0:
    i=0
    dist=np.sqrt((coordenadas[i][0]-coordenadas[i][2])**2+(coordenadas[i][1]-coordenadas[i][3])**2)
    fs = 512
    nsegs = 60
    file=files[i]
    stream= obspy.read(path + file, format='MSEED')
    stream.decimate(1)
    ND = [x.stats.station for x in stream]

    #f7, coh7, cohp = fl.cross_coherence(stream, fs, 30, ND, norm='1bit')
    f8, coh8, cohp = fl.cross_coherence(stream, fs, 45, ND, norm='1bit',onesided=True)
    f82, coh82, cohp2 = fl.cross_coherence(stream, fs, 45, ND, norm='1bit',onesided=False)
    #f9, coh9, cohp = fl.cross_coherence(stream, fs, 60, ND, norm='1bit')


    #plt.plot(f9[:5000], coh9[:5000])
    plt.plot(f8[:5000], coh8[:5000],'r')

    plt.plot(f82[:5000], coh82[:5000],'g')
    #plt.plot(f7[:5000], coh7[:5000])

    print('dist',dist)

n=0
if step1:
    nf = 100

    fs = 512

    nsegs = 60
    files = files[:nf]
    n=0
    for file in files:
        dist = np.sqrt((coordenadas[n][0] - coordenadas[n][2]) ** 2 + (coordenadas[n][1] - coordenadas[n][3]) ** 2)

        stream = obspy.read(path + file, format='MSEED')
        stream.decimate(1)
        print(len(stream[0]))
        ND = [x.stats.station for x in stream]
        f,coh,cohp=fl.cross_coherence(stream,fs,nsegs,ND,norm='1bit',onesided=False)
        #xd=fl.cross_coherence(stream,fs,nsegs,ND,norm='1bit',onesided=False)

        idx = np.where((f <= 64))[0]
        plt.plot(f[idx],coh[idx],'k',linewidth=2)
        plt.xlabel('Frequency [Hz] r='+str(dist))
        plt.ylabel('Amplitude')
        plt.savefig(path+file[:-5]+'png',format='png', dpi=300, bbox_inches='tight')
        plt.close()
        np.savetxt(path+file[:-5]+'txt',np.vstack((f,coh)).T)
        n+=1
if step2:
    # read i-th coherence file
    i = 18
    f1 =4
    f2 =26
    f3 =26
    smooth_factor = 100

    cps = []
    fxs = []
    dist=np.sqrt((coordenadas[i][0]-coordenadas[i][2])**2+(coordenadas[i][1]-coordenadas[i][3])**2)
    data=np.loadtxt(path+cohs[i])
    f = data[:, 0]
    coh = data[:, 1]
    fnew=np.linspace(f[0],f[-1],len(f)*5)
    interp=interpolate.interp1d(f,coh,kind='cubic')
    coh2=interp(fnew)
    ## define proper frequency limits and smooth factor

    CC_smooth=cc.smooth(coh,smooth_factor)
    #plt.plot(f,coh)

    #plt.plot(f,CC_smooth)

    p1 = coordenadas[i][0:2]
    p2=coordenadas[i][2:4]
    cms,zcs,wzero,r,jn,c_est,A_est,rho_pre,fx,CCX,CCX_smooth,fxp,CCXp,CCXp_smooth=cc.get_dispersion_curve2(f, coh, dist, smooth_factor, f1, f2,
                                                                                      i, coord_uni, np.vstack((p1,p2)),
                                                                                      create_figures=False)

    idx2=np.where(fx <=f3)
    fx2=fx[idx2]
    CCX2_smooth=CCX_smooth[idx2]
    rho_pre2=rho_pre[idx2]
    c_est2=c_est[idx2]
    zcs2=[x for x in zcs if x<= idx2[0][-1]]
    wzero2=wzero[:len(zcs2)]



    sigmad=1
    sigma=1
    alpha=1
    rhop, cp, Ap, H2, CM, err = cc.get_dispersion_curve_newton(fx2 * 2 * np.pi, dist, rho_pre2, CCX2_smooth, fxp, CCXp,
                                                               CCXp_smooth, c_est2, A_est, sigmad, sigma,
                                                               alpha, i, coord_uni, np.vstack((p1, p2)), wzero2, zcs2)
    #rhop, cp, Ap, H2, CM,err=cc.get_dispersion_curve_newton(fx*2*np.pi, dist,rho_pre, CCX_smooth, fxp,CCXp,CCXp_smooth, c_est, A_est, sigmad, sigma,
    #                                                    alpha,i, coord_uni, np.vstack((p1,p2)),wzero,zcs)

    lembdas=cp/fx2
    print('dist',dist)
    print('lambdamax',lembdas[0])
    print('lambdamin',lembdas[-1])
    print('lambda/r',lembdas[0]/dist)
    print('lambda/r',lembdas[-1]/dist)




if step3:
    cps=[]
    fxs=[]
    wnews=[]
    nf=100
    #lims = np.loadtxt(path[:-7] + 'lims_cc_new_2.txt') ## dic21
    lims = np.loadtxt(path[:-7] + 'lims_cc_new.txt') ## mar22,hs21

    idx=np.where(lims[:,3]> 1e-8)

    cohs_new=['par'+str(x).zfill(2)+'.txt' for x in idx[0] ]
    lims=lims[:nf]
    lims=np.array([x for x in lims if int(x[2]) != 0])
    n=0

    coordenadas=coordenadas[idx]
    #coord_uni=
    for lim, file in zip(lims, cohs_new):
        data=np.loadtxt(path+file)
        f=data[:,0]
        coh=data[:,1]
        #idx = np.where((f < lim[1]) & (f >= lim[0]))[0]
        #fx=f[idx]
        Cx=coh[idx]
        f1=lim[1]
        f2=lim[2]
        f3=lim[3]
        dist=lim[4]
        smooth_factor=lim[5]
        #smooth_factor
        sigma_d=1 ## varianza de los datos
        #sigma_D=np.logspace(-4,4,80)
        sigma_D=1 ## varianza de la suavizacion
        sigma_A=1 ## varianza de la solucion inicial
        p1=coordenadas[n][0:2]
        p2=coordenadas[n][2:4]

        cms, zcs, wzero, r, jn, c_est, A_est, rho_pre, fx, CCX, CCX_smooth,fxp,CCXp,CCXp_smooth = cc.get_dispersion_curve2(f, coh, dist,
                                                                                                      int(smooth_factor), f1,
                                                                                                      f2, n, coord_uni,
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

        sigmad = 1
        sigma = 1
        alpha = 1
        rhop, cp, Ap, H2, CM, err = cc.get_dispersion_curve_newton(fx2 * 2 * np.pi, dist, rho_pre2, CCX2_smooth, fxp,
                                                                   CCXp,
                                                                   CCXp_smooth, c_est2, A_est, sigmad, sigma,
                                                                   alpha, n, coord_uni, np.vstack((p1, p2)), wzero2,
                                                                   zcs2)
        # rhop, cp, Ap, H2, CM,err=cc.get_dispersion_curve_newton(fx*2*np.pi, dist,rho_pre, CCX_smooth, fxp,CCXp,CCXp_smooth, c_est, A_est, sigmad, sigma,
        #                                                    alpha,i, coord_uni, np.vstack((p1,p2)),wzero,zcs)

        lembdas = cp / fx2
        # rhop, cp, Ap, H2, CM,err = cc.get_dispersion_curve_newton(fx * 2 * np.pi, dist, rho_pre, CCX_smooth, fxp,
        #                                                       CCXp, CCXp_smooth, c_est, A_est,sigmad, sigma, alpha,
        #                                                        n, coord_uni, np.vstack((p1, p2)),wzero,zcs)
        print('dist', dist)
        n+=1
        cps.append(cp)
        wnews.append(fx2*2*np.pi)


    fig,ax=plt.subplots()
    for i in range(len(cps)):
        ax.plot(wnews[i]/(2*np.pi),cps[i],'k',linewidth=0.5)
        ax.text(wnews[i][-1]/(2*np.pi),cps[i][-1],'c'+str(i).zfill(2),fontsize=10)
    plt.show()
    plt.savefig('AJUSTES/ajuste_todas.png', format='png', dpi=300, bbox_inches='tight')

if step2d5:
    # read i-th coherence file
    i = 0
    f1 =6
    f2 =16

    cps = []
    fxs = []
    dist=np.sqrt((coordenadas[i][0]-coordenadas[i][2])**2+(coordenadas[i][1]-coordenadas[i][3])**2)
    data=np.loadtxt(path+cohs[i])
    f = data[:, 0]
    coh = data[:, 1]
    fnew=np.linspace(f[0],f[-1],len(f)*5)
    interp=interpolate.interp1d(f,coh,kind='cubic')
    coh2=interp(fnew)
    ## define proper frequency limits and smooth factor
    smooth_factor = 100

    CC_smooth=cc.smooth(coh,smooth_factor)
    plt.plot(f,coh)

    plt.plot(f,CC_smooth)

    p1 = coordenadas[i][0:2]
    p2=coordenadas[i][2:4]
    A_best,rho_best,rho_obs,err_best,ccombs_best,c_best,wx=cc.get_dispersion_curve22(f,coh,dist,smooth_factor,f1,f2,i, coord_uni, np.vstack((p1,p2)))
    plt.figure(33)
    plt.plot(wx/(2*np.pi),rho_obs)
    plt.plot(wx/(2*np.pi),rho_best)
    # cms,zcs,wzero,r,jn,c_est,A_est,rho_pre,fx,CCX,CCX_smooth,fxp,CCXp,CCXp_smooth=cc.get_dispersion_curve2(f, coh, dist, smooth_factor, f1, f2,
    #                                                                                   i, coord_uni, np.vstack((p1,p2)),
    #                                                                                   create_figures=False)

    vmin=50
    vmax=800
    zero_crossings = np.where(np.diff(np.sign(CCX_smooth)))[0]

    vv=np.arange(50,810,10)
    feval=jv(0,fx/(2*np.pi)*r/vv[0])
    sigmad=1
    sigma=1
    alpha=1
    rhop, cp, Ap, H2, CM,err=cc.get_dispersion_curve_newton(fx*2*np.pi, dist,rho_pre, CCX_smooth, fxp,CCXp,CCXp_smooth, c_est, A_est, sigmad, sigma,
                                                        alpha,i, coord_uni, np.vstack((p1,p2)),wzero,zcs)

    lembdas=cp/fx
    print('dist',dist)
    print('lambdamax',lembdas[0])
    print('lambdamin',lembdas[-1])
    print('lambda/r',lembdas[0]/dist)
    print('lambda/r',lembdas[-1]/dist)

