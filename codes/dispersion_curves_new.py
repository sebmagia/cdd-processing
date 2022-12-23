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
step2=False ## runn all coherences with correct freq. limits
step5=True## test each coherence to find proper freq. limits
step7=False
Dic21=True
Hs21=False
Mar22=False
Sept22=False

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

files=sorted(os.listdir(path))
coordenadas=np.loadtxt(path[:-7]+'coords_all.txt')
coord_split=np.split(coordenadas,2,axis=1)
coord_uni=np.unique(np.vstack((coord_split[0],coord_split[1])),axis=0)
cohs=[x for x in files if '.txt' in x]
files=[x for x in files if '.mseed' in x]

#stream=obspy.read(path+files[3])
#stream[0].stats.starttime=obspy.UTCDateTime(2021,12,21,15,00,00)
#stream[1].stats.starttime=obspy.UTCDateTime(2021,12,21,15,00,00)
#stream.plot()
#stream.plot(outfile='GIASIS_FIGURES/ruido03.png',size=(2000, 1000))

if step0:
    i=13
    dist=np.sqrt((coordenadas[i][0]-coordenadas[i][2])**2+(coordenadas[i][1]-coordenadas[i][3])**2)
    fs = 512
    nsegs = 60
    file=files[i]
    stream= obspy.read(path + file, format='MSEED')
    stream.decimate(1)
    ND = [x.stats.station for x in stream]
    #f1, coh1, cohp = fl.cross_coherence(stream, fs, 1, ND, norm='1bit')
    #f2, coh2, cohp = fl.cross_coherence(stream, fs, 5, ND, norm='1bit')
    #f3, coh3, cohp = fl.cross_coherence(stream, fs, 10, ND, norm='1bit')
    #f4, coh4, cohp = fl.cross_coherence(stream, fs, 15, ND, norm='1bit')
    #f5, coh5, cohp = fl.cross_coherence(stream, fs, 20, ND, norm='1bit')
    #f6, coh6, cohp = fl.cross_coherence(stream, fs, 25, ND, norm='1bit')
    f7, coh7, cohp = fl.cross_coherence(stream, fs, 30, ND, norm='1bit')
    f8, coh8, cohp = fl.cross_coherence(stream, fs, 45, ND, norm='1bit')
    f9, coh9, cohp = fl.cross_coherence(stream, fs, 60, ND, norm='1bit')

    #f5, coh4, cohp = fl.cross_coherence(stream, fs, 90, ND, norm='1bit')

    #plt.plot(f5[:5000],coh4[:5000])
    plt.plot(f9[:5000], coh9[:5000])
    plt.plot(f8[:5000], coh8[:5000])
    plt.plot(f7[:5000], coh7[:5000])
    # plt.plot(f6[:5000], coh6[:5000])
    # plt.plot(f5[:5000], coh5[:5000])
    # plt.plot(f4[:5000],coh4[:5000])
    # plt.plot(f3[:5000],coh3[:5000])
    # plt.plot(f2[:5000],coh2[:5000])
    # plt.plot(f1[:5000],coh1[:5000])

    print('dist',dist)

n=0
if step1:
    nf = 100

    fs = 512

    nsegs = 45
    files = files[:nf]
    n=0
    for file in files:
        dist = np.sqrt((coordenadas[n][0] - coordenadas[n][2]) ** 2 + (coordenadas[n][1] - coordenadas[n][3]) ** 2)

        stream = obspy.read(path + file, format='MSEED')
        stream.decimate(1)
        print(len(stream[0]))
        ND = [x.stats.station for x in stream]
        f,coh,cohp=fl.cross_coherence(stream,fs,nsegs,ND,norm='1bit')
        idx = np.where((f <= 64))[0]
        plt.plot(f[idx],coh[idx],'k',linewidth=2)
        plt.xlabel('Frequency [Hz] r='+str(dist))
        plt.ylabel('Amplitude')
        plt.savefig(path+file[:-5]+'png',format='png', dpi=300, bbox_inches='tight')
        plt.close()
        np.savetxt(path+file[:-5]+'txt',np.vstack((f,coh)).T)
        n+=1
if step5:

    i =1
    cps = []
    fxs = []
    nf = 100
    dist=np.sqrt((coordenadas[i][0]-coordenadas[i][2])**2+(coordenadas[i][1]-coordenadas[i][3])**2)
    # lims = np.loadtxt(path[:-7] + 'lims_cc.txt')
    # lim=lims[i]
    data=np.loadtxt(path+cohs[i])
    f = data[:, 0]
    coh = data[:, 1]
    #idx = np.where((f >=4) & (f <=30))[0]
    #fx=f[idx]
    #Cx=coh[idx]
    smooth_factor = 100
    f1=2
    f2=20
    p1 = coordenadas[i][0:2]
    p2=coordenadas[i][2:4]
    cms,zcs,wzero,r,jn,c_est,A_est,rho_pre,fx,CCX,CCX_smooth=cc.get_dispersion_curve2(f, coh, dist, smooth_factor, f1, f2, i, coord_uni, np.vstack((p1,p2)),create_figures=False)

    sigmad=1
    sigma=1
    alpha=1
    rhop, cp, Ap, H2, CM=cc.get_dispersion_curve_newton(fx*2*np.pi, dist,rho_pre, CCX_smooth, c_est, A_est, sigmad, sigma, alpha,i, coord_uni, np.vstack((p1,p2)))
    print('dist',dist)

    #jnneg=-1*jn[::-1]
    #jnnp=np.hstack((jnneg,jn))

    # for i in range(100):
    #     cm1 = wzero * r / (jn[2 * i:2 * i + len(wzero)])
    #     #cm2 = wzero * r / (jn[2 * i:2 * i + len(wzero)][::-1])
    #
    #     print(2 * i, 2 * i + len(wzero))
    #     cms.append(cm1)
    #     #cms.append(cm2)
    #     plt.plot(wzero / (2 * np.pi), cm1, 'o')
    #     #plt.plot(wzero / (2 * np.pi), cm2, 'o')


    # #csm1=smooth(Cx,6)
    # #csm2=smooth(Cx,12)
    # #csm3=smooth(Cx,18)
    # #csm4=smooth(Cx,32)
    # csm4=smooth(Cx,128)
    #
    # plt.plot(fx,Cx)
    # #plt.plot(fx,csm1)
    # #plt.plot(fx,csm2)
    # #plt.plot(fx,csm3)
    # plt.plot(fx,csm4)
    #
    # #dist=lim[2]
    # sigma_d=1 ## varianza de los datos
    # #sigma_D=np.logspace(-4,4,80)
    # sigma_D=0.01 ## varianza de la suavizacion
    # sigma_A=1 ## varianza de la solucion inicial
    # p1=coordenadas[i][0:2]
    # p2=coordenadas[i][2:4]

    #ccombs, cgrid= cc.get_dispersion_curve(fx, Cx, dist, sigma_d, sigma_D, sigma_A, i,coord_uni,np.vstack((p1,p2)),True)
    # rhop,cp,wnew,Ap,H2,CM,rho_grid=cc.get_dispersion_curve(fx, csm4, dist, sigma_d, sigma_D, sigma_A, i,coord_uni,np.vstack((p1,p2)),True)
    # print('distancia',dist)
    # zero_crossings = np.where(np.diff(np.sign(csm4)))[0]
    # wzero = fx[zero_crossings]*2*np.pi
    # jn = jn_zeros(0, 400)
    # cms=[]
    # cms2=[]
    # for i in range(5):
    #     cm = wzero * dist / (jn[2 * i:2 * i + len(wzero)])
    #     print(2 * i,2 * i + len(wzero))
    #     cms.append(cm)
    # plt.figure(1)
    # plt.plot(wzero / (2 * np.pi), cms[0], 'go')
    # plt.plot(wzero / (2 * np.pi), cms[1], 'ro')
    # plt.plot(wzero / (2 * np.pi), cms[2], 'ko')
    # plt.plot(wzero / (2 * np.pi), cms[3], 'bo')
    #
    # plt.plot(wnew/(2*np.pi),cp)
    # plt.figure(2)
    # plt.plot(fx,csm4)
    # plt.plot(wnew/(2*np.pi),rhop)
    # plt.plot(wzero / (2 * np.pi),csm4[zero_crossings],'ro')


    # popt,pcov=curve_fit(cc.tanhfit,wgrids,cnodes_best)
    # cbest_int=cc.tanhfit(wgrids,*popt)
    # cbest_curve=cc.tanhfit(wnew,*popt)
    # plt.figure(1)
    # plt.plot(wgrids, cnodes_best,'go')
    # plt.plot(wnew,cbest_curve,'ko')
    # plt.plot(wgrids, cbest_int,'ro')
    #
    # popt2, pcov2 = curve_fit(cc.tanhfit, wnew, c_best)
    # cbest_int2 = cc.tanhfit(wnew, *popt2)
    # #plt.figure(2)
    # #plt.plot(wnew, c_best, 'ko')
    # plt.plot(wnew, cbest_int2, 'mo')
    #rhop, cp, wnew, Ap, H2, CM, rho_grid = cc.get_dispersion_curve(fx, Cx, dist, sigma_d, sigma_D, sigma_A, i,coord_uni,np.vstack((p1,p2)),True)

    # idx = np.where((f < 26) & (f >= (1)))[0]
    # fx = f[idx]
    # Cxx = coh[idx]
    # w = 2 * np.pi * fx
    # wmin = w[0]
    # wmax = w[-1]
    # wnew = np.linspace(wmin, wmax, 10000)
    # # w+=1e-08
    # # wnew+=1e-08
    # interp = interpolate.interp1d(w, Cxx, kind='cubic')
    # rho_obs = interp(wnew)
    # plt.plot(w/(2*np.pi),Cxx,'r')
    # plt.plot(wnew / (2 * np.pi), rho_obs,'k')
    # print('distancia', dist)
    # zero_crossings = np.where(np.diff(np.sign(rho_obs)))[0]
    # wzero = wnew[zero_crossings]
    # jn = jn_zeros(0, 400)
    # cms=[]
    # for i in range(5):
    #     cm = wzero * dist / (jn[2 * i:2 * i + len(wzero)])
    #     cms.append(cm)
    # m = 0
    # #cms = []
    # # for i in range(5):
    # #     cm = wzero * r / (jn[2 * i:2 * i + len(wzero)])
    # #     cms.append(cm)
    # plt.figure(19)
    # plt.plot(wnew/(2*np.pi),rho_obs)
    # plt.figure(20)
    # #plt.plot(wpn / (2 * np.pi), cms[0][0::2], 'go-')
    # #plt.plot(wnp / (2 * np.pi), cms[0][1::2], 'ro-')
    # plt.plot(wzero / (2 * np.pi), cms[0], 'ko-')
#export/home/otros/senunez/data/MSEEDS/

if step7:
    cps=[]
    fxs=[]
    wnews=[]
    nf=100
    lims = np.loadtxt(path[:-7] + 'lims_cc_new_2.txt')
    idx=np.where(lims[:,4]> 1e-8)

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
        dist=lim[3]
        smooth_factor=lim[4]
        #smooth_factor
        sigma_d=1 ## varianza de los datos
        #sigma_D=np.logspace(-4,4,80)
        sigma_D=1 ## varianza de la suavizacion
        sigma_A=1 ## varianza de la solucion inicial
        p1=coordenadas[n][0:2]
        p2=coordenadas[n][2:4]
        cms, zcs, wzero, r, jn, c_est, A_est, rho_pre, fx, CCX, CCX_smooth = cc.get_dispersion_curve2(f, coh, dist,
                                                                                                      int(smooth_factor), f1,
                                                                                                      f2, n, coord_uni,
                                                                                                      np.vstack(
                                                                                                          (p1, p2)),
                                                                                                      create_figures=False)


        sigmad = 1
        sigma = 1
        alpha = 0.1
        rhop, cp, Ap, H2, CM = cc.get_dispersion_curve_newton(fx * 2 * np.pi, dist, rho_pre, CCX_smooth, c_est, A_est,
                                                              sigmad, sigma, alpha, n, coord_uni, np.vstack((p1, p2)),wzero,zcs)
        print('dist', dist)
        #rhop, cp, wnew, Ap, H2, CM, rho_grid = cc.get_dispersion_curve(fx, Cx, dist, sigma_d, sigma_D, sigma_A, n,coord_uni,np.vstack((p1,p2)),True)
        n+=1
        cps.append(cp)
        wnews.append(fx*2*np.pi)


    fig,ax=plt.subplots()
    for i in range(len(cps)):
        ax.plot(wnews[i]/(2*np.pi),cps[i],'k',linewidth=0.5)
        ax.text(wnews[i][-1]/(2*np.pi),cps[i][-1],'c'+str(i).zfill(2),fontsize=10)
    plt.show()
    plt.savefig('AJUSTES/ajuste_todas.png', format='png', dpi=300, bbox_inches='tight')

if step2:
    cps=[]
    fxs=[]
    wnews=[]
    nf=100
    lims = np.loadtxt(path[:-7] + 'lims_cc_new.txt')
    idx=np.where(lims[:,3]> 1e-8)

    cohs_new=['par'+str(x).zfill(2)+'.txt' for x in idx[0] ]
    lims=lims[:nf]
    lims=np.array([x for x in lims if int(x[3]) != 0])
    n=0

    coordenadas=coordenadas[idx]
    #coord_uni=
    for lim, file in zip(lims, cohs_new):
        data=np.loadtxt(path+file)
        f=data[:,0]
        coh=data[:,1]
        idx = np.where((f < lim[1]) & (f >= lim[0]))[0]
        fx=f[idx]
        Cx=coh[idx]
        dist=lim[2]
        sigma_d=1 ## varianza de los datos
        #sigma_D=np.logspace(-4,4,80)
        sigma_D=0.25 ## varianza de la suavizacion
        sigma_A=1 ## varianza de la solucion inicial
        p1=coordenadas[n][0:2]
        p2=coordenadas[n][2:4]
        rhop, cp, wnew, Ap, H2, CM, rho_grid = cc.get_dispersion_curve(fx, Cx, dist, sigma_d, sigma_D, sigma_A, n,coord_uni,np.vstack((p1,p2)),True)
        n+=1
        cps.append(cp)
        wnews.append(wnew/(2*np.pi))


    fig,ax=plt.subplots()
    for i in range(len(cps)):
        ax.plot(wnews[i],cps[i],'k',linewidth=0.5)
        ax.text(wnews[i][-1],cps[i][-1],'c'+str(i).zfill(2),fontsize=10)
    plt.show()
    plt.savefig('AJUSTES/ajuste_todas.png', format='png', dpi=300, bbox_inches='tight')
    #np.savetxt(path[:-7]+'frecuencias_SA.txt',(np.array(fxs)/(2*np.pi)).T)
    #np.savetxt(path[:-7]+'curvas_SA.txt',np.array(cps).T)

    # fig,ax=plt.subplots(1)
    # ax.plot(np.mean(cps,axis=0),'r',linewidth=2)
    # ax.plot(np.mean(cps,axis=0)+2*np.std(cps,axis=0),'r--',linewidth=2)
    # ax.plot(np.mean(cps,axis=0)-2*np.std(cps,axis=0),'r--',linewidth=2)
    # ax.set_xlabel('Frequency [Hz]')
    # ax.set_ylabel('$C_R(f)$')
    # plt.show()
    # plt.savefig('AJUSTES/ajuste_mean.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()
