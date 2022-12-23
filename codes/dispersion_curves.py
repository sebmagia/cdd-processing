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
    path+='Hs21/STACKS/'
if Mar22:
    path+='Mar22/STACKS/'
if Sept22:
    path+='Sept22/STACKS/'

files=sorted(os.listdir(path))
nf=2
coords=np.loadtxt(path[:-7]+'coords_all.txt',max_rows=nf)
files=files[:nf]
#file=files[0]
fs=512
nsegs=1
dists=[]
#stream=obspy.read(path+file,format='MSEED')
for pos,file in zip(coords,files):
    dist= ((pos[0]-pos[2])**2 + (pos[1]-pos[3])**2)**0.5
    dists.append(dist)
    stream = obspy.read(path + file, format='MSEED')
    ND = [x.stats.station for x in stream]
    f,coh,cohp=fl.cross_coherence(stream,fs,nsegs,ND,norm='1bit')
    try:
        Cij=np.vstack((Cij,coh))
    except:
        Cij=coh


i=0
C=Cij[i,:]
idx = np.where((f < 40.7) & (f >= (2)))[0]
fx=f[idx]
Cx=C[idx]
w=2*np.pi*fx
wmin=w[0]
wmax=w[-1]
wnew=np.linspace(wmin,wmax,120)
#w+=1e-08
#wnew+=1e-08
interp=interpolate.interp1d(w,Cx,kind='cubic')
rho_obs=interp(wnew)
plt.plot(wnew/(2*np.pi),rho_obs)
print('distancia',dists[i])
bandera = True
if bandera:
    zero_crossings=np.where(np.diff(np.sign(rho_obs)))[0]
    wzero=wnew[zero_crossings]
    wpn=wzero[0::2]
    wnp=wzero[1::2]
    jn=jn_zeros(0,40)
    r=dists[0]

    m=0
    cms=[]
    for i in range(5):
        cm=wzero*r/(jn[2*i:2*i+len(wzero)])
        cms.append(cm)
    plt.figure(2)
    plt.plot(wpn/(2*np.pi),cms[0][0::2],'go-')
    plt.plot(wnp/(2*np.pi),cms[0][1::2],'ro-')
    plt.plot(wzero/(2*np.pi),cms[0],'ko-')
    cbounds=[85,225]

    clow=cbounds[0]
    chigh=cbounds[1]
    ## numero de nodos
    NN = 4
    ## numero de frecuencias
    NC = 40
    N=len(wnew)
    cgrid = np.linspace(chigh, clow, NC)
    from itertools import permutations,combinations,product
    #cnodes=np.array(list(product(cgrid,repeat=NN)))
    #print(len(list(cnodes)))
    #cnudes=np.array(list(permutations(cgrid,NN)))
    cnoodles=np.array(list(combinations(cgrid,NN)))
    A_best, rho_best, err_best, cnodes_best, c_best, wgrid = cc.grid_search_cw4(N, wnew, rho_obs, r, NN, NC, cbounds,
                                                                                cnoodles, True)
    aa=cc.tanhfit(wnew, 0.5, 0.5, 0.5)

    b1=0.01
    b2=np.inf
    bounds=([b1,-1/wnew[-1],b1],[b2,1/wnew[-1],b2])
    popt,pcov=curve_fit(cc.tanhfit,wnew[1:],c_best[1:])
    #popt2,pcov2=curve_fit(cc.tanh2fit,wnew[1:],c_best[1:])
    plt.figure(8)
    plt.plot(wnew/(2*np.pi),c_best)
    plt.plot(wnew/(2*np.pi),cc.tanhfit(wnew,*popt))

    cbest_int=cc.tanhfit(wnew,*popt)
    cbest_grid=c_best
    rho_0 = A_best * jv(0, wnew * r / cbest_int)
    rho_grid=rho_best
    sigmad=1
    alpha=1
    sigma=1

    rhop, cp, Ap,H2,CM = cc.newton_search(wnew,rho_obs,rho_0,rho_grid,cbest_int,cbest_grid,A_best,r,sigmad,sigma,alpha)

    #plt.plot(wnew[1:],cc.tanh2fit(wnew[1:],*popt2))


    #cbest_int=cc.tanhfit(wnew[1:],*popt)
    #small=0.01
    #smooth=50
    #rhop, cp, Ap,H2 = cc.newton_search(wnew[1:],cnew,rho_best,cbest_int[1:],A_best,r,small,smooth)


    # rhop, cp, Ap,H2 = cc.newton_search(wnew,cnew,rho_best,c_best,A_best,r,small,smooth)

    #A_best, rho_best, err_best, cnodes_best, c_best, wgrid= cc.grid_search_cw3(N, wmin, wmax, cnew, r, NN, NC, cbounds,cnoodles,True)
    #ipw= cc.grid_search_cw3(N, wmin, wmax, cnew, r, NN, NC, cbounds,cnoodles,True)


    # #wgrid,J0,rho_obs,jargmin,jargmax=cc.grid_search_cw3(N, wmin, wmax, cnew, r, NN, NC, cbounds,cnoodles,True)
    #
    # iw = np.arange(0, N, 1)
    # igrid = np.floor((N - 1) * np.arange(0, NN) / (NN - 1)).astype(int)
    # wgrid = wnew[igrid]
    # cgrid=np.linspace(clow,chigh,NC)
    # NJ = 1001
    # jargmin = 0.5 * wmin * r / chigh
    # jargmax = 2.0 * wmax * r / clow
    # LJ = jargmax - jargmin
    # DJ = (NJ - 1) / LJ
    # jarg = np.linspace(jargmin, jargmax, NJ)
    # nu = 0
    # J0 = jv(nu, jarg)
    # rho_obs=cnew
    # ip1=interpolate.interp1d(wgrid,cnoodles[0,:],'linear',fill_value='extrapolate')
    # ip2=interpolate.interp1d(wgrid,cnoodles,'linear',fill_value='extrapolate')
    # #
    # ipw1=ip1(wnew)
    # ipw2=ip2(wnew)
    # #
    # rho_pre1 = J0[np.floor(DJ * ((wnew * r / ipw1) - jargmin)).astype(int)]
    # rho_pre2 = J0[np.floor(DJ * ((wnew * r / ipw2) - jargmin)).astype(int)]
    #
    # A_est1 = (rho_pre1 @ rho_obs)/(rho_pre1@rho_pre1)
    # A_est2 = (rho_pre2 @ rho_obs)/(np.linalg.norm(rho_pre2,axis=1))**2
    # A_est2t=A_est2.reshape(len(A_est2),1)
    # rho_pre2=np.multiply(A_est2t, rho_pre2)
    # Drho1=rho_obs-rho_pre1
    # Drho2=rho_obs-rho_pre2
    #
    # err1=Drho1@Drho1
    # err2=np.linalg.norm(Drho2,axis=1)**2
    #
    # agmin=np.argmin(err2)
    # A_best2=A_est2[agmin]
    # rho_best2=rho_pre2[agmin]
    # err_best2=err2[agmin]
    # c_best2=ipw2[agmin]
    # cnodes_best2=cnoodles[agmin]
    # plt.plot(wnew,rho_best,'g')
    # plt.plot(wnew,rho_best2,'r')
    # plt.plot(wnew,rho_obs,'k')
    #/ (rho_pre2 @ rho_pre2)

    #norm2=rho_pre2@rho_pre2
    #rho_pre22=rho_pre2[:500,:]
    #dot22=rho_pre22@rho_pre22.T
    #rho_obss=rho_obs.reshape(len(rho_obs),1)
    #norm1=np.linalg.norm(rho_pre[0,:],axis=-1)
    #norm3=np.linalg.norm(rho_pre2)


# plt.waitforbuttonpress()

# w=2*np.pi*fx
# wmin=w[0]
# wmax=w[-1]
#
# #rs = RunningStats()
#

# wnew=np.linspace(wmin,wmax,10000)
# dw=wnew[1]-wnew[0]
# df=dw/(2*np.pi)
# cnew=interp(wnew)
# ZC=np.where(np.diff(np.sign(cnew)))[0]
# ZCC=np.hstack((np.array([0]),ZC))
# ZC=np.hstack((ZCC,np.array(len(cnew)-1)))
# #fzc = ZCC[0]
# thresh = 5.0

# for i in range(0,len(ZC)-2,2):
#     #wi=c
#     wi=wnew[ZC[i]:ZC[i+2]]
#
#     arri=cnew[ZC[i]:ZC[i+2]]
#     div=np.max(cnew)/np.max(np.abs(arri))
#     mean=np.mean(arri)
#     std=np.std(arri)
#     amp=np.max(np.abs(arri))
#     #if div >= 5.0:
#         #break
#     #print(np.std(arri))
#     print('len',len(arri))
#     print('ratio',np.max(cnew)/np.max(np.abs(arri)))
    #print('meanstd',mean,std)
    #print('rstd',np.abs(mean/std))

    #print('stdr',(np.std(arri)))
    #print(np.max(np.abs(arri))/(np.std(arri)))
# i=1
# plt.plot(wnew/(2*np.pi),cnew)
# plt.plot(wnew[ZC]/(2*np.pi),cnew[ZC],'go')
# plt.plot(wi/(2*np.pi),arri)
# plt.plot(wnew[ZC[i]:ZC[i+2]]/(2*np.pi),cnew[ZC[i]:ZC[i+2]])
# window = 5  # size of the window
# Aw = np.lib.stride_tricks.sliding_window_view(Cx, window)
# Avar = np.var(Aw, axis=-1)
# Astd = np.std(Aw, axis=-1)
# Amean = np.mean(Aw, axis=-1)

# means=[]
# variances=[]
# stds=[]
#
# for i in cnew:
#     rs.push(i)
#     means.append(rs.mean())
#     variances.append(rs.variance())
#     stds.append(rs.standard_deviation())


# plt.figure(1)
# plt.plot(f,C)
# plt.plot(fx,Cx,'go')
# plt.plot(wnew/(2*np.pi),cnew)

#print(len(list(cnudes)))
# # ns = np.arange(0, 30)
# func = partial(product, repeat=NN)
# g = func(cgrid)
# #cnodes=np.array(list(g))
# #cnodes = np.split(np.array(list(g)), 256)
# cnodes=np.array(list(g))
# iw = np.arange(0, N, 1)
# w = np.linspace(wmin, wmax, N)
# #return w
# #T = 1 / w
# ## w-index of nodes, where we search for phase velocities
# igrid = np.floor((N - 1) * np.arange(0, NN) / (NN - 1)).astype(int)
# wgrid=w[igrid]
# clow=cbounds[0]
# chigh=cbounds[1]
# cgrid=np.linspace(clow,chigh,NC)
#  ## number of trial models
# # Nsearch=NC**NN
# ## for lstsq
# rho_obs=cnew
# denom=rho_obs@rho_obs
#
# NJ=1001
# jargmin=0.5*wmin*r/chigh
# jargmax=2.0*wmax*r/clow
# LJ=jargmax-jargmin
# Dj=(NJ-1)/LJ
# jarg = np.linspace(jargmin, jargmax, NJ)
# nu = 0
# J0 = jv(nu, jarg)
#
# ip = interpolate.interp1d(wgrid,cnodes,'linear', fill_value='extrapolate')
# ipw=ip(w)
# rho_pre = J0[np.floor(Dj * ((w * r / ipw) - jargmin)).astype(int)]
# rho_obss=rho_obs.reshape(len(rho_obs),1)
# A_est = (rho_pre @ rho_obss) #/ (rho_pre @ rho_pre)
#
# tic = time.time()
# #np.savetxt(path+'Cx.txt',np.vstack((wnew,cnew)).T)
#
#
#WW= cc.grid_search_cw2(N, wmin, wmax, cnew, r, NN, NC, cbounds,cnodes,True)
#
# A_best, rho_best, err_best, cnodes_best, c_best, wgrid,w= cc.grid_search_cw2(N, wmin, wmax, cnew, r, NN, NC, cbounds,cnodes,
#                                                                                  True)

# A_best, rho_best, err_best, cnodes_best, c_best, wgrid,w= cc.grid_search_cw(N, wmin, wmax, Cx, r, NN, NC, cbounds,
#                                                                                  cnodes, True)
    #m+=1

# plt.plot(w/(2*np.pi),Cx)
# plt.plot(wnew/(2*np.pi),cnew)
# plt.plot(wnew[zero_crossings]/(2*np.pi),cnew[zero_crossings],'go')
# np.sign(cnew)
#
# #zb=np.array([2.4048,5.5201,8.6537,11.7915,14.9309,18.0710,21.2116,24.35247,27.4934,30.6346])
#
# zcs=wzero*r/zb
# plt.plot(wnew,cnew,'go')
# plt.plot(w,Cx)

#plt.plot(f[idx],C[idx],'go')
#plt.xlabel('Frequency')
#plt.ylabel('Amplitude')
# #plt.savefig('C1.png')
# plt.savefig('C1.png', format='png', dpi=300, bbox_inches='tight')

#stream2=obspy.read(path2,format='MSEED')
#stream3=obspy.read(path3,format='MSEED')

#ND = [x.stats.station for x in stream]
#f,coh,cohp=fl.cross_coherence(stream,fs,nsegs,ND,norm='1bit')
#f2,coh2,cohp=fl.cross_coherence(stream2,fs,nsegs,['E1','E2'],norm='1bit')
#f3,coh3,cohp=fl.cross_coherence(stream3,fs,nsegs,['E1','E2'],norm='1bit')

#plt.plot(f,coh,'k')
#plt.plot(f2,coh2,'r')
#plt.plot(f3,coh3,'b')

# data1=stream[0]
# data2=stream[1]
#
# t=np.linspace(0,len(data1)*data1.stats.delta,len(data1))
# plt.plot(t,data1)
# plt.plot(t,data2)

#stream2=obspy.read(path2,format='MSEED')
#stream3=obspy.read(path3,format='MSEED')

#t=np.linspace(0,len(stream[0])*stream[0].stats.delta,len(stream[0]))
#t1=np.linspace(0,len(stream2[0])*stream2[0].stats.delta,len(stream2[0]))
#t2=np.linspace(0,len(stream3[0])*stream3[0].stats.delta,len(stream3[0]))

#plt.plot(t1,2*stream2[0].data)
#plt.plot(t,2*stream[0].data)
