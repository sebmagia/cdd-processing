import datetime
import itertools
import obspy
import funciones_linda as fl
import os
from datetime import datetime
import multiprocessing
from scipy.special import jv,jn_zeros
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
## step_1, fixing data delays
MED='Rect_1_27'
nw='HS'
path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/'+MED+'/'
path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/Los_Presidentes_HS/' + MED + '/'
fs = 512

step_1=False
step_2=False
step_3=False
step_4=True

check_coh=True
if step_1:
    ## load data
    data=np.genfromtxt(path+'SynchronizedZ.dat')
    ## plot coherence to check delayed stations
    f,C=fl.raw_coherence(data,fs)

    print('coherence 1 done')
    #fix stations
    data1=data[fs*3:-1,0]
    data2=data[:-1-fs*3,1]
    data3=data[:-1-fs*3,2]
    data4=data[:-1-fs*3,3]
    data5=data[:-1-fs*3,4]
    # # #
    # data1 = data[:, 0]
    # data2 = data[:, 1]
    # data3 = data[:, 2]
    # data4 = data[:, 3]
    # data5 = data[:, 4]

    # solo si una estacion no midio
    #data_new=np.vstack((data1,data2,data3,data4)).T

    data_new=np.vstack((data1,data2,data3,data4,data5)).T

    if check_coh:
        f2,C2=fl.raw_coherence(data_new,fs)
        print('coherence 2 done')
    # ## save fixed
        np.savetxt(path+'SynchronizedZ_fix.dat',data_new)

## step 2, generate mseed
#step_2=False
if step_2:
    data=np.genfromtxt(path+'SynchronizedZ_fix.dat')
    data1=data[:,0]
    data2=data[:,1]
    data3=data[:,2]
    data4=data[:,3]
    data5=data[:,4]

    trace1=obspy.Trace(data=data1,header={'sampling_rate':fs,'network':nw,'station':'E1'})
    trace2=obspy.Trace(data=data2,header={'sampling_rate':fs,'network':nw,'station':'E2'})
    trace3=obspy.Trace(data=data3,header={'sampling_rate':fs,'network':nw,'station':'E3'})
    trace4=obspy.Trace(data=data4,header={'sampling_rate':fs,'network':nw,'station':'E4'})
    trace5=obspy.Trace(data=data5,header={'sampling_rate':fs,'network':nw,'station':'E5'})
    # solo si una estacion no mitio
    #stream=obspy.Stream(traces=[trace1,trace2,trace3,trace4])

    stream=obspy.Stream(traces=[trace1,trace2,trace3,trace4,trace5])
    stream.write(path+MED+'.mseed',format='MSEED')

## generate cross_correlations
#step_3=False
nsegs=2
if step_3:
    stream=obspy.read(path+MED+'.mseed',format='MSEED')
    ND = [x.stats.station for x in stream]
    combs = list(itertools.combinations(ND, 2))
    Rij=[]
    for x in combs:
        print(x)
        corr,stack=fl.correlateNoise2(stream,x,nsegs)
        try:
            Rij=np.vstack((Rij,corr))
        except:
            Rij=corr
## generate cross-coherences
#step_4=True
nsegs=1
if step_4:
    stream=obspy.read(path+MED+'.mseed',format='MSEED')
    coords=np.loadtxt(path+'coords_'+MED+'.txt')
    #tr1=stream[0]
    #tr2=stream[1]
    #trnew1=tr1.copy()
    #trnew2=tr1.copy()
    #trnew3=tr1.copy()
    #trnew4=tr1.copy()

    #trnew1=fl.normalize(trnew1,norm_win=2,norm_method='ramn')
    #trnew2=fl.normalize(trnew2,norm_win=4,norm_method='ramn')
    #trnew3=fl.normalize(trnew3,norm_win=8,norm_method='ramn')
    #trnew4=fl.normalize(trnew4,norm_win=60,norm_method='ramn')

    ND = [x.stats.station for x in stream]
    combs = list(itertools.combinations(ND, 2))
    combs_coords=list(itertools.combinations(coords, 2))
    fig,ax=plt.subplots(len(combs)//2,2,figsize=(10,8))
    m=70
    for n,x in enumerate(combs):
        print(x)
        dist=np.linalg.norm(combs_coords[n][0]-combs_coords[n][1])
        #f,coh,cohp=fl.cross_coherence(stream,fs,nsegs,x,norm='1bit')
        #x,y=fl.cross_coherence(stream, fs, nsegs, x, norm='ramn')
        #f2,cxy2,cxyp2=fl.cross_coherence(stream,fs,nsegs,x,norm='ramn')
        #f3,cxy3,cxyp3=fl.cross_coherence(stream,fs,nsegs,x,norm=None)

        #plt.plot(f[:300],cxy[:300])
        #plt.plot(f2[:300],cxy2[:300])
        #plt.plot(f3[:300],cxy3[:300])

        f,coh,cohp=fl.cross_coherence(stream,fs,nsegs,x,norm='1bit')
        try:
            Cij=np.vstack((Cij,coh))
        except:
            Cij=coh

        if n < len(combs) // 2:
            ax[n, 0].plot(f[:m], coh[:m],linewidth='2')
            ax[n,0].text(30, 0.25, 'Intersta_Dist:' +str(np.round(dist,2))  + '[m]',fontsize='8')
            #ax[n,0].set_xlabel('Frequency')
        else:
            ax[n - 5, 1].plot(f[:m], coh[:m],linewidth='2')
            ax[n - 5, 1].text(30, 0.25, 'Intersta_Dist:' +str(np.round(dist,2))  + '[m]',fontsize='8')
            #ax[n-5,1].set_xlabel('Frequency')
    plt.savefig('Coh_1.png',format='png', dpi=300)
    ax[-1, 0].set_xlabel('Frequency (Hz)')
    ax[-1, 1].set_xlabel('Frequency (Hz)')
    ff=f.reshape((len(f),1))
    np.savetxt(path+'Cxy_'+MED+'.txt',np.hstack((ff,Cij.T)))

#fl.data_resample(path,1)
#path2='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/ARR1A/SynchronizedZ_Resamp.dat'
#fix_correlations_new(path2,1,512)
# path3='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/ARR1A/SynchronizedZ_Resamp_Fixed.dat'
# # data=np.genfromtxt(path3)
# # data1=data[:,0]
# # data2=data[:,1]
# # data3=data[:,2]
# # data4=data[:,3]
# # data5=data[:,4]
# # trace1=obspy.Trace(data=data1,header={'sampling_rate':512,'network':'A1','station':'E1'})
# # trace2=obspy.Trace(data=data2,header={'sampling_rate':512,'network':'A1','station':'E2'})
# # trace3=obspy.Trace(data=data3,header={'sampling_rate':512,'network':'A1','station':'E3'})
# # trace4=obspy.Trace(data=data4,header={'sampling_rate':512,'network':'A1','station':'E4'})
# # trace5=obspy.Trace(data=data5,header={'sampling_rate':512,'network':'A1','station':'E5'})
# path4='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/ARR1A/MED1A.mseed'
# cTAS=['E1','E2','E3','E4','E5']
#
# st1=obspy.read(path4)
# for tr in st1:
#     tr=fl.normalize(tr,norm_method='1bit')
#     #tr=whiten(tr,0.6,40)
# #tr1=fl.normalize(st1[0],norm_win=2,norm_method='ramn')
# corr,stack1=fl.correlateNoise(st1,['E1','E2'],3)
# #plt.plot(stack1,'k')
#
# st2=obspy.read(path4)
# for tr in st2:
#     tr=fl.normalize(tr,norm_method='ramn',norm_win=2)
#     #tr=whiten(tr,0.6,40)
#
# st3=obspy.read(path4)
#
# st4=obspy.read(path4)
# for tr in st4:
#     tr=fl.normalize(tr,norm_method='ramn',norm_win=20)
#     #tr=whiten(tr,0.6,40)
#
# #tr1=fl.normalize(st1[0],norm_win=2,norm_method='ramn')
#
#
# corr,stack1=fl.correlateNoise(st1,['E1','E2'],5)
# corr,stack2=fl.correlateNoise(st2,['E1','E2'],5)
# corr,stack3=fl.correlateNoise(st3,['E1','E2'],5)
# corr,stack4=fl.correlateNoise(st4,['E1','E2'],5)
#
# plt.plot(stack1,'k')
# plt.plot(stack2,'r')
# plt.plot(stack3,'g')
# plt.plot(stack4,'b')
#

#plt.plot(stack2,'r')
#stream=obspy.Stream(traces=[trace1,trace2,trace3,trace4,trace5])
#stream.write(path4,format='MSEED')
#xd=obspy.read('example.mseed')
#stream.plot()
#n=len(data)
#freqmin=0.1
#freqmax=25
#R=correlations(data,512,2,norm=False)
# for x in data.T:
#     wh=whiten(x,freqmin,freqmax)
#     try:
#         datawh=np.vstack((datawh,wh))
#     except:
#         datawh=wh
#datawh=datawh.T
#R2=correlations(data,512,2,norm=True)
#R3=correlations(data,512,2,norm=False)
#R2=correlations(datawh,512,2)
#from scipy.signal import coherence

#f,cxy=coherence(datawh[:,0],datawh[:,1],fs=512,nperseg=512*1,noverlap=128,detrend='linear')
# frange=float(freqmax)-float(freqmin)
# n=len(data[:,0])
# nsmo = int(np.fix(min(0.01, 0.5 * (frange)) * float(n) / nsamp))
# f = np.arange(n) * nsamp / (n - 1.)
#
# JJ = ((f > float(freqmin)) & (f < float(freqmax))).nonzero()[0]
# FFTs = np.fft.fft(data[:,0])
# FFTsW = np.zeros(n) + 1j * np.zeros(n)
#
# # Apodization to the left with cos^2 (to smooth the discontinuities)
# smo1 = (np.cos(np.linspace(np.pi / 2, np.pi, nsmo + 1)) ** 2)
# FFTsW[JJ[0]:JJ[0] + nsmo + 1] = smo1 * np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0] + nsmo + 1]))
#
# # boxcar
# FFTsW[JJ[0] + nsmo + 1:JJ[-1] - nsmo] = np.ones(len(JJ) - 2 * (nsmo + 1)) \
#                                         * np.exp(1j * np.angle(FFTs[JJ[0] + nsmo + 1:JJ[-1] - nsmo]))
#
# # Apodization to the right with cos^2 (to smooth the discontinuities)
# smo2 = (np.cos(np.linspace(0., np.pi / 2., nsmo + 1)) ** 2.)
# espo = np.exp(1j * np.angle(FFTs[JJ[-1] - nsmo:JJ[-1] + 1]))
# FFTsW[JJ[-1] - nsmo:JJ[-1] + 1] = smo2 * espo
#
# whitedata = 2. * np.fft.ifft(FFTsW).real
#
# datawh = np.require(whitedata, dtype="float32")
# plt.figure(1)
# plt.plot(fft(data[:,0]),'k')
# plt.figure(2)
# plt.plot(fft(datawh),'g')


#from scipy.signal import coherence
#data1=data[:,0]
#data2=data[:,3]

#f,cxy=coherence(data1,data2,fs=512,nperseg=512*1,noverlap=128,detrend='linear')
#f2,cxy2=csd(data1,data2,fs=512,nperseg=512*1,noverlap=128,detrend='linear')
#plt.figure(1)
#plt.plot(f,cxy)
# d1,d2=correlations(data,512)
#
# RS,CS=correlations(data,512,2)
# plt.figure(1)
# plt.plot(CS[-1])
# plt.figure(2)
# plt.plot(RS[-1])
# data1=data[:,0]
# data2=data[:,1]
# s=4
# fs=512
# size=int(fs*s)
# step=int(size/2)
# W=tukey(size,alpha=0.1)
# Pxy = csd(data1, data2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
#
# Pxx = welch(data1, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
# Pyy = welch(data2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
# Cxyp = (1 / len(Pxy[0])) * np.real(Pxy[1]) / np.sqrt(((1 / len(Pxy[0])) ** 2) * np.real(Pxx[1]) * np.real(Pyy[1]))
#data_norm=normalize(data[:,0],512,norm_method='1bit')
#plt.figure(1)
#plt.plot(Cxyp)
#plt.figure(2)
#plt.plot(ifft(Cxyp))


def fix_correlations(pathf,stage):
    print(pathf)
    data = np.genfromtxt(pathf)
    #return data
    data = np.genfromtxt(pathf)[:1280]
    datas=np.split(data,data.shape[1],axis=1)
    ND = [x for x in range(data.shape[1])]
    combs = list(itertools.combinations(ND, 2))

    R=[]
    tmaxes=[]
    #lags=correlation_lags(data.shape[0],data.shape[0],mode='full')
    lags=correlation_lags(data.shape[0],data.shape[0])
    #t = np.linspace(-data.shape[0] / 128, data.shape[0] / 128, len(lags))
    t = np.arange(-data.shape[0] / 128, data.shape[0] / 128, 1. / 128)[1:]

    idx=np.where((t>= -4) & (t<= 4) )
    # data to use for correlation
    #data_corr=data[np.where((t>= -4) & (t<= 4) )]
    #print(data_corr.shape)
    if stage==1:
        problematic_stations=[]
        for n, x in enumerate(combs):
            rij = correlate(data[:,x[0]], data[:, x[1]])
            t = np.arange(-data.shape[0] / 128, data.shape[0] / 128, 1. / 128)[1:]

            #rij=rij[idx]
            #rij = correlate(data_corr[:, x[0]], data_corr[:, x[1]])
            print(rij.shape,t.shape)
            R.append(rij)
            tmax = t[np.argmax(np.abs(rij))]
            if tmax > 2.5:
                problematic_stations.append(combs[n][0])
            if tmax < -2.5:
                problematic_stations.append(combs[n][1])

            tmaxes.append(tmax)

        datas_new=[]
        for n in range(len(datas)):
            if n in problematic_stations:
                datas_new.append(datas[n][:128*3-1])
            if n not in problematic_stations:
                datas_new.append(datas[n][:-1-128*3])

        datas_new=np.squeeze(np.array(datas_new)).T
        #np.savetxt(pathf + 'SynchronizedZ_Resamp_Fix.dat', datas_new)
        return tmaxes,R,t
    if stage==2:
        for n, x in enumerate(combs):
            rij = correlate(data_corr[:, x[0]], data_corr[:, x[1]])
            R.append(rij)
            tmax = t[np.argmax(np.abs(rij))]
            tmaxes.append(tmax)
        return tmaxes,R


# #path='/home/doctor/Doctor/Magister/Tesis/databases/Vs_Los_Presidentes/Lin_1_20/'
# path='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/ARR1A/'
# data_resample(path,1)
#
#
# s=10
# fs=512
# size=int(fs*s)
# step=int(size/2)
# W=tukey(size,alpha=0.1)
# data1=data[fs*3:-1,0]
# data2=data[:-1-fs*3,1]
# data3=data[:-1-fs*3,2]
# data4=data[:-1-fs*3,3]
# data5=data[:-1-fs*3,4]
#
# xd=correlate(data1[:100000],data2[:100000])
# t=np.arange(-100000/fs,100000/fs,1./fs)[1:]
# plt.plot(t,xd)
# Pxy = csd(data1, data2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
#
# Pxx = welch(data1, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
# Pyy = welch(data2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
# Cxyp = (1 / len(Pxy[0])) * np.real(Pxy[1]) / np.sqrt(((1 / len(Pxy[0])) ** 2) * np.real(Pxx[1]) * np.real(Pyy[1]))
# #data_sp = np.split(Cxyp, 2)
# #Cxy = (data_sp[0] + data_sp[1][::-1]) / 2
# plt.plot(Cxyp)
#
# #lags=correlation_lags(1000,1000,mode='same')
#
# nmin=int(len(data1)/fs/60)
# data1_s=np.split(data1[:fs*60*nmin],nmin)
# data2_s=np.split(data2[:fs*60*nmin],nmin)
# #data1_s=np.array_split(data1,473)
# #data2_s=np.array_split(data2,473)
# R=[]
# i=0
# W=tukey(len(data1_s[0]),alpha=0.1)
# #plt.plot(x)
# #plt.plot(W*x)
# for x,y in zip(data1_s,data2_s):
#     x=detrend(x)
#     y=detrend(x)
#     x=W*x
#     y=W*y
#     rij=correlate(x,y)
#     print(rij.shape)
#     R.append(rij)
#     i+=1
# C=[]
# for x,y in zip(data1_s,data2_s):
#     x=detrend(x)
#     y=detrend(x)
#     x=W*x
#     y=W*y
#     cij=(fft(x)*np.conjugate(fft(y)))
#     C.append(cij)
#     #print(rij.shape)
#     #R.append(rij)
#     i+=1
# Rsum=np.cumsum(R,axis=0)
# Csum=np.cumsum(C,axis=0)
# #plt.plot(Cxy)

########
#ff=welch(data[:,2],128)
#pathf='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/ARR1A/SynchronizedZ_Resamp.dat'
#data=fix_correlations(pathf,1)

#tmaxes1,R1,t=fix_correlations(pathf,1)
#pathf='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/ARR1A/SynchronizedZ_Resamp_Fix.dat'
#tmaxes2,R=fix_correlations(pathf,2)

# data=np.genfromtxt(path+'SynchronizedZ_resamp.dat')
# datas=np.split(data,5,axis=1)
# ND = [x for x in range(data.shape[1])]
# combs = list(itertools.combinations(ND, 2))
#
# R=[]
# #sigs2=[]
# #sigs3=[]
#
# tmaxes=[]
# lags=correlation_lags(data.shape[0],data.shape[0],mode='full')
# lags2=correlation_lags(data.shape[0],data.shape[0],mode='same')
#
# t = np.linspace(-data.shape[0] / 128, data.shape[0] / 128, len(lags))
# t=t[np.where((t>= -4) & (t<= 4) )]
# data_corr=data[np.where((t>= -4) & (t<= 4) )]
#
# #t2 = np.linspace(-data.shape[0] / 128, data.shape[0] / 128, len(lags2))
# #t2=t2[np.where((t2>= -4) & (t2<= 4) )]
#
# #t = np.arange(-data.shape[0] / 128, data.shape[0] / 128, 1. / 128)[1:]
# problematic_stations=[]
# for n,x in enumerate(combs):
#     rij=correlate(data_corr[:,x[0]],data_corr[:,x[1]],mode='same')
#     R.append(rij)
#     tmax=t[np.argmax(rij)]
#     if tmax > 2.5:
#         problematic_stations.append(combs[n][0])
#     if tmax < -2.5:
#         problematic_stations.append(combs[n][1])
#     # sig2=correlate(data[:,x[0]],data[:,x[1]],mode='valid')
#     # sigs2.append(sig2)
#
#     #sig3=correlate(data[:,x[0]],data[:,x[1]],mode='same')
#     #sigs3.append(sig3)
#     tmaxes.append(tmax)
#
# datas_new=[]
# for n in range(len(datas)):
#     if n in problematic_stations:
#         datas_new.append(datas[n][128*3:-1])
#     if n not in problematic_stations:
#         datas_new.append(datas[n][:-1-128*3])
#
# datas_new=np.squeeze(np.array(datas_new)).T
#
# tmaxes=[]
# lags=correlation_lags(datas_new.shape[0],datas_new.shape[0],mode='full')
# lags2=correlation_lags(datas_new.shape[0],datas_new.shape[0],mode='same')
#
# t = np.linspace(-datas_new.shape[0] / 128, datas_new.shape[0] / 128, len(lags))
# t=t[np.where((t>= -4) & (t<= 4) )]
# data_corr=datas_new[np.where((t>= -4) & (t<= 4) )]
# R=[]
# #sigs2=[]
# #sigs3=[]
#
# tmaxes=[]
# for n,x in enumerate(combs):
#     rij=correlate(data_corr[:,x[0]],data_corr[:,x[1]],mode='same')
#     R.append(rij)
#     tmax=t[np.argmax(rij)]
#     tmaxes.append(tmax)
# fig,ax=plt.subplots(5,2)
# ax[0,0].plot(t,sigs[0])
# ax[1,0].plot(t,sigs[1])
# ax[2,0].plot(t,sigs[2])
# ax[3,0].plot(t,sigs[3])
# ax[4,0].plot(t,sigs[4])
#
# ax[0,1].plot(t,sigs[5])
# ax[1,1].plot(t,sigs[6])
# ax[2,1].plot(t,sigs[7])
# ax[3,1].plot(t,sigs[8])
# ax[4,1].plot(t,sigs[9])


# fig,ax=plt.subplots(5,2)
# ax[0,0].plot(t2,sigs3[0])
# ax[1,0].plot(t2,sigs3[1])
# ax[2,0].plot(t2,sigs3[2])
# ax[3,0].plot(t2,sigs3[3])
# ax[4,0].plot(t2,sigs3[4])
#
# ax[0,1].plot(t2,sigs3[5])
# ax[1,1].plot(t2,sigs3[6])
# ax[2,1].plot(t2,sigs3[7])
# ax[3,1].plot(t2,sigs3[8])
# ax[4,1].plot(t2,sigs3[9])
# CC1=np.genfromtxt(path+'SynchronizedZ.dat')
# CC1_B=resample(CC1,CC1.shape[0]//2)
# CC1_C=resample(CC1,CC1.shape[0]//4)
# np.savetxt(path+'SynchronizedZ_resamp.dat',CC1_C)
#
#
# t1=np.arange(0,CC1.shape[0]/512,1/512)
# t2=np.arange(0,CC1_B.shape[0]/256,1/256)
# t3=np.arange(0,CC1_C.shape[0]/128,1/128)
#
# ff1=rfft(CC1[:,0])
# ff2=rfft(CC1_B[:,0])
# ff3=rfft(CC1_C[:,0])
#
# freq1=rfftfreq(CC1.shape[0],d=1./512)
# freq2=rfftfreq(CC1_B.shape[0],d=1./256)
# freq3=rfftfreq(CC1_C.shape[0],d=1./128)
# plt.figure(1)
# plt.plot(t1,CC1[:,0])
# plt.plot(t2,CC1_B[:,0])
# plt.plot(t3,CC1_C[:,0])
# plt.figure(2)
# plt.plot(freq1,ff1)
# plt.plot(freq2,ff2)
# plt.plot(freq3,ff3)

# #f= filter(os.path.isdir, os.listdir(path))
# carpetas = [path + x for x in os.listdir(path)]
# carpetas = sorted(list(filter(os.path.isdir,carpetas)))
# CC1=np.genfromtxt(carpetas[0]+'/EqualizedFile.dat',skip_header=100,encoding='ISO-8859-1')
# CC2=np.genfromtxt(carpetas[0]+'/EqualizedFile.dat',skip_header=100,encoding='ISO-8859-1')

#data=np.loadtxt()