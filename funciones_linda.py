# Functions
# collection of functions used in noise correlation processing
import numpy as np
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

def repeticiones(combinaciones,mediciones):
    dicc={}
    for x,y in zip(combinaciones,mediciones):
            
        dicc.setdefault(x,[]).append(y)
        
    # for x,y in zip(combinaciones,mediciones):
    #     if x not in dicc.keys():
    #         print(x,y)
    #         dicc.update({x:y})

    return dicc



def flatten(l):
    return [item for sublist in l for item in sublist]

def data_resample(path,resamp):
    data = np.genfromtxt(path + 'SynchronizedZ.dat')
    data=resample(data,data.shape[0]//resamp)

    #lags=correlation_lags(data.shape[0],data.shape[0],mode='same')
    #data=detrend(data,axis=1)
    #soshp=butter(10,0.1,btype='highpass',output='sos',fs=128)
    #data=sosfilt(soshp,data,axis=1)
    np.savetxt(path+'SynchronizedZ_Resamp.dat',data)
    #return data
    #t = np.linspace(-data.shape[0] / fs, data.shape[0] / fs, len(lags))


def fix_correlations_new(pathf,type,fs):
    data=np.genfromtxt(pathf)
    ND = [x for x in range(data.shape[1])]
    combs = list(itertools.combinations(ND, 2))

    if type==1:
        data1 =data[fs * 3:-1, 0]
        data2=data[:-1-fs*3,1]
        data3=data[:-1-fs*3,2]
        data4=data[:-1-fs*3,3]
        data5=data[:-1-fs*3,4]
        data_fixed=np.array([data1,data2,data3,data4,data5]).T
        np.savetxt(path + 'SynchronizedZ_Resamp_Fixed.dat', data_fixed)
def get_window(N, alpha=0.2):

    window = np.ones(N)
    x = np.linspace(-1., 1., N)
    ind1 = (abs(x) > 1 - alpha) * (x < 0)
    ind2 = (abs(x) > 1 - alpha) * (x > 0)
    window[ind1] = 0.5 * (1 - np.cos(np.pi * (x[ind1] + 1) / alpha))
    window[ind2] = 0.5 * (1 - np.cos(np.pi * (x[ind2] - 1) / alpha))
    return window

def normalize(data,fs, norm_win=None, norm_method="ramn"):

    if norm_method=='ramn':
        # recommended: half the  longest period

        lwin = fs * norm_win
        st = 0  # starting point
        N = lwin  # ending point

        while N < len(data):
            win = data[st:N]

            w = np.mean(np.abs(win)) / (2. * lwin + 1)

            # weight center of window
            print(st+lwin/2,st,lwin)
            #print(data.shape)
            data[int(st + lwin / 2)] /= w

            # shift window
            st += 1
            N += 1

        # taper edges
        taper = get_window(len(data))
        data *= taper

    if norm_method=='1bit':
        data=np.sign(data)
    return data
def raw_coherence(data,fs):
    N=data.shape[1]
    ND = [x for x in range(N)]
    combs = list(itertools.combinations(ND, 2))
    C=[]
    fig,ax=plt.subplots(len(combs)//2,2)
    for n,x in enumerate(combs):
        print('combs',x)
        f,Cij=coherence(data[:,x[0]],data[:,x[1]],fs=fs,nperseg=fs*4,noverlap=fs*2,detrend='linear')
        C.append(Cij)
        if n < len(combs)//2:
            ax[n,0].plot(f[:300],Cij[:300])
        else:
            ax[n-5,1].plot(f[:300],Cij[:300])

    #f,Cxy=coherence(x,y,fs=float(fs),nperseg=float(fs*4),noverlap=float(fs*2),detrend='linear')
    return f,C
    # tsecs=len(data)/fs
    # ## seconds of segments
    # length=len(data1)
    # nsegments=int(length/fs/nsecs)
    # data1=data1[:int(fs*nsecs*nsegments)]
    # data2=data2[:int(fs*nsecs*nsegments)]
    # print(fs,length,nsegments,nsecs)
    # #return data1,data2
    # data1_s=np.split(data1,nsegments)
    # data2_s=np.split(data2,nsegments)
    # #print()
    # R=[]
    # i=0
    # W=tukey(len(data1_s[0]),alpha=0.1)


    # for x,y in zip(data1_s,data2_s):
    #     x=detrend(x)
    #     y=detrend(y)
    #     rij=correlate(x,y)
    #     print(rij.shape)
    #     R.append(rij)
    #     i+=1


def correlations(data,fs,nsecs,norm=True):

    data1=data[:,0]
    data2=data[:,1]

    #nmin = int(len(data1) / fs / 2)
    ## total seconds on data
    tsecs=len(data)/fs
    ## seconds of segments
    length=len(data1)
    nsegments=int(length/fs/nsecs)
    data1=data1[:int(fs*nsecs*nsegments)]
    data2=data2[:int(fs*nsecs*nsegments)]
    print(fs,length,nsegments,nsecs)
    #return data1,data2
    data1_s=np.split(data1,nsegments)
    data2_s=np.split(data2,nsegments)
    #print()
    R=[]
    i=0
    W=tukey(len(data1_s[0]),alpha=0.1)


    for x,y in zip(data1_s,data2_s):
        x=detrend(x)
        y=detrend(y)
        if norm:
            x = normalize(x, 512,norm_win=2, norm_method='ramn')
            y = normalize(y, 512, norm_win=2,norm_method='ramn')

        #x = normalize(x, 512, norm_method='1bit')
        #y = normalize(y, 512, norm_method='1bit')
        #x=W*x
        #y=W*y

        rij=correlate(x,y)
        print(rij.shape)
        R.append(rij)
        i+=1
    Rsum = np.cumsum(R, axis=0)
    return Rsum
    C=[]
    for x,y in zip(data1_s,data2_s):
        x=detrend(x)
        y=detrend(y)

        #x = normalize(x, 512, norm_method='1bit')
        #y = normalize(y, 512, norm_method='1bit')
        x=W*x
        y=W*y
        num=fft(x)*np.conjugate(fft(y))/len(x)
        den=np.sqrt(fft(x)*np.conjugate(fft(x)))*np.sqrt(fft(y)*np.conjugate(fft(y)))*np.sqrt(len(x)**-2)
        #return num,den
        #num=fft(x)*np.conjugate(fft(y))
        #den=np.sqrt(fft(x)*np.conjugate(fft(x)))*np.sqrt(fft(y)*np.conjugate(fft(y)))
        cij=num/den
        #num=np.real((fft(x)*np.conjugate(fft(y))))
        #den=num*den
        #den=np.abs(fft(x))*np.abs(fft(y))
        #den=np.sqrt(((1 / len(x) ** 2)*fft(x)*fft(x)*fft(y)*fft(y)))
        #cij=num/1
        C.append(cij)
        #print(rij.shape)
        #R.append(rij)
        i+=1
    Rsum=np.cumsum(R,axis=0)
    Csum=np.cumsum(C,axis=0)
    return Rsum,Csum
def normalize(tr, clip_factor=6, clip_weight=10, norm_win=None, norm_method="1bit"): 
    print('qbob')
    if norm_method == 'clipping':
        lim = clip_factor * np.std(tr.data)
        tr.data[tr.data > lim] = lim
        tr.data[tr.data < -lim] = -lim

    elif norm_method == "clipping_iter":
        lim = clip_factor * np.std(np.abs(tr.data))
        
        # as long as still values left above the waterlevel, clip_weight
        while tr.data[np.abs(tr.data) > lim] != []:
            tr.data[tr.data > lim] /= clip_weight
            tr.data[tr.data < -lim] /= clip_weight

    elif norm_method == 'ramn':
        lwin = int(tr.stats.sampling_rate * norm_win)
        st = 0                                               # starting point
        N = lwin                                             # ending point

        while N < tr.stats.npts:
            #print('a',st,N)
            win = tr.data[st:N]

            w = np.mean(np.abs(win)) / (2. * lwin + 1)
            
            # weight center of window
            #print(st,lwin)
            tr.data[st + lwin // 2] /= w

            # shift window
            st += 1
            N += 1

        # taper edges
        #taper = get_window(tr.stats.npts)
        #tr.data *= taper

    elif norm_method == "1bit":
        tr.data = np.sign(tr.data)
        tr.data = np.float32(tr.data)

    return tr


def get_window(N, alpha=0.2):

    window = np.ones(N)
    x = np.linspace(-1., 1., N)
    ind1 = (abs(x) > 1 - alpha) * (x < 0)
    ind2 = (abs(x) > 1 - alpha) * (x > 0)
    window[ind1] = 0.5 * (1 - np.cos(np.pi * (x[ind1] + 1) / alpha))
    window[ind2] = 0.5 * (1 - np.cos(np.pi * (x[ind2] - 1) / alpha))
    return window


def whiten(tr, freqmin, freqmax):
    
    nsamp = tr.stats.sampling_rate
    
    n = len(tr.data)
    if n == 1:
        return tr
    else: 
        frange = float(freqmax) - float(freqmin)
        nsmo = int(np.fix(min(0.01, 0.5 * (frange)) * float(n) / nsamp))
        f = np.arange(n) * nsamp / (n - 1.)
        JJ = ((f > float(freqmin)) & (f<float(freqmax))).nonzero()[0]
            
        # signal FFT
        FFTs = np.fft.fft(tr.data)
        FFTsW = np.zeros(n) + 1j * np.zeros(n)

        # Apodization to the left with cos^2 (to smooth the discontinuities)
        smo1 = (np.cos(np.linspace(np.pi / 2, np.pi, nsmo+1))**2)
        FFTsW[JJ[0]:JJ[0]+nsmo+1] = smo1 * np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0]+nsmo+1]))

        # boxcar
        FFTsW[JJ[0]+nsmo+1:JJ[-1]-nsmo] = np.ones(len(JJ) - 2 * (nsmo+1))\
        * np.exp(1j * np.angle(FFTs[JJ[0]+nsmo+1:JJ[-1]-nsmo]))

        # Apodization to the right with cos^2 (to smooth the discontinuities)
        smo2 = (np.cos(np.linspace(0., np.pi/2., nsmo+1.))**2.)
        espo = np.exp(1j * np.angle(FFTs[JJ[-1]-nsmo:JJ[-1]+1]))
        FFTsW[JJ[-1]-nsmo:JJ[-1]+1] = smo2 * espo

        whitedata = 2. * np.fft.ifft(FFTsW).real
        
        tr.data = np.require(whitedata, dtype="float32")

        return tr
# Pxy = csd(data1, data2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
#
# Pxx = welch(data1, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
# Pyy = welch(data2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
# Cxyp = (1 / len(Pxy[0])) * np.real(Pxy[1]) / np.sqrt(((1 / len(Pxy[0])) ** 2) * np.real(Pxx[1]) * np.real(Pyy[1]))
# #data_sp = np.split(Cxyp, 2)
# #Cxy = (data_sp[0] + data_sp[1][::-1]) / 2
def cross_coherence(st,fs,s,stations,norm=None,onesided=True):

    sig1=st.select(station=stations[0])[0]
    #sig1.decimate(4)
    sig1.filter('highpass',freq=0.9,corners=2)
    #sig1.filter('bandpass',freqmin=2.5,freqmax=40.0)

    sig2=st.select(station=stations[1])[0]
    sig2.filter('highpass',freq=0.9,corners=2)
    #sig2.decimate(4)

    #sig2.filter('bandpass',freqmin=2.5,freqmax=40.0)

    print('pre',type(sig1),type(sig2))

    size=int(fs*s)
    step=int(size/2)
    ## alpha=0.1 yields 5% taper at each end
    W=tukey(size,alpha=0.1)
    if norm == '1bit':
        sig1=np.sign(sig1)
        sig2=np.sign(sig2)
        print('1bit',type(sig1),type(sig2))

    if norm == 'ramn':
        print('ramn',type(sig1),type(sig2))

        sig1=normalize(sig1,norm_win=0.25,norm_method='ramn')
        sig2=normalize(sig2,norm_win=0.25,norm_method='ramn')
    #return sig1,sig2

    Pxy = csd(sig1, sig2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=onesided)
    Pxx = welch(sig1, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=onesided)
    Pyy = welch(sig2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=onesided)

    #Pxy = csd(sig1, sig2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
    #Pxx = welch(sig1, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
    #Pyy = welch(sig2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)
    Cxyp = ((1 / len(Pxy[0])) * np.real(Pxy[1]) / np.sqrt(((1 / len(Pxy[0])) ** 2) * np.real(Pxx[1]) * np.real(Pyy[1])))
    print(Cxyp.shape)
    f=Pxy[0]
    if onesided:
        Cxy=Cxyp
    else:
        f=f[:len(f)//2]
        data_sp = np.split(Cxyp, 2)
        Cxy = (data_sp[0] + data_sp[1][::-1]) / 2
    return f,Cxy,Cxyp

    #f=f[:len(f)//2]
    #return Cxyp
    #data_sp = np.split(Cxyp, 2)
    #Cxy = (data_sp[0] + data_sp[1][::-1]) / 2
    #f,Cxy=coherence(sig1,sig2,fs=fs,nperseg=fs*4,noverlap=fs*2,detrend='linear')



    #f,Cxy=coherence(x,y,fs=float(fs),nperseg=float(fs*4),noverlap=float(fs*2),detrend='linear')  
def correlateNoise2(st,stations,corrwin):

    timewin=st.select(station=stations[1])[0].stats.starttime
    #timewin = st.select(station=stations[1])[0].stats.starttime + corrwin

    # loop over timewindows 
    # stop 1 corrwin before the end to account for different stream lengths
    #while timewin < st.select(station=stations[0])[-1].stats.endtime - 1*corrwin:
    while timewin < st.select(station=stations[0])[-1].stats.endtime - 2*corrwin:

        #print(timewin)
        sig1 = st.select(station=stations[0]).slice(timewin, timewin+corrwin)[0].detrend('linear').taper(0.1,type='cosine')
        #sig1.merge(method=0, fill_value=0)
        sig2 = st.select(station=stations[1]).slice(timewin, timewin+corrwin)[0].detrend('linear').taper(0.1,type='cosine')
        #sig2.merge(method=0, fill_value=0)
        #return sig1,sig2 
        #print(sig1.data)
        xcorr = np.correlate(sig1.data, sig2.data, 'same')

        try: 
            # build array with all correlations
            corr = np.vstack((corr, xcorr))
        except: 
            # if corr doesn't exist yet
            corr = xcorr
            
        # shift timewindow by one correlation window length
        timewin += corrwin

        # stack the correlations; normalize
        stack = np.sum(corr, 0)
        stack = stack / float((np.abs(stack).max()))    
    print ("...done")

    return corr, stack
def correlateNoise(st, stations, corrwin):

    print ('correlating stations', (stations[0], stations[1]))

    # initialize sliding timewindow (length = corrwin) for correlation
    # start 1 corrwin after the start to account for different stream lengths
    timewin = st.select(station=stations[1])[0].stats.starttime + corrwin
    timewin=st.select(station=stations[1])[0].stats.starttime

    # loop over timewindows 
    # stop 1 corrwin before the end to account for different stream lengths
    #while timewin < st.select(station=stations[0])[-1].stats.endtime - 2*corrwin:
    while timewin < st.select(station=stations[0])[-1].stats.endtime - 1*corrwin:

        print(timewin)
        sig1 = st.select(station=stations[0]).slice(timewin, timewin+corrwin)
        sig1.merge(method=0, fill_value=0)
        sig2 = st.select(station=stations[1]).slice(timewin, timewin+corrwin)
        sig2.merge(method=0, fill_value=0)
        xcorr = np.correlate(sig1[0].data, sig2[0].data, 'same')

        try: 
            # build array with all correlations
            corr = np.vstack((corr, xcorr))
        except: 
            # if corr doesn't exist yet
            corr = xcorr
            
        # shift timewindow by one correlation window length
        timewin += corrwin

        # stack the correlations; normalize
        stack = np.sum(corr, 0)
        stack = stack / float((np.abs(stack).max()))    
    print ("...done")

    return corr, stack


def plotStack(st, stack, maxlag, figurename=None):

    # define the time vector for the correlation (length of corr = corrwin + 1)
    limit = (len(stack) / 2.) * st[0].stats.delta
    timevec = np.arange(-limit, limit, st[0].stats.delta)

    plt.plot(timevec, stack, 'k')
    stations = list(set([_i.stats.station for _i in st]))
    plt.title("Stacked correlation between %s and %s" % (stations[0], stations[1]))
    plt.xlim(-maxlag, maxlag)
    plt.xlabel('time [s]')

    if figurename is not None:
        fig.savefig(figurename, format="pdf")
    else:
        plt.show()
        
        
def plotXcorrEvent(st, stn, stack, maxlag, acausal=False, figurename=None):

    eventtime = UTCDateTime(1998,7,15,4,53,21,0)                 # event near MLAC

    # station locations
    latP, lonP = 35.41, -120.55                                  # station PHL
    latM, lonM = 37.63, -118.84                                  # station MLAC
    latE, lonE = 37.55, -118.809                                 # event 1998
    
    # calculate distance between stations
    dist = gps2DistAzimuth(latP, lonP, latM, lonM)[0]            # between PHL and MLAC
    distE = gps2DistAzimuth(latP, lonP, latE, lonE)[0]           # between event and PHL
                                                                 #
    # CROSSCORRELATION
    # reverse stack to plot acausal part (= negative times of correlation)
    if acausal:
        stack = stack[::-1]
    
    # find center of stack
    c = int(np.ceil(len(stack)/2.) + 1)
    
    #cut stack to maxlag
    stack = stack[c - maxlag * int(np.ceil(stn[0].stats.sampling_rate)) : c + maxlag * int(np.ceil(stn[0].stats.sampling_rate))]
    
    # find new center of stack
    c2 = int(np.ceil(len(stack)/2.) + 1)

    # define time vector for cross correlation
    limit = (len(stack) / 2.) * stn[0].stats.delta
    timevec = np.arange(-limit, limit, stn[0].stats.delta)
    # define timevector: dist / t
    timevecDist = dist / timevec
    
    # EVENT
    ste = st.copy()
    st_PHL_e = ste.select(station='PHL')
    
    # cut down event trace to 'maxlag' seconds
    dt = len(stack[c2:])/stn[0].stats.sampling_rate                  #xcorrlength
    st_PHL_e[0].trim(eventtime, eventtime + dt)
    
    # create time vector for event signal
    # extreme values:
    limit = st_PHL_e[0].stats.npts * st_PHL_e[0].stats.delta
    timevecSig = np.arange(0, limit, st_PHL_e[0].stats.delta)

    # PLOTTING
    fig = plt.figure(figsize=(12.0, 6.0))
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    # plot noise correlation
    ax1.plot(timevecDist[c2:], stack[c2:], 'k')
    ax1.set_title('Noise correlation between MLAC and PHL')

    # plot event near MLAC measured at PHL
    ax2.plot(distE/timevecSig, st_PHL_e[0].data / np.max(np.abs(st_PHL_e[0].data)), 'r')
    ax2.set_title('Event near MLAC observed at PHL')

    ax2.set_xlim((0, 8000))
    ax1.set_xlim((0, 8000))

    ax2.set_xlabel("group velocity [m/s]")
 
    if figurename is not None:
        fig.savefig(figurename, format="pdf")
    else:
        plt.show()