import datetime
import itertools
import os
from datetime import datetime
import multiprocessing
from scipy.special import jv,jn_zeros
from scipy import interpolate
import matplotlib.pyplot as plt
from pandas import DataFrame,read_csv
import numpy as np
from scipy.fft import fft, ifft, fftshift, fftfreq
from scipy.signal import detrend, correlate, correlation_lags, tukey, butter, sosfilt, \
    welch,csd
import datetime
from functools import partial
from itertools import product
import time
from datetime import datetime
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,LinearOperator,aslinearoperator

def preprocess_data(datas,fs,s,filter,combs):
    "yields the NCF"    
    "data1: first time series" \
    "data2: second time series" \
    "s:lenght of segments in seconds" \
    "filter:whether we apply a filter or not"
    #print('combo',combs)
    data1=datas[combs[0]]
    data2=datas[combs[1]]
    if filter:
            soshp = butter(10, 0.9, btype='highpass', output='sos', fs=512)
            data1 = sosfilt(soshp, data1)
            data2 = sosfilt(soshp, data2)

    size = int(fs * s)
    step = int(size / 2)
    W = tukey(size, alpha=0.1)
    #print('ad',fs,s,W.shape,data1.shape,data2.shape)
    Pxy = csd(data1, data2, fs=fs, window=W, noverlap=step, detrend='linear',return_onesided=False)

    Pxx = welch(data1, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)

    Pyy = welch(data2, fs=fs, window=W, noverlap=step, detrend='linear', return_onesided=False)

    Cxyp = (1 / len(Pxy[0])) * np.real(Pxy[1]) / np.sqrt(((1 / len(Pxy[0])) ** 2) * np.real(Pxx[1]) * np.real(Pyy[1]))
    data_sp = np.split(Cxyp, 2)

    Cxy = (data_sp[0] + data_sp[1][::-1]) / 2
    return Cxy,Pxy[0]



def grid_search_cw2(N,wmin,wmax,rho_obs,r,NN,NC,cbounds,ccombs,create_figures):
    "initial estimate of dispersion curve"    

    #wmin=2*np.pi/Tmax
    #wmax=2*np.pi/Tmin
    iw = np.arange(0, N, 1)
    w = np.linspace(wmin, wmax, N)
    #T = 1 / w
    ## w-index of nodes, where we search for phase velocities
    igrid = np.floor((N - 1) * np.arange(0, NN) / (NN - 1)).astype(int)
    wgrid=w[igrid]
    clow=cbounds[0]
    chigh=cbounds[1]
    cgrid=np.linspace(clow,chigh,NC)
    ## number of trial models
    Nsearch=NC**NN
    ## for lstsq
    denom=rho_obs@rho_obs

    NJ=1001
    jargmin=0.5*wmin*r/chigh
    jargmax=2.0*wmax*r/clow
    LJ=jargmax-jargmin
    Dj=(NJ-1)/LJ
    jarg = np.linspace(jargmin, jargmax, NJ)
    nu = 0
    J0 = jv(nu, jarg)

    for n,x in enumerate(ccombs):
        #print(n,len(ccombs))
        #print('as',len(wgrid),len(x))
        ip = interpolate.interp1d(wgrid, x, 'linear', fill_value='extrapolate')
        #rho_pre=jv(0,w*r/ip(w))
        rho_pre = J0[np.floor(Dj * ((w * r / ip(w)) - jargmin)).astype(int)]
        A_est = (rho_pre @ rho_obs.reshape(len(rho_obs), 1)) / denom
        rho_pre = A_est * rho_pre
        Drho = rho_obs - rho_pre
        err = np.sum(Drho * Drho, axis=1)
        idx = np.argmin(err)
        if n == 0:
            A_best = A_est[idx]
            rho_best = rho_pre[idx]
            err_best = err[idx]
            ccombs_best = ccombs[n][idx]
            c_best = ip(w)[idx]
        elif n > 0 and np.min(err) < err_best:
            # print('foo',np.min(err),cnodes_best)
            A_best = A_est[idx]
            rho_best = rho_pre[idx]
            err_best = err[idx]
            ccombs_best = ccombs[n][idx]
            c_best = ip(w)[idx]
    if create_figures:
        fig,ax=plt.subplots(2,1)
        #ax[0].plot(w, c_obs, 'k', linewidth=2, label='target phase vel')
        ax[0].plot(w, c_best, 'r', linewidth=2, label='est phase vel')
        ax[0].plot(wgrid, ccombs_best, 'go', linewidth=2)
        ax[0].set_xlabel('w')
        ax[0].set_ylabel('c(w)')
        ax[0].legend()

        ax[1].plot(w, rho_obs, 'k--', linewidth=2, label='target corr func')
        ax[1].plot(w, rho_best, 'r:', linewidth=2, label='estimated corr func')
        ax[1].legend()
        ax[1].set_xlabel('w')

    return A_best,rho_best,err_best,ccombs_best,c_best

def wlstsq(m,G,H):
    return G.T @ G @ m + H.T @ H @ m
    
def newton_search(w,rho_obs,rho0,c0,A0,r):
    " newton search from initial estimate of c(w)" 
    wmin=w[0]
    wmax=w[-1]
    N=len(w)
    Drho=rho_obs-rho0
    E=Drho@Drho
    cp=c0
    Ap=A0
    rhop=rho0
    Ntop = N + 1;
    Nbot = N - 2
    small = 0.01
    H1 = 0.01 * np.eye(Ntop)
    smooth = 10
    H2 = np.zeros((N, N))
    H = np.zeros((Ntop + Nbot, N + 1))
    for i in range(1, len(H2) - 1):
        H2[i, i - 1] = smooth
        H2[i, i] = -2 * smooth
        H2[i, i + 1] = smooth
    H[:N + 1, :N + 1] = H1
    H3 = H2[1:-1]
    H[N + 1:, :-1] = H3
    h = np.zeros((Ntop + Nbot))
    niter = 125
    H = scsp.csr_matrix(H)
    f = lambda x: wlstsq(x, G, H)
    HH = H.T @ h
    for iter in range(niter):
        #print(iter)
        # mydiag=np.diag(Ap*jv(1,w*r/cp)*w*r*cp**(-2))
        # G1=scsp.csr_matrix((N,N+1))
        G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
        G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
        G = scsp.hstack((G1, G2)).tocsr()
        ## delta rho, data vector
        Drho = rho_obs - rhop
        # E=Drho@Drho
        L = LinearOperator((G.shape[1], G.shape[1]), matvec=f, rmatvec=f)
        Dm = bicg(L, G.T @ Drho + HH, tol=1e-05, maxiter=4 * N)
        cp = cp + Dm[0][0:N]
        Ap = Ap + Dm[0][N]
        rhop = Ap * jv(0, w * r / cp)
