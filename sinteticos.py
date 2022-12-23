import itertools
import funciones_linda as fl
import os
from scipy.special import jv,jn_zeros
from scipy.optimize import curve_fit
from scipy import interpolate
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy.signal import detrend, correlate, correlation_lags, tukey, butter, sosfilt, \
    welch,csd,resample,coherence
import datetime
from functools import partial
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,bicgstab,cg,cgs,gmres,LinearOperator,aslinearoperator,gcrotmk
from collections import Counter
import cc_processing as cc
from itertools import combinations

N=100

fmin=0.01
fmax=0.1
wmin=2*np.pi*fmin
wmax=2*np.pi*fmax

w=np.linspace(wmin,wmax,N)

r=200

chigh=4
clow=3
c_true=chigh + (clow-chigh)*np.arange(0,N,1)/(N)
c_true+=0.1*np.sin(np.pi*(w-wmin)/(wmax-wmin))
A=1.4
NJ=1001
jargmin=0.5*wmin*r/chigh
jargmax=2.0*wmax*r/clow
LJ=jargmax-jargmin
DJ=(NJ-1)/LJ
jarg = np.linspace(jargmin, jargmax, NJ)
nu = 0
J0 = jv(nu, jarg)
# plt.figure(1)
# plt.plot(w,c_true)
# plt.xlabel('w')
# plt.ylabel('c')

rho_true=A*jv(0,w*r/c_true)

#make correlated noise
s_rho=0.05
noise=np.random.normal(0,s_rho,len(rho_true))
rho_obs=rho_true+noise

# plt.figure(2)
# plt.plot(w,rho_true)
# plt.plot(w,rho_obs)

NN = 3
## numero de frecuencias
NC = 40
N = len(w)
cgrid = np.linspace(chigh, clow, NC)

cbounds=(chigh,clow)
cnoodles = np.array(list(combinations(cgrid, 3)))
A_best, rho_best, err_best, cnodes_best, c_best, wgrid = cc.grid_search_cw4(N, w, rho_obs, r, NN, NC, cbounds,
                                                                         cnoodles, True)

## varianza datos
sigma=1
## varianza minimizacion
small=1000
## varianza laplaciano
smooth=10
#c0=c_best
p=np.polyfit(w,c_best,1)
c0=w*p[0]+p[1]
popt,pcov=curve_fit(cc.tanhfit,w,c_best)


plt.figure(12)
plt.plot(w,c0)
plt.plot(w,c_best)


A0=A_best
rho0=A0*jv(0,w*r/c0)
Drho=rho_obs-rho0
E=Drho@Drho
E0=E
cp=c0
Ap=A0
rhop=rho0
Ntop=N+1
Nbot=N-2

#std datos
sigmad=1
# std smooth
alpha= 1
#std prior
sigma=1000
H = np.zeros((Ntop + Nbot, N + 1))
H1 = (sigma**-1) * np.eye(Ntop)
H2 = np.zeros((N, N))
def wlstsq(m,G,H):
    return G.T @ G @ m + H.T @ H @ m

f = lambda x: wlstsq(x, G, H)

for i in range(1, len(H2) - 1):
    H2[i, i - 1] = 1 # *smooth
    H2[i, i] = -2 # *smooth
    H2[i, i + 1] = 1 #*smooth


H2=alpha*H2

H[:N + 1, :N + 1] = H1
H3 = H2[1:-1]
H[N + 1:, :-1] = H3
H = scsp.csr_matrix(H)

h = np.zeros((Ntop + Nbot))
HH = H.T @ h

niter = 125
#H = scsp.csr_matrix(H)
G1 = np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2))
Cps=[]
Cps.append(cp)
Aps=[]
Aps.append(Ap)

#H=alpha*

for iter in range(niter):
        #print(iter)
        # mydiag=np.diag(Ap*jv(1,w*r/cp)*w*r*cp**(-2))
        # G1=scsp.csr_matrix((N,N+1))
        tic = time.time()

        G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
        G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
        G = sigmad**(-1)*scsp.hstack((G1, G2)).tocsr()
        ## delta rho, data vector
        Drho = rho_obs - rhop
        # E=Drho@Drho
        toc = time.time()
        #print('a Done in {:.4f} seconds'.format(toc - tic))

        F=G.T@Drho+HH

        tic = time.time()
        L = LinearOperator((G.shape[1], G.shape[1]), matvec=f, rmatvec=f)
        Dm = cg(L, F, tol=1e-05, maxiter=4 * N)

        #Aps.append(Dm[0][N])
        toc = time.time()
        #print('b Done in {:.4f} seconds'.format(toc - tic))

        tic = time.time()
        cpold=cp
        cp = cp + Dm[0][0:N]
        dcp=cp-cpold
        Cps.append(cp)
        Apold=Ap
        Ap = Ap + Dm[0][N]
        dAp=Ap-Apold
        ratio=100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold))
        Aps.append(Ap)
        print('deltas',np.sum(np.abs(dcp)),np.sum(np.abs(dAp)))
        print('ratio',100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold)))
        rhop = Ap * jv(0, w * r / cp)
        if ratio < 0.01:
            #G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
            #G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
            #G = sigma * scsp.hstack((G1, G2)).tocsr()
            CM=scsp.linalg.inv(scsp.vstack((G,H)).T@scsp.vstack((G,H))).todense()
            break
        toc = time.time()
        #print('c Done in {:.4f} seconds'.format(toc - tic))

        #print('c')
        #print(np.linalg.norm(Drho))
plt.figure(5)
plt.plot(w,cp,'g')
plt.plot(w,c_true,'r')
plt.plot(w,c0,'k')

plt.figure(6)
plt.plot(w,rhop)
plt.plot(w,rho0)
plt.plot(w,rho_obs)