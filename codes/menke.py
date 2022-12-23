import numpy as np
from scipy.special import jv
import matplotlib.pyplot as plt

## Set up frequencies
N=101
Tmin=10
Tmax=100
wmin=2*np.pi/Tmax
wmax=2*np.pi/Tmin
iw=np.arange(0,N,1)
w=np.linspace(wmin,wmax, (N))
T=1/w
## true phase velocity
r=200
c1=4
c2=3
c_true=np.linspace(c1,c2,N)
c_true+=0.1*np.sin(np.pi*(w-wmin)/(wmax-wmin))
A=1.0

## this version uses lookup into a precomputed version of J0
## argument of J array
Nj=1001
jargmin=0.5*wmin*r/c1
jargmax=2.0*wmax*r/c2
Lj=jargmax-jargmin
jarg=np.linspace(jargmin,jargmax,Nj)
nu=0
J0=jv(nu,jarg)

## plot true phase velocity
plt.figure(1)
plt.plot(w,c_true,'k',linewidth=2)
plt.xlabel('w')
plt.ylabel('c')

# true correlation function

rho=np.zeros(N)
nu=0
rho_true=A*jv(nu,w*r/c_true)

# make correlated noise
s_rho=0.1
uncorrelated_noise=np.random.normal(0.0,s_rho,N)
L=10
correlated_noise=np.convolve(uncorrelated_noise,np.ones(L)/L)
correlated_noise=correlated_noise[:N]
sc=np.std(correlated_noise)

## synthethic observed correlogra, equals true plus noise

rho_obs=rho_true+correlated_noise
plt.figure(2)
plt.plot(w,rho_true,'k',linewidth=2)
plt.plot(w,rho_obs,'r',linewidth=2)
plt.xlabel('w')
plt.ylabel('rho')

## grid search over phase velocities, with least squares fit of A
## frequency nodes
Ngrid=4
Ngrid=3
## w-index of nodes
igrid=np.floor((N-1)*np.arange(0,Ngrid)/(Ngrid-1)).astype(int)
## values of nodes

wgrid=w[igrid]
## bounds
clow=2
chigh=5
## number of trial phase velocities between bounds
Nc=40
Nc=5
cgrid=np.linspace(clow,chigh,Nc)
## number of trial models
Nsearch=Nc**Ngrid
## indices of c at nodes
p=np.zeros(Ngrid).astype(int)
## for lstsq
denom=rho_obs@rho_obs
Dj=(Nj-1)/Lj
cnodes=np.zeros(Ngrid)

for it in range(Nsearch):
    #print(it)
    ## build phase velocity of nodes
    print('lol',cgrid[p])
    cnodes=cgrid[p]
    for ig in range(Ngrid):
        cnodes[ig]=cgrid[p[ig]]
        print('q',cgrid[p[ig]])
    #print(p)
    ## interpolated phase velocities between nodes
    ## a lot of effort here to avoid using interp1
    ## c_trial= interp1(wgrid,cnodes,w,'linear','extrap') in matlab
    c_trial=np.zeros(N)
    for iv in range(Ngrid-1):
        igl=igrid[iv]
        igr=igrid[iv+1]
        Lg=igr-igl+1
        c_trial[igl:igr+1]=cnodes[iv]+(cnodes[iv+1]-cnodes[iv])*(iw[igl:igr+1]-igl)/Lg
    ## evaluate model by lookup into precomputed J0
    imyarg=np.floor(Dj*((w*r/c_trial)-jargmin)).astype(int)
    rho_pre=J0[imyarg]

    ## error
    # determine A by least squares
    A_est=(rho_pre@rho_obs)/denom
    rho_pre*=A_est
    Drho=rho_obs-rho_pre
    E=(rho_obs-rho_pre)@(rho_obs-rho_pre)
    #print(E)
    #E=Drho@Drho

    # save the model if it is the current best
    if it == 0:
        E_best=E
        c_best=c_trial
        p_best=p
        A_best=A_est
        rho_best=rho_pre
        cnodes_best=cnodes
    if E < E_best:
        E_best = E
        c_best = c_trial.copy()
        p_best = p.copy()
        A_best = A_est
        rho_best = rho_pre
        cnodes_best = cnodes.copy()
    ## prepare for next model
    for ig in range(Ngrid):
        #print(p[ig])
        carry=0
        p[ig]=p[ig]+1
        if p[ig]>=Nc:
            p[ig]=0
            carry=1
        if carry==0:
            break
