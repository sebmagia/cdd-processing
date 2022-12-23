import numpy as np
from scipy.special import jv,jn_zeros
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,LinearOperator,aslinearoperator
import matplotlib.pyplot as plt
from scipy import interpolate
from itertools import combinations_with_replacement,permutations,product
#bicg(3)
import time
import TimeSeries as ts
start_time=time.time()
## Set up frequencies
N=101
Tmin=10
Tmax=100
#Tmin=0.1
#Tmax=10
wmin=2*np.pi/Tmax
wmax=2*np.pi/Tmin
#wmin=50
#wmax=200
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
Nj=10001
jargmin=0.5*wmin*r/c1
jargmax=2.0*wmax*r/c2
Lj=jargmax-jargmin
jarg=np.linspace(jargmin,jargmax,Nj)
nu=0
J0=jv(nu,jarg)

## plot true phase velocity
# plt.figure(1)
# plt.plot(w,c_true,'k',linewidth=2)
# plt.xlabel('w')
# plt.ylabel('c')
# plt.title('true phase velocity')
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
wmin=2*np.pi/Tmax
wmax=2*np.pi/Tmin
NI=1000000
ww=np.linspace(wmin,wmax,NI)

rho_obs=rho_true+correlated_noise
ip=interpolate.interp1d(w,rho_obs,'linear', fill_value='extrapolate')
rho_int=ip(ww)
zc = np.where(np.diff(np.sign(rho_int)))[0]
pos_zc=zc[0::2]
neg_zc=zc[1::2]
wpos=ww[pos_zc]
wneg=ww[neg_zc]
wtod=ww[zc]
plt.plot(ww,rho_int)
#plt.plot(ww[zc],rho_int[zc],'go')
plt.plot(wpos,rho_int[pos_zc],'ro')
plt.plot(wneg,rho_int[neg_zc],'mo')

k=0
#zeros=np.loadtxt('bessel_zeros.txt')
zeros=jn_zeros(0,2000)
zeros_pos=zeros[::2]
zeros_neg=zeros[1::2]
#zeros=np.hstack((-zeros_pos[::-1],zeros_pos))
cc=np.zeros((25,zc.size))
ccpos=np.zeros((13,6))
ccneg=np.zeros((13,6))

for m in range(25):
    print(2*m,2*m+12)
    cc[m,:]=wtod*200/(zeros[m:m+zc.size])
# for m in range(13):
#     print(2*m,2*m+6)
#     ccpos[m,:]=wpos*200/(zeros_pos[2*m:2*m+6])
#     ccneg[m,:]=wneg*200/(zeros_neg[2*m:2*m+6])

    #print(cc[0,:])
#cc1=wpos*200/(zeros[0::2])
#cc2=wneg*200/(zeros[1::2])
plt.figure(33)
plt.plot(w,c_true,'g--')
plt.plot(wtod,cc[0,:])



# xx=np.linspace(0,40,NI)
# yy=jv(0,xx)
# plt.plot(ww,xx)
A_best,rho_best,err_best,cnodes_best,c_best=ts.grid_search_cw(N,wmin,wmax,rho_obs,r,3,40,(2,5),create_figures=True)
print("--- %s seconds ---" % (time.time() - start_time))
## compute error
c0=c_best
A0=A_best[0]
rho0=A0*jv(nu,w*r/c0)
Drho0=rho_obs-rho0
E=Drho0@Drho0
# y = wr c^-1
# d/dy Jo(y) = - J1(y)
# so d(rho)/dc = d(rho)/dy  dy/dc
# = -J1(wr/c) * -wr c^-2
# = J1(wr/c) wr/c^2


#starting values of linearized inversion
cp = c0;
Ap = A0;
rhop = rho0;

# matrix H of prior constraints
# first half, smallness
# second half, smoothness
def wlstsq(m,G,H):
    return G.T@G@m + H.T@H@m


#SECONDDERIV = 1;  % second derivative smoothing on 1, first deviv on 0
## Ntop is applied over length(c)+length(A)=N+1 parameters
## Nbot is Laplacian Smoothing applied over N-1 c(w) parameters
Ntop = N+1;
Nbot=N-2
small=0.01
H1=0.01*np.eye(Ntop)
smooth=50
H2=np.zeros((N,N))
H=np.zeros((Ntop+Nbot,N+1))
for i in range(1,len(H2)-1):
    H2[i,i-1]=smooth
    H2[i,i]=-2*smooth
    H2[i,i+1]=smooth
H[:N+1,:N+1]=H1
H3=H2[1:-1]
H[N+1:,:-1]=H3
h=np.zeros((Ntop+Nbot))
niter=120
H=scsp.csr_matrix(H)
f = lambda x: wlstsq(x, G, H)
HH=H.T@h
for iter in range(niter):
    print(iter)
    #mydiag=np.diag(Ap*jv(1,w*r/cp)*w*r*cp**(-2))
    #G1=scsp.csr_matrix((N,N+1))
    G1=scsp.lil_matrix(np.diag(Ap*jv(1,w*r/cp)*w*r*cp**(-2)))
    G2=scsp.lil_matrix(jv(nu,w*r/c0)).transpose()
    G=scsp.hstack((G1,G2)).tocsr()
    ## delta rho, data vector
    Drho=rho_obs-rhop
    #E=Drho@Drho
    L=LinearOperator((G.shape[1],G.shape[1]), matvec=f, rmatvec=f)
    Dm=bicg(L,G.T@Drho+HH,tol=1e-05,maxiter=4*N)
    cp=cp+Dm[0][0:N]
    Ap=Ap+Dm[0][N]
    rhop=Ap*jv(nu,w*r/cp)
# #m=np.hstack((cp,A0))
# #xd=wlstsq(m,G,H)
print("--- %s seconds ---" % (time.time() - start_time))
#
# #L=LinearOperator((G.shape[1],G.shape[1]),matvec=lambda x: wlstsq(x,G,H))
#
# # Hsp=scsp.csr_matrix(H1)
# # H = scsp.csr_matrix((Ntop,Nbot))
#     #spalloc(Ntop+Nbot,N+1,5*N);
fig, ax = plt.subplots(2, 1)
# ax[0].plot(w, c_obs, 'k', linewidth=2, label='target phase vel')
ax[0].plot(w, cp, 'r', linewidth=2, label='est phase vel.')
ax[0].plot(w, c_true, 'g--', linewidth=2, label='true phase vel.')
ax[0].plot(wtod,cc[1,:],'k--',linewidth=2,label='zero crossings')

ax[0].set_xlabel('w')
ax[0].set_ylabel('c(w)')
ax[0].legend()

ax[1].plot(w, rho_obs, 'k--', linewidth=2, label='target corr func')
ax[1].plot(w, rhop, 'r:', linewidth=2, label='estimated corr func')
ax[1].legend()
ax[1].set_xlabel('w')
plt.savefig('synthetic.png', format='png', dpi=300, bbox_inches='tight')
plt.figure(2)
plt.plot(w,rho_true,'k',linewidth=2,label='true correl')
plt.plot(w,rho_obs,'r',linewidth=2,label='noisy correl')
plt.xlabel('w')
plt.ylabel('rho')
plt.legend()
#plt.close()
# ## grid search over phase velocities, with least squares fit of A
# ## frequency nodes
# #Ngrid=4
# Ngrid=3
# ## w-index of nodes
# igrid=np.floor((N-1)*np.arange(0,Ngrid)/(Ngrid-1)).astype(int)
# ## values of nodes
#
# wgrid=w[igrid]
# ## bounds
# clow=2
# chigh=5
# ## number of trial phase velocities between bounds
# Nc=40
# #Nc=5
# cgrid=np.linspace(clow,chigh,Nc)
# ## number of trial models
# Nsearch=Nc**Ngrid
# ## indices of c at nodes
# p=np.zeros(Ngrid).astype(int)
# ## for lstsq
# denom=rho_obs@rho_obs
# Dj=(Nj-1)/Lj
# cnodes=np.zeros(Ngrid)
# #perm = list(permutations([0, 1, 2,3,4],3))
# x=np.arange(0,40)
# lol=[np.array(p) for p in product(x,repeat=Ngrid)]
# ps=[]
# cnodes=np.array([[cgrid[x[y]] for y in range(Ngrid)] for x in lol])
# c_nodes_sp=np.split(cnodes,256)
# del cnodes
# # ip = interpolate.interp1d(wgrid, cnodes, 'linear', fill_value='extrapolate')
# # #c_trial=ip(w)
# # #imyarg = np.floor(Dj * ((w * r / c_trial) - jargmin)).astype(int)
# # #rho_pre = J0[imyarg]
# # rho_pre = (J0[np.floor(Dj * ((w * r / ip(w)) - jargmin)).astype(int)])
# # A_est = (rho_pre @ rho_obs.reshape(len(rho_obs),1)) / denom
# # rho_pre = A_est*rho_pre
# # Drho = rho_obs - rho_pre
# # err=np.sum(Drho*Drho,axis=1)
# # idx=np.argmin(err)
# # A_best=A_est[idx]
# # rho_best=rho_pre[idx]
# # err_best=err[idx]
# # cnodes_best=cnodes[idx]
# for n,x in enumerate(c_nodes_sp):
#
#     ip=interpolate.interp1d(wgrid,x, 'linear', fill_value='extrapolate')
#     rho_pre = J0[np.floor(Dj * ((w * r / ip(w)) - jargmin)).astype(int)]
#     A_est = (rho_pre @ rho_obs.reshape(len(rho_obs),1)) / denom
#     rho_pre = A_est*rho_pre
#     Drho = rho_obs - rho_pre
#     err=np.sum(Drho*Drho,axis=1)
#     idx=np.argmin(err)
#     if n==0:
#         A_best=A_est[idx]
#         rho_best = rho_pre[idx]
#         err_best = err[idx]
#         cnodes_best = c_nodes_sp[n][idx]
#         c_best=ip(w)[idx]
#     elif n>0 and np.min(err)< err_best:
#         #print('foo',np.min(err),cnodes_best)
#         A_best = A_est[idx]
#         rho_best = rho_pre[idx]
#         err_best = err[idx]
#         cnodes_best = c_nodes_sp[n][idx]
#         c_best=ip(w)[idx]
# ## plot estimated and true phase velocity curve
# plt.figure(3)
# plt.plot(w,c_true,'k',linewidth=2,label='target phase vel')
# plt.plot(w,c_best,'r',linewidth=2,label='est phase vel')
# plt.plot(wgrid,cnodes_best,'go',linewidth=2)
# plt.legend()
# ## plot observed an estimated correlation function
# plt.figure(4)
# plt.plot(w,rho_obs,'k--',linewidth=2,label='target corr func')
# plt.plot(w,rho_true,'r:',linewidth=2,label='estimated corr func')
# plt.legend()
#c_trial=ip(w)
#imyarg = np.floor(Dj * ((w * r / c_trial) - jargmin)).astype(int)
#rho_pre = J0[imyarg]
#c_trial=ip(w)
#imyarg = np.floor(Dj * ((w * r / c_trial) - jargmin)).astype(int)
#rho_pre=np.memmap('memmapped.dat',dtype=np.float32,mode='w+',shape=(len(cnodes),N))
#rho_pre = (J0[np.floor(Dj * ((w * r / ip(w)) - jargmin)).astype(int)])
#rho_pre=rho_prex.astype(np.float32)
#del rho_prex
# error
# determine A by least squares
# A_est = (rho_pre @ rho_obs.reshape(len(rho_obs),1)) / denom
# rho_pre = A_est*rho_pre
# Drho = rho_obs - rho_pre
# err=np.sum(Drho*Drho,axis=1)
# idx=np.argmin(err)
# A_best=A_est[idx]
# rho_best=rho_pre[idx]
# err_best=err[idx]
# cnodes_best=cnodes[idx]

#dr1=Drho[0]@Drho[idx]


#Drho2=np.sum((rho_obs-rho_pre)**2)
#E = (rho_obs - rho_pre) @ (rho_obs - rho_pre)

    # print(E)
    # E=Drho@Drho
#cmag=np.array([np.array([cgrid[x[0]],cgrid[x[1]],cgrid[x[2]]]) for x in lol ])
# for it in range(Nsearch):
#     #print(it)
#     ## build phase velocity of nodes
#     #print('lol',cgrid[p])
#     cnodes=cgrid[p]
#     #print(p)
#     ## interpolated phase velocities between nodes
#     ## a lot of effort here to avoid using interp1
#     ## c_trial= interp1(wgrid,cnodes,w,'linear','extrap') in matlab
#     c_trial=np.zeros(N)
#     ip=interpolate.interp1d(wgrid, cnodes, 'linear', fill_value='extrapolate')
#     c_trial=ip(w)
#     imyarg = np.floor(Dj * ((w * r / c_trial) - jargmin)).astype(int)
#     rho_pre = J0[imyarg]
#     ## error
#     # determine A by least squares
#     A_est = (rho_pre @ rho_obs) / denom
#     rho_pre *= A_est
#     Drho = rho_obs - rho_pre
#     E = (rho_obs - rho_pre) @ (rho_obs - rho_pre)
#     # print(E)
#     # E=Drho@Drho
#
#     # save the model if it is the current best
#     if it == 0:
#         E_best = E
#         c_best = c_trial
#         p_best = p
#         A_best = A_est
#         rho_best = rho_pre
#         cnodes_best = cnodes
#     if E < E_best:
#         E_best = E
#         c_best = c_trial.copy()
#         p_best = p.copy()
#         A_best = A_est
#         rho_best = rho_pre
#         cnodes_best = cnodes.copy()
#     ## prepare for next model
#     #print('q')
#     ps.append(p.tolist())
#     for ig in range(Ngrid):
#         #print('a',p[ig])
#         print('b',p)
#         carry = 0
#         p[ig] = p[ig] + 1
#         if p[ig] >= Nc:
#             p[ig] = 0
#             carry = 1
#         if carry == 0:
#             break

    # for iv in range(Ngrid-1):
    #     igl=igrid[iv]
    #     igr=igrid[iv+1]
    #     Lg=igr-igl+1
    #     c_trial[igl:igr+1]=cnodes[iv]+(cnodes[iv+1]-cnodes[iv])*(iw[igl:igr+1]-igl)/(Lg-1)
    #     print('duh',iv,c_trial)
    ## evaluate model by lookup into precomputed J0



