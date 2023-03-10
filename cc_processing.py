import numpy as np
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
import numpy as np
from scipy.fft import rfft, fft, ifft, fftshift, fftfreq,rfftfreq
from scipy.signal import detrend, correlate, correlation_lags, tukey, butter, sosfilt, \
    welch,csd,resample,coherence
import datetime
from functools import partial
from itertools import permutations,combinations,product
import time
from matplotlib.gridspec import GridSpec
from datetime import datetime
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,cg,LinearOperator,aslinearoperator
from numba import jit
def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#def zero_crossings():


def get_dispersion_curve_newton22(w,r,rho_0,rho_obs,f_full,rhofull_noise,
                                rhofull_smooth,c0,A0,sigma_D,sigma_A,alpha,
                                ii,nfold,coord_uni,p12,wzero,zcs,CDINV,folder,med):
    wmin=w[0]
    wmax=w[-1]
    N=len(w)
    Drho=rho_obs-rho_0
    E=Drho@Drho
    cp=c0
    Ap=A0
    rhop=rho_0
    Ntop = N + 1
    Nbot = N - 2
    #std datos
    #sigmad=1
    # std smooth
    #alpha= 1
    #std prior
    #sigma=1000
    ## H Matrix with smallness and smoothness constraints
    H = np.zeros((Ntop + Nbot, N + 1))
    #H1 = sigma_A**(-1) * np.eye(Ntop)
    #smallness
    H1 = sigma_A * np.eye(Ntop)
    #H1=sigma
    #smooth = 1000
    # smoothness
    dw=w[1]-w[0]
    H2 = np.zeros((N, N))
    for i in range(1, len(H2) - 1):
        H2[i, i - 1] = 1
        H2[i, i] = -2 
        H2[i, i + 1] = 1
    #H2=sigma_D**(-1)*H2
    H2=alpha*H2
    H[:N + 1, :N + 1] = H1
    H3 = H2[1:-1]
    H[N + 1:, :-1] = H3
    h = np.zeros((Ntop + Nbot))
    niter = 25
    H = scsp.csr_matrix(H)
    #h=np.zeros(Nbot)
    #H=H2
    #print(H.shape,cp.shape,h.shape)
    f = lambda x: wlstsq(x, G, H)
    HH = H.T @ h
    CDINV=scsp.csr_matrix(CDINV)
    print('error_inicial', Drho@Drho)
    for iter in range(niter):
        #print('error',Drho@Drho)
        Drho0=Drho

        G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
        G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
        #G = sigma_d**(-1)*scsp.hstack((G1, G2)).tocsr()
        #G = sigma_D*scsp.hstack((G1, G2)).tocsr()
        G = CDINV@scsp.hstack((G1, G2)).tocsr()
        #print('xd',G.shape)
        ## delta rho, data vector
        #print(np.sum(Drho**2))
        F=G.T@Drho+HH

       
        L = LinearOperator((G.shape[1], G.shape[1]), matvec=f, rmatvec=f)
        Dm = cg(L, G.T @ Drho + HH, tol=1e-08, maxiter=4 * N)
        #Dm = cg(L, G.T@CD@G @ Drho + HH, tol=1e-08, maxiter=4 * N)

        cpold=cp
        cp = cp + Dm[0][0:N]
        #print('norma',np.linalg.norm(cp-cpold))
        #print(np.linalg.norm(cp-cpold)/())
        dcp=cp-cpold
        #print('norm',np.linalg.norm(dcp)/np.linalg.norm(cp))
        Ap = Ap + Dm[0][N]
        rhop = Ap * jv(0, w * r / cp)
        Drho = rho_obs - rhop
        #ratio=100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold))
        errpor=np.abs(100*(np.sum(Drho0)-np.sum(Drho))/(np.sum(Drho0)))

        #print('errpor',errpor)
        #print( np.sum(Drho),np.sum(Drho0))
        #if errpor <1:
            #break
        #print(( (np.sum((Drho)**2)-np.sum((Drho0)**2))/(np.sum((Drho0)**2)))*100)
        #print(ratio)
        #print('ratio',100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold)))

        #if ratio < 0.0001:
            
         #   break
    ## update G to generate posterior covariance matrix
    G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
    G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
    #G = (sigma_d)**-1 * scsp.hstack((G1, G2)).tocsr()
    G=CDINV * scsp.hstack((G1, G2)).tocsr()
    CM=scsp.linalg.inv(scsp.vstack((G,H)).T@scsp.vstack((G,H))).todense()
    err=np.abs(rho_obs-rhop)
    create_figures=True
    print('error_final', Drho@Drho)

    if create_figures:
        fig = plt.figure(tight_layout=True,figsize=(10,8))
        gs = GridSpec(2, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])


        ax1.plot(f_full, rhofull_noise, 'k', linewidth=2, label='tar noise')
        ax1.plot(f_full, rhofull_smooth, 'r', linewidth=2, label='tar smooth')
        ax1.plot(w/(2*np.pi), rho_0, 'g', linewidth=2, label='1st estimate')
        ax1.plot(w/(2*np.pi), rhop, 'y', linewidth=2, label='2nd estimate')
        #ax1.plot(wzero/(2*np.pi),rho_obs[zcs],'bo',label='ZC')
        ax1.set_ylabel('CC Amplitude')
        ax1.legend()
        #ax[0].set_xlim([min(w/(2*np.pi)),max(w/(2*np.pi))])
        ax2 = fig.add_subplot(gs[1, 0])

        ax2.plot(w/(2*np.pi), c0, 'g', linewidth=2, label='1st estimate')
        ax2.plot(w/(2*np.pi),cp,'r', linewidth=2, label='2nd estimate')
        ax2.plot(w/(2*np.pi),cp+2*np.diag(CM)[0:-1],'k--', linewidth=2, label='vel $\pm 2 \sigma$')
        ax2.plot(w/(2*np.pi),cp-2*np.diag(CM)[0:-1],'k--', linewidth=2)
        ax2.sharex(ax1)
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('$C_R(f)$')
        ax2.legend()

        ax3 = fig.add_subplot(gs[:, 1])
        #ax3 = fig.add_subplot(gs[-1, 0])
        ax3.plot(coord_uni[:,0],coord_uni[:,1],'g^')
        ax3.plot(p12[:,0],p12[:,1],'r^-')
        ax3.set_xlabel('X distance [m]')
        ax3.set_ylabel('Y distance [m]')
        folder2='AJUSTES_NUEVOS_GS_FINAL_RE/'+folder
        fig.suptitle(folder+'_'+str(nfold).zfill(2)+ ' ' +'Distance: '+str(np.round(r,2)) + ' [m] ' + med)
        plt.savefig(folder2+'/ajuste'+str(ii).zfill(3)+'.png', format='png', dpi=300, bbox_inches='tight',pad_inches=0.1)
        #plt.close()
        print('save',folder2+'/ajuste'+str(ii).zfill(3)+'.png')
        e1=(p12[0,0],p12[0,1])
        e2=(p12[1,0],p12[1,1])
        np.savetxt(folder2+'/cdd'+str(ii).zfill(3)+'.csv',np.vstack((w/(2*np.pi),cp,np.diag(CM)[0:-1])).T,
            fmt=['%3.3f','%3.3f','%3.3f'],header='Frequency,Velocity,Velstd \n coords,'+str(e1)+','+str(e2),delimiter=',')
        return rhop,cp,Ap,H2,CM,err,G,G1,G2

      


def get_dispersion_curve2_gs(f,CC,r,smooth_factor,f1,f2,chigh,clow,ii,coord_uni,p12,create_figures=False):
    w=2*np.pi*f
    wmin=w[0]
    wmax=w[-1]
    #cbounds=[500,1000]
    ## smooth CC to find reliable zero_crossings
    CC_smooth=smooth(CC,smooth_factor)
    ## Define window where calculate zero_crossings
    idx=np.where((f >=f1) & (f <=f2))[0]
    fx=f[idx]
    wx=fx*2*np.pi
    ## smooth CC to find reliable zero_crossings
    #CCX=smooth(CC,smooth_factor)
    CCX_smooth=CC_smooth[idx]
    CCX=CC[idx]
    wmin=wx[0]
    wmax=wx[-1]
    dw=wmax-wmin
    NN=4
    wtrial=np.logspace(np.log10(wmin+dw/6),np.log10(wmax-dw/6),NN)
    wtrial=np.logspace(np.log10(wmin),np.log10(wmax),NN)
    ftrial=wtrial/(2*np.pi)
    #chigh=400
    #clow=100
    NC=40
    cgrid=np.linspace(chigh,clow,NC)
    c_combs=np.array(list(combinations(cgrid,NN)))
    print(c_combs.shape)
    ccombs2=np.split(c_combs,10,axis=0)
    best_sol=1000000

    f1p=1
    f2p=35
    idxp=np.where((f >=f1p) & (f <=f2p))[0]
    fxp=f[idxp]
    CCXp=CC[idxp]
    CCXp_smooth=CC_smooth[idxp]
    
    for x in ccombs2:
            #print('a')
            #print(wtrial.shape,x.shape,wx.shape)
            print('b')
            ip=interpolate.interp1d(wtrial,x,'linear',fill_value='extrapolate')
            #ipw=np.interp(wx,wtrial,x)
            #
            print('c')

            ipw=ip(wx).astype(np.float32)
            #print('ipw',ipw.shape)

            #ipw2=ip2(wnew)
            #
            #print(np.floor(DJ * ((wx * r / ipw) - jargmin)).shape)
            print('d')

            #rho_pre=jv(0,wx*r/ipw)
            rho_pre=j0(wx*r/ipw)

            #rho_pre  = J0[np.floor(DJ * ((wx * r / ipw) - jargmin)).astype(int)]
            print('e')
            A_est = (rho_pre @ CCX_smooth)/(np.linalg.norm(rho_pre,axis=1))**2
            print('f')
            A_est_t=A_est.reshape(len(A_est),1)
            print('g')
            rho_pre=np.multiply(A_est_t, rho_pre)
            print('h')
            Drho=CCX_smooth-rho_pre
            print('i')
            err=np.linalg.norm(Drho,axis=1)**2

            agmin=np.argmin(err)
            Asol=A_est[agmin]
            rhosol=rho_pre[agmin]
            errsol=err[agmin]
            csol=ipw[agmin]
            ccombsol=x[agmin]
            print('best',best_sol)
            if err[agmin]<best_sol:
                best_sol=err[agmin]
                agmin=np.argmin(err)
                A_best=Asol
                rho_best=rhosol
                err_best=errsol
                c_best=csol
                ccombs_best=ccombsol

    lambdar=0.45
    cond=np.where( ccombs_best/ftrial/r >= lambdar)
    print('lambdar',ccombs_best/ftrial/r)
    print('wfor',wtrial)
    #wtrial=wtrial[cond]
    print('wafter',wtrial)
    popt,pcov=curve_fit(tanhfit,wtrial,ccombs_best)

    cbest_int=tanhfit(wx,*popt)
    print('parametros',popt)
    print('covarianzas',pcov)
    rho_best=jv(0, wx*r/(cbest_int))
    print('rho',rho_best.shape)
    A_best = (np.sum(rho_best *CCX_smooth))/(np.linalg.norm(rho_best))**2
    print('A',A_best.shape)
    rho_best=A_best*rho_best
   
   
    return fx,CCX,CCX_smooth,c_combs,ftrial,wtrial,A_best,rho_best,err_best,ccombs_best,cbest_int,fxp,CCXp,CCXp_smooth


def get_dispersion_curve2(f,CC,r,smooth_factor,f1,f2,ii,coord_uni,p12,create_figures=False):
    w=2*np.pi*f
    wmin=w[0]
    wmax=w[-1]
    #cbounds=[500,1000]
    ## smooth CC to find reliable zero_crossings
    CC_smooth=smooth(CC,smooth_factor)
    ## Define window where calculate zero_crossings
    idx=np.where((f >=f1) & (f <=f2))[0]
    fx=f[idx]
    wx=fx*2*np.pi
    ## smooth CC to find reliable zero_crossings
    #CCX=smooth(CC,smooth_factor)
    CCX_smooth=CC_smooth[idx]
    CCX=CC[idx]
    #CCX_smooth=smooth(CCX,smooth_factor)
    ## lims for plotting at the end
    f1p=1
    f2p=27
    idxp=np.where((f >=f1p) & (f <=f2p))[0]
    fxp=f[idxp]
    CCXp=CC[idxp]
    CCXp_smooth=CC_smooth[idxp]
    zero_crossings = np.where(np.diff(np.sign(CCX_smooth)))[0]

    #zero_crossings = np.where(np.diff(np.sign(CC_smooth)))[0]

    #wzero = f[zero_crossings]*2*np.pi
    print(np.diff(fx[zero_crossings]))
    print('frecuencias de cero crossing', fx[zero_crossings])
    if np.any( ( np.diff(fx[zero_crossings]) < 1.2 ) ):
        print('zero crossings are too close, try other smoothing factor')
        return
    wzero = fx[zero_crossings]*2*np.pi
    jn = jn_zeros(0, 400)
    cms=[]
    cms2=[]
    #plt.figure(1)
    m=np.arange(-2,3,1)
    for i in range(len(m)):
        #if m[i]>=0:
            cm = wzero * r / (jn[2 * i:2 * i + len(wzero)])
            #cm = wzero * r / (jn[i:i + len(wzero)])

            #print( i:i + len(wzero))
            print(2 * i,2 * i + len(wzero))
            #print(i,i + len(wzero))

            cms.append(cm)
            #plt.plot(wzero/(2*np.pi),cm,'o')
        #elif m[i]<0:

    
    print('minimos',(wzero[0]/(2*np.pi)),cms[0][0],cms[0][0]/(wzero[0]/(2*np.pi)))
    print('maximos',(wzero[-1]/(2*np.pi)),cms[0][-1],cms[0][-1]/(wzero[-1]/(2*np.pi)))

    print(wzero,cms[0])

    ##eliminar los zero crossings en donde la relacion lambda/r >= 0.45 no se cumpla
    lambdar=0.45
    cond=np.where(cms[0]/(wzero/(2*np.pi))/r >= lambdar)
    print('lambdar',cms[0]/(wzero/(2*np.pi))/r)
    wzero=wzero[cond]
    cmc=cms[0][cond]
    zero_crossings=zero_crossings[cond]
    print(cmc)
    ## condicion para los arrays

    try:
        popt,pcov=curve_fit(tanhfit,wzero,cmc)

        cbest_int=tanhfit(wx,*popt)
        print('parametros',popt)
        print('covarianzas',pcov)


    except TypeError:
        print('error: cant reliably build a dispersion curve from two zero crossings')

       
    
    rho_pre=jv(0, wx*r/(cbest_int))
    A_est = (np.sum(rho_pre *CCX_smooth))/(np.linalg.norm(rho_pre))**2
    rho_pre2=A_est*rho_pre
   
    return cmc,zero_crossings,wzero,r,jn,cbest_int,A_est,rho_pre2,fx,CCX,CCX_smooth,fxp,CCXp,CCXp_smooth



       
    

def wlstsq(m,G,H):
    return G.T @ G @ m + H.T @ H @ m

def newton_iteration(Ap,w,r,cp,c0,rho_obs,rhop,HH,N,f):

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
    return rhop,cp,Ap



def tanhfit(omega,d,e,f):
    c=d/np.tanh(e*omega)+f/np.sqrt(omega)
    return c

def tanh2fit(omega,d,e):
    c=d/np.tanh(e*omega)
    return c
