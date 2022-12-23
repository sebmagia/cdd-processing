import numpy as np
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

def newton_search2(w,rho_obs,rho_0,c0,A0,r,sigma_d,sigma_D,sigma_A,create_figures=False):
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
    H2 = np.zeros((N, N))
    for i in range(1, len(H2) - 1):
        H2[i, i - 1] = 1
        H2[i, i] = -2 
        H2[i, i + 1] = 1
    #H2=sigma_D**(-1)*H2
    H2=sigma_D*H2
    H[:N + 1, :N + 1] = H1
    H3 = H2[1:-1]
    H[N + 1:, :-1] = H3
    h = np.zeros((Ntop + Nbot))
    niter = 200
    H = scsp.csr_matrix(H)
    f = lambda x: wlstsq(x, G, H)
    HH = H.T @ h
    for iter in range(niter):
        Drho0=Drho

        G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
        G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
        #G = sigma_d**(-1)*scsp.hstack((G1, G2)).tocsr()
        G = sigma_d*scsp.hstack((G1, G2)).tocsr()

        ## delta rho, data vector
        #print(np.sum(Drho**2))
        F=G.T@Drho+HH

       
        L = LinearOperator((G.shape[1], G.shape[1]), matvec=f, rmatvec=f)
        Dm = cg(L, G.T @ Drho + HH, tol=1e-08, maxiter=4 * N)

        cpold=cp
        cp = cp + Dm[0][0:N]
        dcp=cp-cpold
        Ap = Ap + Dm[0][N]
        rhop = Ap * jv(0, w * r / cp)
        Drho = rho_obs - rhop
        #ratio=100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold))
        errpor=np.abs(100*(np.sum(Drho0)-np.sum(Drho))/(np.sum(Drho0)))
        print('errpor',errpor)
        #print( np.sum(Drho),np.sum(Drho0))
        #if errpor <1:
         #   break
        #print(( (np.sum((Drho)**2)-np.sum((Drho0)**2))/(np.sum((Drho0)**2)))*100)
        #print(ratio)
        #print('ratio',100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold)))

        #if ratio < 0.0001:
            
         #   break
    ## update G to generate posterior covariance matrix
    G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
    G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
    #G = (sigma_d)**-1 * scsp.hstack((G1, G2)).tocsr()
    G=(sigma_d) * scsp.hstack((G1, G2)).tocsr()
    CM=scsp.linalg.inv(scsp.vstack((G,H)).T@scsp.vstack((G,H))).todense()
    return rhop,cp,Ap,H2,CM

def get_dispersion_curve_newton(w,r,rho_0,rho_obs,f_full,rhofull_noise,rhofull_smooth,c0,A0,sigma_D,sigma_A,alpha,ii,coord_uni,p12,wzero,zcs):
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
    niter = 50
    H = scsp.csr_matrix(H)
    f = lambda x: wlstsq(x, G, H)
    HH = H.T @ h
    for iter in range(niter):
        Drho0=Drho

        G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
        G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
        #G = sigma_d**(-1)*scsp.hstack((G1, G2)).tocsr()
        G = sigma_D*scsp.hstack((G1, G2)).tocsr()

        ## delta rho, data vector
        #print(np.sum(Drho**2))
        F=G.T@Drho+HH

       
        L = LinearOperator((G.shape[1], G.shape[1]), matvec=f, rmatvec=f)
        Dm = cg(L, G.T @ Drho + HH, tol=1e-08, maxiter=4 * N)

        cpold=cp
        cp = cp + Dm[0][0:N]
        dcp=cp-cpold
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
    G=(sigma_D) * scsp.hstack((G1, G2)).tocsr()
    CM=scsp.linalg.inv(scsp.vstack((G,H)).T@scsp.vstack((G,H))).todense()
    err=np.abs(rho_obs-rhop)
    create_figures=True
    if create_figures:
        fig,ax=plt.subplots()
        ax.plot(w/(2*np.pi),cp,'r', linewidth=2, label='newton phase vel.')
        ax.plot(w/(2*np.pi),cp+2*np.diag(CM)[0:-1],'k--', linewidth=2, label='vel $\pm 2 \sigma$')
        ax.plot(w/(2*np.pi),cp-2*np.diag(CM)[0:-1],'k--', linewidth=2)
        ax.plot(w/(2*np.pi), c0, 'g', linewidth=2, label='grid serarch phase vel. inter')
        #plt.close()
        #ax.plot(w/(2*np.pi), cbest_grid, 'b', linewidth=2, label='grid serarch phase vel.')

        fig,ax=plt.subplots()
        #plt.plot(wnew/(2*np.pi),rho_obs,'r')
        plt.plot(w/(2*np.pi),rhop,'g')
        plt.close()
        #plt.plot(wnew/(2*np.pi),rho_grid,'k')

        fig,ax=plt.subplots()
        plt.semilogx((w/(2*np.pi)),cp)
        #ax.plot(wnew/(cp[::-1]*2*np.pi),cp[::-1],'r', linewidth=2, label='newton phase vel.')
        #ax.plot(wnew/(cp*2*np.pi),cp+2*np.diag(CM)[0:-1],'k--', linewidth=2, label='vel $\pm 2 \sigma$')
        #ax.plot(wnew/(cp*2*np.pi),cp-2*np.diag(CM)[0:-1],'k--', linewidth=2)
        #ax.plot(wnew/(cbest_int*2*np.pi), cbest_int, 'g', linewidth=2, label='grid serarch phase vel. inter')
        #ax.plot(wnew/(cbest_grid*2*np.pi), cbest_grid, 'b', linewidth=2, label='grid serarch phase vel.')

        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('$C_R(f)$')
        plt.close()

        #ax.set_ylim([85,350])
        ax.legend()
        fig = plt.figure(constrained_layout=True)
        gs = GridSpec(3, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, :])
        ax1.plot(w/(2*np.pi),cp,'r', linewidth=2, label='newton phase vel.')
        ax1.plot(w/(2*np.pi),cp+2*np.diag(CM)[0:-1],'k--', linewidth=2, label='vel $\pm 2 \sigma$')
        ax1.plot(w/(2*np.pi),cp-2*np.diag(CM)[0:-1],'k--', linewidth=2)
        ax1.plot(w/(2*np.pi), c0, 'g', linewidth=2, label='grid serarch phase vel. inter')

        #ax1.plot(wnew/(2*np.pi), cbest_grid, 'b', linewidth=2, label='grid serarch phase vel.')
        #ax1.set_ylim([85,450])
        # ax[0].plot(w, c_true, 'g--', linewidth=2, label='true phase vel.')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('$C_R(f)$')
        ax1.legend(loc=3,prop={'size':4})
        # # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
        ax2 = fig.add_subplot(gs[1, :])
        ax2.plot(f_full, rhofull_noise, 'k', linewidth=2, label='target corr func')
        ax2.plot(f_full, rhofull_smooth, 'r', linewidth=2, label='target corr func')
        #ax2.plot(w/(2*np.pi), rho_obs, 'k--', linewidth=2, label='target corr func')
        #ax2.plot(w/(2*np.pi), rho_grid, 'g-', linewidth=2, label='estimated corr func grid search')
        ax2.plot(w/(2*np.pi), rho_0, 'g', linewidth=2, label='estimated corr func interp grid')
        ax2.plot(w/(2*np.pi), rhop, 'g--', linewidth=2, label='estimated corr func newton')
        ax2.plot(wzero/(2*np.pi),rho_obs[zcs],'bo',label='Zero Crossings')
        #ax2.legend(loc=1,prop={'size':6})
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('Amplitude')
        ax2.legend()
        ax3 = fig.add_subplot(gs[-1, 0])
        ax3.plot(coord_uni[:,0],coord_uni[:,1],'g^')
        ax3.plot(p12[:,0],p12[:,1],'r^-')
        texto='Interstation Distance ' +str(r)+'[m]'
        ax3.text(12,15,texto,fontsize='8')
        ax3.set_xlabel('X  [m]')
        ax3.set_ylabel('Y  [m]')

        ax3.set_aspect('equal')
        ax4 = fig.add_subplot(gs[-1, 1])
        mat=ax4.imshow(CM,cmap='jet')
        fig.colorbar(mat,ax=ax4)
        plt.savefig('AJUSTES/ajuste'+str(ii).zfill(2)+'.png', format='png', dpi=300, bbox_inches='tight')
        plt.close()
        e1=(p12[0,0],p12[0,1])
        e2=(p12[1,0],p12[1,1])

        np.savetxt('AJUSTES/cdd'+str(ii).zfill(2)+'.csv',np.vstack((w/(2*np.pi),cp,np.diag(CM)[0:-1])).T,fmt=['%3.3f','%3.3f','%3.3f'],header='Frequency,Velocity,Velstd \n coords,'+str(e1)+','+str(e2),delimiter=',')
        plt.close()
        fig,ax=plt.subplots()
        ax.plot(f_full, rhofull_noise, 'k', linewidth=2, label='target noisy')
        ax.plot(f_full, rhofull_smooth, 'r', linewidth=2, label='target smoothed')
        #ax2.plot(w/(2*np.pi), rho_obs, 'k--', linewidth=2, label='target corr func')
        #ax2.plot(w/(2*np.pi), rho_grid, 'g-', linewidth=2, label='estimated corr func grid search')
        #ax.plot(w/(2*np.pi), rho_0, 'g', linewidth=2, label='estimated corr func interp grid')
        ax.plot(w/(2*np.pi), rhop, 'g', linewidth=2, label='estimated corr func newton')
        ax.plot(wzero/(2*np.pi),rho_obs[zcs],'bo',label='Zero Crossings')
        #ax2.legend(loc=1,prop={'size':6})
        ax.set_xlabel('Frequency (Hz)  ' +texto)
        ax.set_ylabel('Amplitude')
        plt.savefig('AJUSTES/ajuste_cc'+str(ii).zfill(2)+'.png', format='png', dpi=300, bbox_inches='tight')
        plt.close()
        ax.legend()

        #return rhop,cp,wnew,Ap,H2,CM,rho_grid

    return rhop,cp,Ap,H2,CM,err

def get_dispersion_curve22(f,CC,r,smooth_factor,f1,f2,ii,coord_uni,p12,create_figures=False):
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
    f2p=35
    idxp=np.where((f >=f1p) & (f <=f2p))[0]
    fxp=f[idxp]
    CCXp=CC[idxp]
    CCXp_smooth=CC_smooth[idxp]
    zero_crossings = np.where(np.diff(np.sign(CCX_smooth)))[0]

    wzero = fx[zero_crossings]*2*np.pi

    #rho_obs=CC
    cbounds=[60,500]
    #cbounds=[500,1000]

    clow=cbounds[0]
    chigh=cbounds[1]

    ## numero de nodos
    NN = len(zero_crossings)
    ## numero de frecuencias
    NC = 40
    N=len(wx)
    cgrid = np.linspace(chigh, clow, NC)
    ## all possible combs of phase velocity  nodes
    c_combs=np.array(list(combinations(cgrid,NN)))
    #print(c_combs.shape)

    ## number of trial models
    Nsearch=NC**NN
    ## for lstsq
    rho_obs=CCX_smooth
    wmin=wx[0]
    wmax=wx[-1]
    denom=rho_obs@rho_obs

    NJ=1001
    jargmin=0.5*wmin*r/chigh
    jargmax=2.0*wmax*r/clow
    LJ=jargmax-jargmin
    DJ=(NJ-1)/LJ
    jarg = np.linspace(jargmin, jargmax, NJ)
    nu = 0
    J0 = jv(nu, jarg)
    print(c_combs.nbytes)
    lenarr=N*(N-1)

    #ccombs2=np.split(ccombs,(int(N*(N-1))))
    print(c_combs.shape)

    #print(len(ccombs2))
    best_sol=1000000
    split_arr=True
    if split_arr:
        ccombs2=np.split(c_combs,10,axis=0)
        for x in ccombs2:
            print(x.shape)
            ip=interpolate.interp1d(wzero,x,'linear',fill_value='extrapolate')
            #
            ipw=ip(wx).astype(np.float32)
            #print('ipw',ipw.shape)

            #ipw2=ip2(wnew)
            #
            #print(np.floor(DJ * ((wx * r / ipw) - jargmin)).shape)
            rho_pre=jv(0,wx*r/ipw)
            #rho_pre  = J0[np.floor(DJ * ((wx * r / ipw) - jargmin)).astype(int)]
            A_est = (rho_pre @ rho_obs)/(np.linalg.norm(rho_pre,axis=1))**2
            A_est_t=A_est.reshape(len(A_est),1)
            rho_pre=np.multiply(A_est_t, rho_pre)
            Drho=rho_obs-rho_pre
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
                #print(ccombs_best)
                #print(c_best)
                #print(err_best)
                #print(A_best)
                #print(rho_best)
                #print(wzero)
        return A_best,rho_best,rho_obs,err_best,ccombs_best,c_best,wx


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
    f2p=35
    idxp=np.where((f >=f1p) & (f <=f2p))[0]
    fxp=f[idxp]
    CCXp=CC[idxp]
    CCXp_smooth=CC_smooth[idxp]
    zero_crossings = np.where(np.diff(np.sign(CCX_smooth)))[0]

    #zero_crossings = np.where(np.diff(np.sign(CC_smooth)))[0]

    #wzero = f[zero_crossings]*2*np.pi
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
    try:
        popt,pcov=curve_fit(tanhfit,wzero,cms[0])
        cbest_int=tanhfit(wx,*popt)

    except TypeError:
        print('just two points, use linear fit')
        y2=cms[0][-1]
        y1=cms[0][0]
        x2=wzero[-1]
        x1=wzero[0]
        m=(y2-y1)/(x2-x1)
        cbest_int=m*(wx-x1)+y1
        #cbest_int=(cms[0][-1]-cms[0][0])/(wzero[-1]-wzero[0])*(wx-cms[0][0])+cms[0][0]
    ## generate a new dispersion curve with a tanh shape.
    #plt.plot(fx,cbest_int)
    #cbest_int=tanhfit(wnew,*popt)
    #cbest_grid=c_best
    ## estimation of CC obtained from 1/tanh interpolation
    rho_pre=jv(0, wx*r/(cbest_int))
    A_est = (np.sum(rho_pre *CCX_smooth))/(np.linalg.norm(rho_pre))**2
    rho_pre2=A_est*rho_pre
    #plt.figure(5)
    #plt.plot(fx,rho_pre)
    #plt.figure(6)
    #plt.plot(f,CC)

    #(rho_pre @ CC)/(np.linalg.norm(rho_pre,axis=1))**2
    #rho_0 = A_best * jv(0, wnew * r / cbest_int)

    ### estimation of CC from grid search
    #rho_grid=rho_best

    # plt.figure(2)
    # plt.plot(fx,CCX,'m',label='raw')
    # plt.plot(fx,CCX_smooth,'c',label='smoothed')
    # plt.plot(fx[zero_crossings],CCX[zero_crossings],'ro',label='zero_raw')
    # plt.plot(fx[zero_crossings],CCX_smooth[zero_crossings],'ko',label='zero_smooth')
    # plt.plot(fx,rho_pre,'b',label='CC')
    # plt.plot(fx,rho_pre2,'k',label='CC')

    # plt.figure(3)
    # lembda=cms[0]/(fx[zero_crossings])
    # plt.plot(lembda,cms[0],'go')
    # print(lembda,cms[0])
    # print('lembda/r',lembda/r)
    #np.savetxt('AJUSTES/cdd'+str(ii).zfill(2)+'.dat',np.vstack((wnew/(2*np.pi),cp,np.diag(CM)[0:-1])).T,fmt=['%3.3f','%3.3f','%3.3f'])

    return cms,zero_crossings,wzero,r,jn,cbest_int,A_est,rho_pre2,fx,CCX,CCX_smooth,fxp,CCXp,CCXp_smooth


def get_dispersion_curve(f,CC,r,sigmad,alpha,sigma,ii,coord_uni,p12,create_figures=False):
    
    w=2*np.pi*f
    print(len(w))
    wmin=w[0]
    wmax=w[-1]
    #wnew=np.linspace(wmin,wmax,len(w)//2)
    wnew=np.linspace(wmin,wmax,len(w)//4)

    interp=interpolate.interp1d(w,CC,kind='cubic')
    rho_obs=interp(wnew)
    #wnew=w
    #rho_obs=CC
    cbounds=[85,450]
    #cbounds=[500,1000]

    clow=cbounds[0]
    chigh=cbounds[1]

    ## numero de nodos
    NN = 5
    ## numero de frecuencias
    NC = 35
    N=len(wnew)
    cgrid = np.linspace(chigh, clow, NC)
    ## all possible combs of phase velocity  nodes
    c_combs=np.array(list(combinations(cgrid,NN)))
    print(c_combs.shape)
    #c_combs=np.array(list(permutations(cgrid,NN)))

    #c_combs
    print(len(c_combs))
    print(c_combs.nbytes)
    #return c_combs,cgrid
    print(wnew,'xd')
    print(c_combs.nbytes)
    A_best, rho_best, err_best, cnodes_best, c_best, wgrid = grid_search_cw4(N, wnew, rho_obs, r, NN, NC, cbounds,
                                                                                c_combs, False)
    
    #return A_best,rho_best,err_best,cnodes_best,c_best,wgrid,wnew
    ## fit 1/tanh function to the best grid search solution
    #popt,pcov=curve_fit(tanhfit,wnew,c_best)
    ## fit tanh curve from nodes
    popt,pcov=curve_fit(tanhfit,wgrid,cnodes_best)
    ## generate a new dispersion curve with a tanh shape.
    cbest_int=tanhfit(wnew,*popt)

    #cbest_int=tanhfit(wnew,*popt)
    cbest_grid=c_best
    ## estimation of CC obtained from 1/tanh interpolation
    rho_0 = A_best * jv(0, wnew * r / cbest_int)

    ## estimation of CC from grid search
    rho_grid=rho_best

    print('weanding')
    rhop, cp, Ap,H2,CM= newton_search(wnew,rho_obs,rho_0,rho_grid,cbest_int,cbest_grid,A_best,r,sigmad,sigma,alpha)
    print('initial inter',np.sum(np.abs((rho_0-rho_obs))))
    print('initial grid ',np.sum(np.abs((rho_best-rho_obs))))
    print('newton',np.sum((np.abs(rhop-rho_obs))))

    if create_figures:
        fig,ax=plt.subplots()
        ax.plot(wnew/(2*np.pi),cp,'r', linewidth=2, label='newton phase vel.')
        ax.plot(wnew/(2*np.pi),cp+2*np.diag(CM)[0:-1],'k--', linewidth=2, label='vel $\pm 2 \sigma$')
        ax.plot(wnew/(2*np.pi),cp-2*np.diag(CM)[0:-1],'k--', linewidth=2)
        ax.plot(wnew/(2*np.pi), cbest_int, 'g', linewidth=2, label='grid serarch phase vel. inter')
        ax.plot(wnew/(2*np.pi), cbest_grid, 'b', linewidth=2, label='grid serarch phase vel.')

        fig,ax=plt.subplots()
        #plt.plot(wnew/(2*np.pi),rho_obs,'r')
        plt.plot(wnew/(2*np.pi),rhop,'g')
        #plt.plot(wnew/(2*np.pi),rho_grid,'k')

        fig,ax=plt.subplots()
        plt.semilogx((wnew/(2*np.pi)),cp)
        #ax.plot(wnew/(cp[::-1]*2*np.pi),cp[::-1],'r', linewidth=2, label='newton phase vel.')
        #ax.plot(wnew/(cp*2*np.pi),cp+2*np.diag(CM)[0:-1],'k--', linewidth=2, label='vel $\pm 2 \sigma$')
        #ax.plot(wnew/(cp*2*np.pi),cp-2*np.diag(CM)[0:-1],'k--', linewidth=2)
        #ax.plot(wnew/(cbest_int*2*np.pi), cbest_int, 'g', linewidth=2, label='grid serarch phase vel. inter')
        #ax.plot(wnew/(cbest_grid*2*np.pi), cbest_grid, 'b', linewidth=2, label='grid serarch phase vel.')

        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('$C_R(f)$')
        ax.set_ylim([85,350])
        ax.legend()
        fig = plt.figure(constrained_layout=True)
        gs = GridSpec(3, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, :])
        ax1.plot(wnew/(2*np.pi),cp,'r', linewidth=2, label='newton phase vel.')
        ax1.plot(wnew/(2*np.pi),cp+2*np.diag(CM)[0:-1],'k--', linewidth=2, label='vel $\pm 2 \sigma$')
        ax1.plot(wnew/(2*np.pi),cp-2*np.diag(CM)[0:-1],'k--', linewidth=2)
        ax1.plot(wnew/(2*np.pi), cbest_int, 'g', linewidth=2, label='grid serarch phase vel. inter')

        ax1.plot(wnew/(2*np.pi), cbest_grid, 'b', linewidth=2, label='grid serarch phase vel.')
        ax1.set_ylim([85,450])
        # ax[0].plot(w, c_true, 'g--', linewidth=2, label='true phase vel.')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('$C_R(f)$')
        ax1.legend(loc=3,prop={'size':4})
        # # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
        ax2 = fig.add_subplot(gs[1, :])
        ax2.plot(wnew/(2*np.pi), rho_obs, 'k--', linewidth=2, label='target corr func')
        ax2.plot(wnew/(2*np.pi), rho_grid, 'g-', linewidth=2, label='estimated corr func grid search')
        ax2.plot(wnew/(2*np.pi), rho_0, 'm', linewidth=2, label='estimated corr func interp grid')
        ax2.plot(wnew/(2*np.pi), rhop, 'r:', linewidth=2, label='estimated corr func newton')
        #ax2.legend(loc=1,prop={'size':6})
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('Amplitude')
        ax2.legend()
        ax3 = fig.add_subplot(gs[-1, 0])
        ax3.plot(coord_uni[:,0],coord_uni[:,1],'g^')
        ax3.plot(p12[:,0],p12[:,1],'r^-')
        texto='Interstation Distance ' +str(r)+'[m]'
        ax3.text(12,15,texto,fontsize='8')
        ax3.set_xlabel('X  [m]')
        ax3.set_ylabel('Y  [m]')

        ax3.set_aspect('equal')
        ax4 = fig.add_subplot(gs[-1, 1])
        mat=ax4.imshow(CM,cmap='jet')
        fig.colorbar(mat,ax=ax4)
        plt.savefig('AJUSTES/ajuste'+str(ii).zfill(2)+'.png', format='png', dpi=300, bbox_inches='tight')
        #plt.close()
        np.savetxt('AJUSTES/cdd'+str(ii).zfill(2)+'.dat',np.vstack((wnew/(2*np.pi),cp,np.diag(CM)[0:-1])).T,fmt=['%3.3f','%3.3f','%3.3f'])
        return rhop,cp,wnew,Ap,H2,CM,rho_grid

       
    

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

def newton_search(w,rho_obs,rho_0,rho_grid,c0,cgrid,A0,r,sigma_d,sigma_D,sigma_A,create_figures=False):
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
    H2 = np.zeros((N, N))
    for i in range(1, len(H2) - 1):
        H2[i, i - 1] = 1
        H2[i, i] = -2 
        H2[i, i + 1] = 1
    #H2=sigma_D**(-1)*H2
    H2=sigma_D*H2
    H[:N + 1, :N + 1] = H1
    H3 = H2[1:-1]
    H[N + 1:, :-1] = H3
    h = np.zeros((Ntop + Nbot))
    niter = 200
    H = scsp.csr_matrix(H)
    f = lambda x: wlstsq(x, G, H)
    HH = H.T @ h
    for iter in range(niter):
        Drho0=Drho

        G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
        G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
        #G = sigma_d**(-1)*scsp.hstack((G1, G2)).tocsr()
        G = sigma_d*scsp.hstack((G1, G2)).tocsr()

        ## delta rho, data vector
        #print(np.sum(Drho**2))
        F=G.T@Drho+HH

       
        L = LinearOperator((G.shape[1], G.shape[1]), matvec=f, rmatvec=f)
        Dm = cg(L, G.T @ Drho + HH, tol=1e-08, maxiter=4 * N)

        cpold=cp
        cp = cp + Dm[0][0:N]
        dcp=cp-cpold
        Ap = Ap + Dm[0][N]
        rhop = Ap * jv(0, w * r / cp)
        Drho = rho_obs - rhop
        #ratio=100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold))
        errpor=np.abs(100*(np.sum(Drho0)-np.sum(Drho))/(np.sum(Drho0)))
        print('errpor',errpor)
        #print( np.sum(Drho),np.sum(Drho0))
        if errpor <1:
            break
        #print(( (np.sum((Drho)**2)-np.sum((Drho0)**2))/(np.sum((Drho0)**2)))*100)
        #print(ratio)
        #print('ratio',100*np.sum(np.abs(dcp))/np.sum(np.abs(cpold)))

        #if ratio < 0.0001:
            
         #   break
    ## update G to generate posterior covariance matrix
    G1 = scsp.lil_matrix(np.diag(Ap * jv(1, w * r / cp) * w * r * cp ** (-2)))
    G2 = scsp.lil_matrix(jv(0, w * r / c0)).transpose()
    #G = (sigma_d)**-1 * scsp.hstack((G1, G2)).tocsr()
    G=(sigma_d) * scsp.hstack((G1, G2)).tocsr()
    CM=scsp.linalg.inv(scsp.vstack((G,H)).T@scsp.vstack((G,H))).todense()
    return rhop,cp,Ap,H2,CM

       
    if create_figures:
        fig = plt.figure(constrained_layout=True)
        gs = GridSpec(3, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, :])
        ax1.plot(w/(2*np.pi), c0, 'k', linewidth=2, label='grid interp phase vel.')
        #ax[0].plot(w/(2*np.pi), cp, 'r', linewidth=2, label='newton phase vel.')
        ax1.plot(w/(2*np.pi), cgrid, 'b', linewidth=2, label='grid serarch phase vel.')
        ax1.errorbar(w/(2*np.pi),cp,np.diag(CM)[0:-1],'r', linewidth=2, label='newton phase vel.')

        ax1.set_xlabel('w')
        ax1.set_ylabel('c(w)')
        ax1.legend()        #ax1.plot(fx,Cx)
        # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
        ax2 = fig.add_subplot(gs[1, :])
        ax2.plot(w/(2*np.pi), rho_obs, 'k--', linewidth=2, label='target corr func')
        ax2.plot(w/(2*np.pi), rho_grid, 'g-', linewidth=2, label='estimated corr func grid search')
        ax2.plot(w/(2*np.pi), rho_0, 'm', linewidth=2, label='estimated corr func interp grid')
        ax2.plot(w/(2*np.pi), rhop, 'r:', linewidth=2, label='estimated corr func newton')
        ax2.legend()
        ax2.set_xlabel('w')
        #ax3 = fig.add_subplot(gs[1:, -1])
        ax3 = fig.add_subplot(gs[-1, 0])
        ax3.remove()
        ax3.imshow(CM,cmap='jet')
        ax3.set_colorbar()
        ax4 = fig.add_subplot(gs[-1, 1])

        plt.savefig('ajuste10.png', format='png', dpi=300, bbox_inches='tight')

        
def tanhfit(omega,d,e,f):
    c=d/np.tanh(e*omega)+f/np.sqrt(omega)
    return c

def tanh2fit(omega,d,e):
    c=d/np.tanh(e*omega)
    return c
def grid_search_cw4(N,wnew,rho_obs,r,NN,NC,cbounds,ccombs,create_figures,split_arr=False):
    #wmin=2*np.pi/Tmax
    #wmax=2*np.pi/Tmin
    iw = np.arange(0, N, 1)
    #w = np.linspace(wmin, wmax, N)
    w=wnew
    wmin=w[0]
    wmax=w[-1]
    #return w
    #T = 1 / w
    ## w-index of nodes, where we search for phase velocities
    igrid = np.floor((N - 1) * np.arange(0, NN) / (NN - 1)).astype(int)
    print(igrid)
    igrid2 = np.floor((N - 1) * np.arange(0, 2*NN) / (2*NN - 1)).astype(int)
    print(igrid2)
    #ii=np.linspace(0,len(w),10).astype(int)
    wgrid=w[igrid]
    #print(ii)
    #wgrid=np.zeros(5)

    #wgrid[0]=w[igrid2[0]]
    #wgrid[1]=w[igrid2[2]]
    #wgrid[2]=w[igrid2[3]]
    #wgrid[3]=w[igrid2[6]]
    #wgrid[4]=w[igrid2[9]]
    #print(igrid2[0],igrid2[2],igrid2[3],igrid2[6],igrid2[9])
    # print(wgrid,igrid,ii)
    #wgrid=w[igrid]
    #print(w[igrid])

    #igrid=w[]
    #wgrid=w[]
    #wgrid=w[igrid]
    #print(len(w),len(w)//6)
    #print(igrid,wgrid)
    #return wgrid,igrid
    clow=cbounds[0]
    chigh=cbounds[1]
    print(NC)
    cgrid=np.linspace(clow,chigh,NC)
    ## number of trial models
    Nsearch=NC**NN
    ## for lstsq
    denom=rho_obs@rho_obs

    NJ=1001
    jargmin=0.5*wmin*r/chigh
    jargmax=2.0*wmax*r/clow
    LJ=jargmax-jargmin
    DJ=(NJ-1)/LJ
    jarg = np.linspace(jargmin, jargmax, NJ)
    nu = 0
    J0 = jv(nu, jarg)
    print(ccombs.nbytes)
    lenarr=N*(N-1)

    #ccombs2=np.split(ccombs,(int(N*(N-1))))
    print(ccombs.shape)

    #print(len(ccombs2))
    best_sol=1000000
    split_arr=False
    if split_arr:
        ccombs2=np.split(ccombs,220,axis=0)
        for x in ccombs2:
            print(x.shape)
            ip=interpolate.interp1d(wgrid,x,'linear',fill_value='extrapolate')
            #
            ipw=ip(w).astype(np.float32)
            print('ipw',ipw.shape)

            #ipw2=ip2(wnew)
            #
            rho_pre  = J0[np.floor(DJ * ((w * r / ipw) - jargmin)).astype(int)]
            A_est = (rho_pre @ rho_obs)/(np.linalg.norm(rho_pre,axis=1))**2
            A_est_t=A_est.reshape(len(A_est),1)
            rho_pre=np.multiply(A_est_t, rho_pre)
            Drho=rho_obs-rho_pre
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
                print(ccombs_best)
                print(c_best)
                print(err_best)
                print(A_best)
                print(rho_best)
                print(wgrid)
        return A_best,rho_best,err_best,ccombs_best,c_best,wgrid

    else:
        ip=interpolate.interp1d(wgrid,ccombs,'linear',fill_value='extrapolate')
            #
        ipw=ip(w).astype(np.float32)
        print('ipw',ipw.shape)

        #ipw2=ip2(wnew)
        #
        rho_pre  = J0[np.floor(DJ * ((w * r / ipw) - jargmin)).astype(int)]
        A_est = (rho_pre @ rho_obs)/(np.linalg.norm(rho_pre,axis=1))**2
        A_est_t=A_est.reshape(len(A_est),1)
        rho_pre=np.multiply(A_est_t, rho_pre)
        Drho=rho_obs-rho_pre
        err=np.linalg.norm(Drho,axis=1)**2

        agmin=np.argmin(err)
        A_best=A_est[agmin]
        rho_best=rho_pre[agmin]
        err_best=err[agmin]
        c_best=ipw[agmin]
        ccombs_best=ccombs[agmin]
        print(ccombs_best)
        print(c_best)
        print(err_best)
        print(A_best)
        print(rho_best)
        print(wgrid)

        if create_figures:
            fig,ax=plt.subplots(2,1)
            ax[0].plot(w, c_best, 'r', linewidth=2, label='est phase vel')
            ax[0].plot(wgrid, ccombs_best, 'go', linewidth=2)
            ax[0].set_xlabel('w')
            ax[0].set_ylabel('c(w)')
            ax[0].legend()

            ax[1].plot(w, rho_obs, 'k--', linewidth=2, label='target corr func')
            ax[1].plot(w, rho_best, 'r:', linewidth=2, label='estimated corr func')
            ax[1].legend()
            ax[1].set_xlabel('w')
        return A_best,rho_best,err_best,ccombs_best,c_best,wgrid


    

def grid_search_cw3(N,wmin,wmax,rho_obs,r,NN,NC,cbounds,ccombs,create_figures):
    #wmin=2*np.pi/Tmax
    #wmax=2*np.pi/Tmin
    iw = np.arange(0, N, 1)
    w = np.linspace(wmin, wmax, N)
    #return w
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
   
    #return wgrid,J0,rho_obs,jargmin,jargmax
    for n,x in enumerate(ccombs):
        #print(x)
        #print(n)
        #print(n,len(ccombs))
        #print('as',len(wgrid),len(x))
        #print('xd',len(wgrid),len(x))
        ip = interpolate.interp1d(wgrid, x, 'linear', fill_value='extrapolate')
        ipw=ip(w)
        #print('ipw',ip(w))
        #print('diff',ipw[0]-ipw[1])
        #return ipw
        #print('diff',np.sign((ip(w[0])-ip(w[-1]))))
        #rho_pre=jv(0,w*r/ip(w))
        rho_pre = J0[np.floor(Dj * ((w * r / ip(w)) - jargmin)).astype(int)]
        A_est = (rho_pre @ rho_obs) / (rho_pre@rho_pre)
        rho_pre = A_est * rho_pre
        Drho = rho_obs - rho_pre
        err = np.sum(Drho * Drho)
        
        #print(A_est,err)

        #idx = np.argmin(err)
        if n == 0:
            A_best = A_est
            rho_best = rho_pre
            err_best = err
            ccombs_best = ccombs[n]
            c_best = ip(w)
        elif n > 0 and np.min(err) < err_best:
            # print('foo',np.min(err),cnodes_best)
            A_best = A_est
            rho_best = rho_pre
            err_best = err
            ccombs_best = ccombs[n]
            c_best = ip(w)
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

    return A_best,rho_best,err_best,ccombs_best,c_best,wgrid
