import TimeSeries as ts
import os
import multiprocessing
from itertools import product
from functools import partial
from scipy.special import jv
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,LinearOperator,aslinearoperator
import time

import itertools
from scipy.signal import detrend,hann,coherence,correlate,periodogram,correlation_lags,tukey,butter,sosfilt,spectrogram,welch,csd
from scipy import interpolate,integrate
from scipy.fft import fft,ifft,fftshift,ifftshift,fftfreq,rfftfreq,rfft,fft2,irfft
if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plt

    print('cordero')
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_DENSEARRAY/M4/'

    path='/home/doctor/Doctor/Magister/Tesis/databases/process_data/Pent_7_30'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/SUR2B/M7R/'
    carpetas=os.listdir(path)
    path2 = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/SUR2B/Outputs/M7R'

    #path='/home/doctor/Doctor/Magister/Tesis/databases/SUR1B/M4/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/SUR2A/M4/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/SUR2B/M3/'


    #name='M4.png'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_SN_P2/M1/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_SN_S1_PILZ/M5/'
    #path='/home/doctor/Doctor/Magister/Tesis/databases/LP_SN_S2_PILZ/M8B/'
    if 'SUR1A' in path:
        name = 'SUR1A/'+path.split('/')[-2] + '.png'
        name2='SUR1A/'+path.split('/')[-2] + '_TS.png'
    if 'SUR1B' in path:
        name = 'SUR1B/' + path.split('/')[-1] + '.png'
    if 'SUR2A' in path:
        name = 'SUR2A/' + path.split('/')[-1] + '.png'
    if 'SUR2B' in path:
        name = 'SUR2B/' + path.split('/')[-1] + '.png'

    path_ele=sorted(os.listdir(path))
    ## check if elements are directories or files
    dire=[os.path.isdir(path+x) for x in path_ele]
    dire=['/Pent_7_30_001part1', '/Pent_7_30_002part1', '/Pent_7_30_003part1', '/Pent_7_30_004part1', '/Pent_7_30_005part1']
    paths=[path+x+'/EqualizedFile.dat' for x in dire]
    #paths=[path+x+'/EqualizedFile.dat' for x,y in zip(path_ele,dire) if y]
    TS=ts.TimeSeries(paths,mp=False)
    # del TS
    ignore_this=False
    if ignore_this:
        TS2=ts.TimeSeries(path2,mp=False,create_data=False)
        TS2.get_NCF(TS2.C, TS2.w, name)
        TS.plot_traces(name2)
    #r1=TS2.datas_eff[0]
    #r2=TS2.datas_eff[1]
    #outliers=np.where(np.abs(r2)>np.mean(r2)+5*np.std(r2))[0]
    #r3=TS2.datas_eff[2]
    #r4=TS2.datas_eff[3]
    #r5=TS2.datas_eff[4]
    # from scipy.signal import medfilt
    # r22=medfilt(r2)
    # xd1=correlate(r1,r2)
    # xd2=correlate(r1,r3)
    # xd3=correlate(r1,r4)
    # xd4=correlate(r1,r5)
    # xd5=correlate(r2,r3)
    # xd6=correlate(r2,r4)
    # xd7=correlate(r2,r5)
    # xd8=correlate(r3,r4)
    # xd9=correlate(r3,r5)
    # xd10=correlate(r4,r5)
    # idx=correlation_lags(r1.size,r2.size)
    # t=np.arange(-r2.size/512,r2.size/512,1./512)[1:]
    # plt.plot(TS.t_eff,r2)
    # plt.plot(TS.t_eff[outliers],r2[outliers])
    #
    # tl1=t[np.argmax(xd1)]
    # tl2=t[np.argmax(xd2)]
    # tl3=t[np.argmax(xd3)]
    # tl4=t[np.argmax(xd4)]
    # tl5=t[np.argmax(xd5)]
    # tl6=t[np.argmax(xd6)]
    # tl7=t[np.argmax(xd7)]
    # tl8=t[np.argmax(xd8)]
    # tl9=t[np.argmax(xd9)]
    # tl10=t[np.argmax(xd10)]


    #import TimeSeries as ts
    # data1=TS.datas_eff[2]
    # data2=TS.datas_eff[4]
    # rij=correlate(data1,data2)
    # t = np.arange(-TS.datas_eff[0].size / 512, TS.datas_eff[1].size / 512, 1. / 512)[1:]
    # tmax = t[np.argmax(rij)]
    #
    # data1b=data1[512 * 3: -1]
    # data2b=data2[0:-1-512*3]
    # t2 = np.arange(-data1b.size / 512, data2b.size / 512, 1. / 512)[1:]
    #
    # rij2=correlate(data1b,data2b)
    # t2max = t2[np.argmax(rij2)]
    #
    #
    # data3=TS.datas_eff[0]
    # data4=TS.datas_eff[2]
    # rij = correlate(data3, data4)
    # t = np.arange(-TS.datas_eff[0].size / 512, TS.datas_eff[1].size / 512, 1. / 512)[1:]
    # tmax = t[np.argmax(rij)]
    # data3b = data3[0:-1-512*3]
    # data4b = data4[512 * 3: -1]
    #
    # t3 = np.arange(-data3b.size / 512, data4b.size / 512, 1. / 512)[1:]
    #
    # rij3 = correlate(data3b, data4b)
    # t3max = t3[np.argmax(rij3)]
    # #tmax2 = t[np.argmax(rij2)]
    # print('wow', t.shape, rij.shape)
    # a,b=TS.correlations(figure=False)
    # plt.plot(TS.distances,np.abs(a),'go')
    # plt.plot(TS.tv[:100000],TS.datas_eff[2][:100000],'r')
    # plt.plot(TS.tv[:100000],TS.datas_eff[4][:100000],'g')
    # corr=correlate(TS.datas_eff[0],TS.datas_eff[1])
    # lags=correlation_lags(TS.datas_eff[0].size,TS.datas_eff[1].size)
    # r1=TS.datas_eff[0]
    # r2=TS.datas_eff[1]
    # r3=TS.datas_eff[2]
    # r4=TS.datas_eff[3]
    # r5=TS.datas_eff[4]
    #
    # xd1=correlate(r1,r2)
    # xd2=correlate(r1,r3)
    # xd3=correlate(r1,r4)
    # xd4=correlate(r1,r5)
    # xd5=correlate(r2,r3)
    # xd6=correlate(r2,r4)
    # xd7=correlate(r2,r5)
    # xd8=correlate(r3,r4)
    # xd9=correlate(r3,r5)
    # xd10=correlate(r4,r5)
    # idx=correlation_lags(r1.size,r2.size)
    # t=np.arange(-r1.size/512,r1.size/512,1./512)[1:]
    # tl1=t[np.argmax(xd1)]
    # tl2=t[np.argmax(xd2)]
    # tl3=t[np.argmax(xd3)]
    # tl4=t[np.argmax(xd4)]
    # tl5=t[np.argmax(xd5)]
    # tl6=t[np.argmax(xd6)]
    # tl7=t[np.argmax(xd7)]
    # tl8=t[np.argmax(xd8)]
    # tl9=t[np.argmax(xd9)]
    # tl10=t[np.argmax(xd10)]
    # plt.plot(t,xd1)
    # plt.plot(t,xd2)

    #TS.get_NCF(TS.C,TS.w,name)
    pata = False
    if pata:
        idx = np.where((TS.w < 29) & (TS.w > (6)))[0]
        wp = 2 * np.pi * TS.w[idx]
        wmin = wp[0]
        wmax = wp[-1]
        data = TS.C[4][idx]
        N = len(wp)
        tic = time.time()
        wf = 2 * np.pi * TS.w[idx]
        r = 33.61
        cbounds = [10, 30]
        tic = time.time()
        clow = cbounds[0]
        chigh = cbounds[1]
        ## numero de nodos
        NN=4
        ## numero de frecuencias
        NC=40
        cgrid = np.linspace(clow, chigh, NC)
        #ns = np.arange(0, 30)
        f=partial(product,repeat=NN)
        g=f(cgrid)
        cnodes=np.split(np.array(list(g)),256)

        tic = time.time()
        A_best, rho_best, err_best, cnodes_best, c_best= ts.grid_search_cw2(N, wmin, wmax, data, r, NN, NC, cbounds,cnodes, True)
        toc = time.time()
        print('gridsearch linear Done in {:.4f} seconds'.format(toc - tic))
        rho0 = A_best[0] * jv(0, wf * r / c_best)

        tic = time.time()
        rho,cp,ap=ts.newton_search(wf, data, rho0, c_best, A_best, r)
        wtod,cc,zc=ts.zero_crossings(wp,data,r)
        toc = time.time()
        print('newtonsearch linear Done in {:.4f} seconds'.format(toc - tic))


    # cc=pool.map(f,cgrid)
    # pool.close()
    # pool.join()
    #g=partial(product,repeat=3)
    # from itertools import product
    #
    # #perms_o = product(ns, repeat=3)
    # import time
    # tic = time.time()
    # pool = multiprocessing.Pool()
    # perms2=np.array(pool.apply(product(cgrid),args=(3)))
    # toc = time.time()
    # print('aDone in {:.4f} seconds'.format(toc - tic))


    # perms = list(product(ns, repeat=3))
    # nnodes=3
    #
    # NP=len(ns)**nnodes
    # #data=np.zeros((len(perms),3))
    #
    #
    # dd=ts.get_nodes_new(cgrid,NP,3,perms_o)
    # pool = multiprocessing.Pool()
    # tic = time.time()
    # R=pool.map(partial(ts.get_nodes_new,cgrid,NP,3),perms)
    # pool.close()
    # pool.join()
    # #out = pool.map(partial(preprocess_data, self.datas_eff, self.fs, self.s, False), self.combs)


    # def myfunc(iter):
    #
    #     return sum(list(iter))
    #
    # def func(perms_o):
    #     for x in perms_o:
    #         print(x)
    # #cn=get_nodes_new(cgrid,product(ns,repeat=3))
    #
    # #A_best, rho_best, err_best, cnodes_best, c_best = ts.grid_search_cw(N, wmin, wmax, data, r, 4, 40, cbounds, True)
    # from itertools import product
    #

    # func(perms_o)
    # perms=np.array(list(product(ns,repeat=3)))
    # arr=np.zeros(perms.shape)
    # for x in product(ns,repeat=3):
    #     print(x)
    #
    # perms = [np.array(p) for p in product(ns, repeat=3)]
    # ## all possible phase velocities
    # c_nodes = np.array([[cgrid[x[y]] for y in range(3)] for x in perms])
    # toc = time.time()
    # print('aDone in {:.4f} seconds'.format(toc - tic))

    # pool = multiprocessing.Pool()
    # tic = time.time()
    # A_best, rho_best, err_best, cnodes_best, c_best = pool.apply(ts.grid_search_cw,args=(N, wmin, wmax, data, r, 4, 40, cbounds, False))
    # toc = time.time()
    # print('aDone in {:.4f} seconds'.format(toc - tic))
    #
    # pool.close()
    # pool.join()
    # print('done')
    # pool = multiprocessing.Pool()
    # tic = time.time()
    # R=pool.map(ts.manage_data_new2,paths)
    # datas=[R[x][0] for x in range(len(R))]
    # metadatas=[R[x][1] for x in range(len(R))]
    # pool.close()
    # pool.join()
    # toc = time.time()
    #print('aDone in {:.4f} seconds'.format(toc - tic))
    cont = False
    if cont:
        #TS=ts.TimeSeries(paths)

        #pool = multiprocessing.Pool()
        #tic = time.time()
        #TS=pool.apply(ts.TimeSeries,args=(paths[0]))
        #TS=ts.TimeSeries(paths)
        mins=len(TS.tt)/512/60
        print(mins)
        #TS.plot_traces()
        #coords=[(x.get('Latitude'),x.get('Longitude')) for x in TS.metadatas]

        Cxy,w,combs=TS.preprocess_data_pilz2(s=8,filter=True)

        TS.get_NCF(Cxy,w)
        wp=np.split(w,2)[0]
        idx=np.where((w<29) & (wp > (6)) )[0]
        plt.figure(1)
        plt.plot(wp[idx],Cxy[4][idx])

        data=Cxy[4][idx]
        wf=2*np.pi*wp[idx]
        r=33.61
        wmin=wf[0]
        ts.zero_crossings(wf,data,r)

        wmax=wf[-1]
        #cbounds=[50,20000] bad fot high
        #cbounds=50-700 for good fit
        #cbounds=10-50 for bad fit low
        cbounds=[10,30]
        N=len(wf)
        pool = multiprocessing.Pool()
        tic = time.time()
        A_best, rho_best, err_best, cnodes_best, c_best=pool.apply(ts.grid_search_cw,args=(N, wmin, wmax, data, r, 4, 40, cbounds,True))
        pool.close()
        pool.join()
        toc = time.time()
        print('aDone in {:.4f} seconds'.format(toc - tic))

        tic = time.time()

        A_best, rho_best, err_best, cnodes_best, c_best = ts.grid_search_cw(N, wmin, wmax, data, r, 4, 40, cbounds,True)
        #                                                                     create_figures=True)
        toc = time.time()
        print('bDone in {:.4f} seconds'.format(toc - tic))

        rho0 = A_best[0] * jv(0, wf * r / c_best)
        pool = multiprocessing.Pool()
        results=pool.apply(ts.newton_search,args=(wf,data,rho0,c_best,A_best,r))
        pool.close()
        pool.join()

        print('alberto')
    #rhop,cp,Ap=ts.newton_search(wf,data,rho0,c_best,A_best,r)

    # c0 = c_best
    # A0 = A_best[0]
    # rho0 = A0 * jv(0, wf * r / c0)
    # Drho0 = data - rho0
    # E = Drho0 @ Drho0
    #
    # # starting values of linearized inversion
    # cp = c0;
    # Ap = A0;
    # rhop = rho0;
    #
    #
    # # matrix H of prior constraints
    # # first half, smallness
    # # second half, smoothness
    # def wlstsq(m, G, H):
    #     return G.T @ G @ m + H.T @ H @ m
    #
    #
    # # SECONDDERIV = 1;  % second derivative smoothing on 1, first deviv on 0
    # ## Ntop is applied over length(c)+length(A)=N+1 parameters
    # ## Nbot is Laplacian Smoothing applied over N-1 c(w) parameters
    # Ntop = N + 1;
    # Nbot = N - 2
    # small = 0.01
    # H1 = 0.01 * np.eye(Ntop)
    # smooth = 50
    # H2 = np.zeros((N, N))
    # H = np.zeros((Ntop + Nbot, N + 1))
    # for i in range(1, len(H2) - 1):
    #     H2[i, i - 1] = smooth
    #     H2[i, i] = -2 * smooth
    #     H2[i, i + 1] = smooth
    # H[:N + 1, :N + 1] = H1
    # H3 = H2[1:-1]
    # H[N + 1:, :-1] = H3
    # h = np.zeros((Ntop + Nbot))
    # niter = 120
    # H = scsp.csr_matrix(H)
    # f = lambda x: wlstsq(x, G, H)
    # HH = H.T @ h
    # for iter in range(niter):
    #     print(iter)
    #     # mydiag=np.diag(Ap*jv(1,w*r/cp)*w*r*cp**(-2))
    #     # G1=scsp.csr_matrix((N,N+1))
    #     G1 = scsp.lil_matrix(np.diag(Ap * jv(1, wf * r / cp) * wf * r * cp ** (-2)))
    #     G2 = scsp.lil_matrix(jv(0, wf * r / c0)).transpose()
    #     G = scsp.hstack((G1, G2)).tocsr()
    #     ## delta rho, data vector
    #     Drho = data - rhop
    #     # E=Drho@Drho
    #     L = LinearOperator((G.shape[1], G.shape[1]), matvec=f, rmatvec=f)
    #     Dm = bicg(L, G.T @ Drho + HH, tol=1e-05, maxiter=4 * N)
    #     cp = cp + Dm[0][0:N]
    #     Ap = Ap + Dm[0][N]
    #     rhop = Ap * jv(0, wf * r / cp)
    #
    # fig, ax = plt.subplots(2, 1)
    # ax[0].plot(wf, c0, 'k', linewidth=2, label='grid search phase vel.')
    # ax[0].plot(wf, cp, 'r', linewidth=2, label='newton phase vel.')
    # #ax[0].plot(w, c_true, 'g--', linewidth=2, label='true phase vel.')
    # ax[0].set_xlabel('w')
    # ax[0].set_ylabel('c(w)')
    # ax[0].legend()
    #
    # ax[1].plot(wf, data, 'k--', linewidth=2, label='target corr func')
    # ax[1].plot(wf, rho0, 'g-', linewidth=2, label='estimated corr func grid search')
    # ax[1].plot(wf, rhop, 'r:', linewidth=2, label='estimated corr func newton')
    # ax[1].legend()
    #ax[1].set_xlabel('w')
    #plt.savefig('good_fit_2.png', format='png', dpi=300, bbox_inches='tight')
    # m=np.hstack((cp,A0))
    # xd=wlstsq(m,G,H)
    #print("--- %s seconds ---" % (time.time() - start_time))

    # K=3
    # L=40
    # w=np.linspace(2*np.pi*0.5,2*np.pi*40.0,K)
    #
    # c=np.linspace(40,800,L)
    # cv = [c[i:i + K] for i in range(0, len(c), 1) if len(c[i:i + K]) == K]
    # j_aju=np.zeros((wf.size,c.size))
    # dw=wf[1]-wf[0]
    # ND = [x for x in range(40)]
    #
    # combs=list(itertools.combinations(ND, 3))
    # j_aju=np.zeros((wf.size,len(combs)))
    # A=np.zeros(len(combs))
    # residual=np.zeros(len(combs))
    #
    # for j,comb in enumerate(combs):
    #         index=np.array(comb)
    #         #jj[:,j]=jv(0,r*2*np.pi*w/c[idx])
    #         #f=interpolate.interp1d(w,jj[:,j])
    #         jj=jv(0,r*w/c[index])
    #         f=interpolate.interp1d(w,jj)
    #         j_aju[:,j]=f(wf)
    #         #j_aju[:,j]=f(wf)
    #         A[j]=(j_aju[:,j].T@data/(j_aju[:,j].T@j_aju[:,j]))
    #         #residual[j]=integrate.simpson((A[j]*j_aju[:,j]-data)**2,x=wf,dx=dw)
    #         residual[j]=np.sum((j_aju[:,j]-data)**2)
    #
    # plt.plot(wf,data)
    # plt.plot(wf,j_aju[400])
    #xD=jv(0,10)
    #f1=interpolate.interp1d(w,jj[:,0])
    #j_aju=f1(wf)
    # r1=TS.datas_eff[0]
    # r2=TS.datas_eff[1]
    # size = int(512 * 60)
    # step = int(size / 2)
    # xd=[r1[i:i + size] for i in range(0, len(r1), step) if len(r1[i:i + size]) == size]
    # fs=512
    # s=10
    # size=fs*s
    # step=int(size/2)
    # W = tukey(size, alpha=0.1)
    # Pxy=csd(r1,r2,fs=fs,window=W,noverlap=step,detrend='linear',return_onesided=False)
    # Pxx=welch(r1,fs=fs,window=W,noverlap=step,detrend='linear',return_onesided=False)
    # Pyy=welch(r2,fs=fs,window=W,noverlap=step,detrend='linear',return_onesided=False)
    # Cxxy=(1/len(Pxy[0]))*np.real(Pxy[1])/np.sqrt(((1/len(Pxy[0]))**2)*np.real(Pxx[1])*np.real(Pyy[1]))
    # data_sp = np.split(Cxxy, 2)
    # # ## obtener la media del causar and acausal branches
    # data_stack = (data_sp[0] + data_sp[1][::-1]) / 2
    # #
    # Cxy,Cxy_acum,Rxy,Rxy_acum,lags=TS.preprocess_data_pilz(s=10,correl=True)
    # Cxy_final=Cxy_acum[:,-1,:]
    # Rxy_final=Rxy_acum[:,-1,:]
    # plt.plot(Cxy_final[0])
    # plt.plot(Cxxy)
    #plt.plot(Pxy[0], Cxy_final[0] / np.max(Cxy_final[0]))

    #plt.plot(fftshift(Pxy[0]), fftshift(np.real(Pxy[1]) / np.max(np.real(Pxy[1]))))
    #
    # wf,yf=TS.get_NCF(Cxy_final)
    # plt.figure(1)
    # plt.plot(wf,yf[0])
    # plt.figure(2)
    # plt.plot(data_stack)
    # dt=1./512
    # t=np.arange(0,10,dt)
    #tf,yf_r=TS.get_NCF_tspace(Rxy_final,t)
    #
    # r1=TS.datas_eff[0]
    # r2=TS.datas_eff[1]
    # r3=TS.datas_eff[2]
    # r4=TS.datas_eff[3]
    # r5=TS.datas_eff[4]
    #
    # xd1=correlate(r1,r2)
    # xd2=correlate(r1,r3)
    # xd3=correlate(r1,r4)
    # xd4=correlate(r1,r5)
    # xd5=correlate(r2,r3)
    # xd6=correlate(r2,r4)
    # xd7=correlate(r2,r5)
    # xd8=correlate(r3,r4)
    # xd9=correlate(r3,r5)
    # xd10=correlate(r4,r5)
    # idx=correlation_lags(r1.size,r2.size)
    # t=np.arange(-r1.size/512,r1.size/512,1./512)[1:]
    # tl1=t[np.argmax(xd1)]
    # tl2=t[np.argmax(xd2)]
    # tl3=t[np.argmax(xd3)]
    # tl4=t[np.argmax(xd4)]
    # tl5=t[np.argmax(xd5)]
    # tl6=t[np.argmax(xd6)]
    # tl7=t[np.argmax(xd7)]
    # tl8=t[np.argmax(xd8)]
    # tl9=t[np.argmax(xd9)]
    # tl10=t[np.argmax(xd10)]
    # plt.plot(t,xd1)
    # plt.plot(t,xd2)


    #TS.plot_spectrogram_pilz(Cxy_acum,lags/512,N=6)
    #TS.get_NCF(Cxy_final)

    # data=Cxy_acum[5,-1]
    # ff=fftfreq(data.size,d=1./512)
    # data_sp = np.split(data, 2)
    # # tp=lagst[len(lagst)//2:]
    # # ## obtener la media del causar and acausal branches
    # data_stack=(data_sp[0]+data_sp[1][::-1])/2
    # w=ff[:len(ff)//2]
    # idx=np.where(w<100.0)[0]
    # w=w[idx]
    # data_stack=data_stack[idx]
    # zc1=np.where(np.diff(np.sign(data_stack)))[0]
    # func = interpolate.interp1d(w, data_stack,kind='linear')
    # ww=np.linspace(w[0],w[-1],20000)
    # ynew=func(ww)
    # zc2=np.where(np.diff(np.sign(ynew)))[0]
    # xx=np.linspace(-1000,1000,1000)
    # J=jv(0.0,xx)
    # zc3=np.where(np.diff(np.sign(J)))[0]
    # zm=[zc3[i:i+len(zc2)] for i in range(0,len(zc3),2) if len(zc3[i:i+len(zc2)]) == len(zc2)]
    # #cm=np.zeros((len(zm),len(zm[0])))
    # #rho_rw=np.zeros((len(zm),len(zm[0])))
    # plt.plot(w, data_stack)
    # plt.plot(ww, ynew)
    # #plt.plot(w[zc1], data_stack[zc1], 'bo')
    #
    # plt.plot(ww[zc2], ynew[zc2], 'go')
    #for i in range(len(zm)):
        #cm[i]=2*np.pi*ww[zc2]*/(xx[zm[i]])
        #rho_rw[i]=jv(0.0,2*np.pi*ww[zc2]*22.8/cm[i])

    #plt.plot(ww,ynew)
    #plt.plot(ww[zc2], rho_rw[100], 'mo')
    #[deff[i:i + size] for i in range(0, len(deff), step) if len(self.datas_eff[0][i:i + size]) == size]
    #cm=np.zeros(len(xx))
    #for i in range(len(zc3)):
        #cm[i]=(zc2*22.8)/(zc3[i+m])
    #plt.plot(xx,J)


    #plt.plot(ww[zc3],J[zc3],'ro')

    #cwn=2*np.pi*ww[zc2[:9]]*22.8/ww[zc3[:9]]
    #j0=jv(0.0,2*np.pi*ww[zc2[0]]/(cwn[0]))
    # Cxy_acum,Rxy_acum,lags,Cxy_complex=TS.preprocess_data(s=20,correl=True)
    # lagst = lags / 512
    # idx = np.argwhere((lagst > -4) & (lagst < 4))
    # Cxy_acum_sh = Cxy_acum[:, :, idx][:, :, :, 0]
    # Rxy_acum_sh = Rxy_acum[:, :, idx][:, :, :, 0]
    # lagst_sh = lagst[idx][:, 0]
    #
    #
    # #TS.plot_spectrogram(Cxy_acum_sh, Rxy_acum_sh, lagst_sh, N=0)
    # #for i in range(10):
    # #    TS.plot_spectrogram(Cxy_acum,Rxy_acum,lagst,N=i)
    #
    # data=np.mean(Cxy_complex, axis=1)[0]
    # ff=fftfreq(data.shape[0],d=1./512)
    # fig,ax=plt.subplots(4,1)
    # ax[0].plot(fftshift(ff),fftshift(np.real(data)))
    # ax[1].plot(lagst,fftshift(np.real(ifft(data))))
    # data_sp=np.split(data,2)
    # tp=lagst[len(lagst)//2:]
    # ## obtener la media del causar and acausal branches
    # data_stack=(data_sp[0]+data_sp[1][::-1])/2
    # ax[2].plot(ff[:len(ff)//2],np.real(data_stack))
    # w=ff[:len(ff)//2]
    # wc = np.logspace(np.log10(0.1), np.log10(50), 300)
    # gamma = 1.0
    # alpha = gamma ** 2 * wc
    # F = np.zeros((len(wc), len(w)))
    # for i in range(len(wc)):
    #      F[i] = np.exp(-alpha[i] * ((w / wc[i]) - 1) ** 2)
    # fig,ax=plt.subplots(3,1)
    # CWF=data_stack*F
    # ax[0].plot(lagst,fftshift(np.real(ifft(data))))
    # ax[1].contourf(tp,wc,np.real(CWF))
    # ax[2].contourf(tp,wc,np.real(ifft(CWF)))

    # ax[0].plot(t, dataspi[0])
    # ax[1].contourf(t, wc, CWF)
    # datar=Rxy_acum[0,-1]
    #
    # datasp=np.split(data,2)
    # dataspi=np.split(ifft(data),2)
    # datat=(datasp[0] + datasp[1][::-1])/2
    # wc=np.logspace(np.log10(0.1),np.log10(250),300)
    # #wc2=-np.logspace(np.log10(0.1),np.log10(50),300)[::-1]
    # #wcc=np.hstack((wc2,wc))
    # gamma = 1.0
    # w=fftshift(fftfreq(data.shape[0],d=1./512))[data.shape[0]//2:]
    #
    # alpha = gamma ** 2 * wc
    # F = np.zeros((len(wc), len(w)))
    # #sosbp = butter(10, (0.5, 60), btype='bandpass', output='sos', fs=512)
    # #filt_bp = sosfilt(sosbp, np.real(datasp[0]))
    # #plt.plot(ifft(datasp[0]))
    # #plt.plot(datar)
    # for i in range(len(wc)):
    #     F[i] = np.exp(-alpha[i] * ((w / wc[i]) - 1) ** 2)
    # t = (lags / 512)[lags.size // 2:]
    # CWF=np.real(ifft(datat * F))
    # t = (lags / 512)[lags.size // 2:]
    # lagst = lags / 512
    #
    # fig, ax = plt.subplots(2, 1)
    # ax[0].plot(t, dataspi[0])
    # ax[1].contourf(t, wc, CWF)
    #ax[0].set_xlim([-0.1, 0.5])
    #ax[1].set_xlim([-0.1, 0.5])
    #plt.plot(filt_bp)
    #plt.plot(datasp[0])
    #plt.plot(datasp[1])
    # datareal=np.real(ifft(data))[:data.shape[0]//2]
    # data=fftshift(data[data.shape[0]//2:])
    # wc=np.logspace(np.log10(0.1),np.log10(80),300)
    # # #w=rfftfreq(data.shape[1],d=1/512)[1:]
    # gamma = 1.0
    # alpha = gamma ** 2 * wc
    # F=np.zeros((len(wc),len(w)))
    # sosbp = butter(10, (80,150.0), btype='bandpass', output='sos', fs=512)
    # filt_bp = sosfilt(sosbp, datareal)
    # plt.plot(datareal)
    # plt.plot(filt_bp)


    # fig,ax=plt.subplots(2,1)
    # ax[0].plot(t,datat)
    # ax[1].contourf(t,wc,CWF)
    # ax[0].set_xlim([-5, 30])
    # ax[1].set_xlim([-5, 30])

    # plt.contourf(t,wc,CWF)
    #
    # #Cxy_acum2,Rxy_acum2,lags2=TS.preprocess_data(s=60,correl=True)


    #fff=np.log10(np.abs(fft(Cxy_acum[0,-1])))**2
    #freq=fftfreq(fff.size,d=1/512)
    #plt.plot(fftshift(freq),fftshift(fff))
    #TS.plot_spectrogram(Cxy_acum, Rxy_acum, lagst, N=0)


    #fw1, CCw = welch(Cxy_acum[0,-1], fs=512, nperseg=512, detrend=False, return_onesided=False)
    #fw2, RRw = welch(Rxy_acum[0,-1], fs=512, nperseg=512, detrend=False, return_onesided=False)
    #plt.figure(1)
    #plt.plot(fftshift(fw1),fftshift(CCw))
    #plt.figure(2)
    #plt.plot(fftshift(fw2), fftshift(RRw))
    #plt.plot(lags[times]/512,Cxy_acum[0,-1,times].T)
    #times=np.where((lags/512>-10) & (10 > lags / 512))

    # fig,ax=plt.subplots(4,figsize=(12,10))
    # ax[0].specgram(Cxy_acum[0,-1],Fs=512,NFFT=256,noverlap=128,sides='twosided')
    # ax[0].set_ylim([-30,30])
    # ax[1].specgram(Rxy_acum[0,-1],Fs=512,NFFT=256,noverlap=64,sides='twosided')
    # ax[1].set_ylim([-30,30])
    # ax[2].specgram(Cxy_acum[4,-1],Fs=512,NFFT=64,noverlap=32,sides='twosided')
    # ax[2].set_ylim([-30,30])
    # ax[3].specgram(Rxy_acum[4,-1],Fs=512,NFFT=64,noverlap=32,sides='twosided')
    # ax[3].set_ylim([-30,30])
    # plt.xlabel('time[seconds]')
    # [x.set_ylabel('Frequency [Hz]') for x in ax]
    # #ax[1].specgram(Cxy_acum[0,-1],Fs=512,sides='onesided')
    #
    # f1, t1, Cxx = spectrogram(Cxy_acum[0, -1,times],nperseg=256, fs=512, detrend=False,return_onesided=False)
    # plt.figure(1)
    # plt.pcolormesh(t1-10, fftshift(f1), fftshift(np.log10(Cxx[0])), shading='gouraud')
    # plt.ylabel('Frequency [Hz]')
    # plt.xlabel('Time [sec]')
    #
    # f2, t2, Rxx = spectrogram(Rxy_acum[0, -1,times],nperseg=256, fs=512, detrend=False)
    #
    # spec=10*np.log10(2*np.abs(fft(Cxy_acum[0,-1,times]))**2/(len(times)))
    #
    # plt.figure(2)
    # plt.pcolormesh(t2-10, f2, np.log10(Rxx[0]), shading='gouraud')
    #
    # #plt.ylim([0,60])
    # #plt.xscale('log')
    # #plt.xlim([0,5])
    # plt.figure(3)
    # plt.plot(spec[0])
    # plt.xscale('log')
    ## media para todas las crosscoherencias
    #media=np.mean(Cxy,axis=1)

    #plt.plot(Rxy_acum[0,-1])
    #data1=np.mean(np.real(fftshift(ifft(Cxy[4]))),axis=0)
    #sosbp = butter(10, (1.0,15.0), btype='bandpass', output='sos', fs=512)
    #filt_bp = sosfilt(sosbp, data1)
    #plt.plot(filt_bp)
    #A=np.zeros((5,60,512*60),dtype='complex128')
    #size=512*60
    #W = tukey(size, alpha=0.25)

    #d1,d2,d3,d4,d5=TS.preprocess_data()

    # ss,ee,a,b,deff=TS.get_eff_signal()
    # #datasync=np.loadtxt(path+'/SynchronizedZ.dat').T
    # fig,ax=plt.subplots(1,1)
    # teff=np.arange(0,1859,1/512)
    # c=0
    # N=80000
    # for data,metadata in zip(deff,TS.metadatas):
    #     ax.plot(teff[:N],data[:N]+c*0.5,label=metadata.get('Instrument')+'_'+metadata.get('Med_and_Node'))
    #     ax.legend()
    #     c+=1
    #datass=[]
    #datass.append(deff[1])
    #datass.append(deff[3])
    #datass.append(deff[2])
    #datass.append(deff[0])
    #datass.append(deff[4])
    #datass=np.array(datass)

    #for d1,d2 in zip(datasync,datass):
    #    ax.plot(teff,d1+0.25*c)
    #    ax.plot(teff,d2+0.25*c)
    #    print(np.allclose(d1[:-1],d2[:-1]))
    #    c+=1
    #plt.plot(deff[0])
    #delta=np.min(ee)-np.max(ss)