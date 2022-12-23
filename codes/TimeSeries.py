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

def tanhfit(omega,d,e,f):
    c=d*np.arctanh(e*omega)+f/np.sqrt(omega)
    return omega
def get_nodes_new(cgrid, NP, NR, iterable):
    print(NP, NR)
    data = np.zeros((NP, NR))
    for n, x in enumerate(iterable):
        print('a',x)
        for i in range(len(x)):
            data[n, i] = cgrid[x[i]]
    return data
def preprocess_data(datas,fs,s,filter,combs):
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


def manage_data_new2(path):

    metadata = []
    with open(path, 'r', encoding='ISO-8859-1') as file:
        c = 0
        for line in file:
            if line[0:6] != '[mm/s]':
                tup = line.partition(':')

                metadata.append([x.strip() for x in tup])

                print(line,c)
                c += 1
            else:

                print(line, c)
                tup = line.partition(':')
                metadata.append([x.strip() for x in tup])
                c += 1
                break
                # print('c', c)
    print('alpaca',path)
    # da, mda = self.create_data(path, metadata, c)

    #try:
    d, mda = create_data_new(path, metadata, c)
    #except:
    #    print('failure')
    #    d, mda = create_data_new_patch(path, metadata, c)

    return d,mda



def create_data_new(path, metadata, c, max=False):
    #print(metadata[0],type(metadata[0]))
    dic = {}
    dic.update({'Site_ID': metadata[0][2].split(',')[0]})
    dic.update({'Med_and_Node': metadata[0][2].split(',')[1].strip()[:-3]})
    dic.update({'Instrument': metadata[0][2].split(',')[1].strip()[-2:]})
    dic.update({'Instrument_Tag': metadata[2][2]})
    #print('as')
    #dic.update({'Node_Instrument': metadata[0][2].split(',')[1].split('_')[1] + '_' + metadata[0][2].split(',')[
    #                                                                                      1].strip()[-2:]})
    dic.update({'Node':metadata[0][2].split(',')[1].split('_')[1].split(' ')[0]})
    dic.update({'Instrument':metadata[0][2].split(' ')[2]})
    dic.update({'Data_Format_Bytes': int(metadata[4][2][:2])})
    print('XD',metadata[5][2])
    dic.update({'Full_scale_mV': int(metadata[5][2][:3])})
    dic.update({'N_Channels': int(metadata[6][2])})
    dic.update({'Sampling_Rate': int(metadata[7][2][:3])})
    dic.update({'Start_Time': datetime.strptime(metadata[9][2], "%d/%m/%y %H:%M:%S")})
    dic.update({'End_Time': datetime.strptime(metadata[10][2], "%d/%m/%y %H:%M:%S")})
    dic.update({'Trace_length': metadata[11][2]})
    if dic['N_Channels'] == 4:
        #print('channel 4')
        dic.update({'Start_Recording_UTC': metadata[18][2].split('\t')[1]})  ## will check this later
        dic.update({'Latitude': metadata[20][2]})
        dic.update({'Longitude': metadata[21][2]})
        dic.update({'Horizontal_Diluition': metadata[23][2]})
        dic.update({'Geoid_Altitude': metadata[24][2]})
        if not max:
            # print('q')
            data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
            print('loaded')
        if max:
            # print('qa')
            data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)
    if dic['N_Channels'] == 7:
        #print('channel 7')
        dic.update({'Start_Recording_UTC': metadata[21][2].split('\t')[1]})  ## will check this later
        dic.update({'Latitude': metadata[23][2]})
        dic.update({'Longitude': metadata[24][2]})
        dic.update({'Horizontal_Diluition': metadata[26][2]})
        dic.update({'Geoid_Altitude': metadata[27][2]})
        if not max:
            data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
        if max:
            data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)
    if dic['N_Channels'] == 8:
        #print('channel 8')
        dic.update({'Start_Recording_UTC': metadata[22][2].split('\t')[1]})  ## will check this later
        dic.update({'Latitude': metadata[24][2]})
        dic.update({'Longitude': metadata[25][2]})
        dic.update({'Horizontal_Diluition': metadata[27][2]})
        dic.update({'Geoid_Altitude': metadata[28][2]})
        if not max:
            data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
            #print('loaded')

        if max:
            data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)

    mdata_dic = dic
    return data, mdata_dic


def dms_to_dd(d, m, s):
    dd = d + float(m)/60 + float(s)/3600
    return dd
def zero_crossings(w,rho_obs,r):
    N=10000
    ww=np.linspace(w[0],w[-1],N)
    ip = interpolate.interp1d(w, rho_obs, 'linear', fill_value='extrapolate')
    rho_interp = ip(ww)
    zc = np.where(np.diff(np.sign(rho_interp)))[0]
    pos_zc=zc[0::2]
    neg_zc=zc[1::2]
    wpos=ww[pos_zc]
    wneg=ww[neg_zc]
    wtod=ww[zc]
    zeros = jn_zeros(0, 2000)
    cc = np.zeros((25, zc.size))
    for m in range(25):
        print(2 * m, 2 * m + 12)
        cc[m, :] = wtod * 200 / (zeros[m:m + zc.size])

    fig,ax=plt.subplots(2,1)
    f=wtod/(2*np.pi)
    ax[0].plot(f, cc[0, :])
    ax[0].plot(f, cc[1, :])
    ax[0].plot(f, cc[2, :])
    ax[0].plot(f, cc[3, :])
    ax[0].plot(f, cc[4, :])

    ax[1].plot(ww/(2*np.pi),rho_interp)
    ax[1].plot(ww[zc]/(2*np.pi),rho_interp[zc],'go')
    return wtod,cc,zc
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

def newton_search(w,rho_obs,rho0,c0,A0,r):
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
    small = 1
    H1 = small * np.eye(Ntop)
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
    create_figures=True
    if create_figures:
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(w/(2*np.pi), c0, 'k', linewidth=2, label='grid search phase vel.')
        ax[0].plot(w/(2*np.pi), cp, 'r', linewidth=2, label='newton phase vel.')
        # ax[0].plot(w, c_true, 'g--', linewidth=2, label='true phase vel.')
        ax[0].set_xlabel('w')
        ax[0].set_ylabel('c(w)')
        ax[0].legend()
    
        ax[1].plot(w/(2*np.pi), rho_obs, 'k--', linewidth=2, label='target corr func')
        ax[1].plot(w/(2*np.pi), rho0, 'g-', linewidth=2, label='estimated corr func grid search')
        ax[1].plot(w/(2*np.pi), rhop, 'r:', linewidth=2, label='estimated corr func newton')
        ax[1].legend()
        ax[1].set_xlabel('w')
        plt.savefig('ajuste10.png', format='png', dpi=300, bbox_inches='tight')
    return rhop,cp,Ap
def get_nodes(cgrid,NG,NC):
    ns = np.arange(0, NC)
    perms = [np.array(p) for p in product(ns, repeat=NG)]
    ## all possible phase velocities
    c_nodes = np.array([[cgrid[x[y]] for y in range(NG)] for x in perms])
    ## split phase velocity vector to reduce size array an avoid memory errors
    c_nodes_sp = np.split(c_nodes, 256)
    return c_nodes_sp
def grid_search_cw2(N,wmin,wmax,rho_obs,r,NN,NC,cbounds,ccombs,create_figures):
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

    return A_best,rho_best,err_best,ccombs_best,c_best,wgrid,w
def grid_search_cw(N,wmin,wmax,rho_obs,r,NN,NC,cbounds,create_figures):
    ## NN: NUMERO DE VELOCIDADES
    ## NC: NUMERO DE NODOS
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
    ## compute all possible permutations with repetitions of subsets with NG elements in NC set
    print('bar')
    #A_best, rho_best, err_best, cnodes_best, c_best = pool.apply(ts.grid_search_cw,args=(N, wmin, wmax, data, r, 4, 40, cbounds, False))
    toc = time.time()
    print('foo')
    ns=np.arange(0,NC)
    perms=[np.array(p) for p in product(ns,repeat=NN)]
    ## all possible phase velocities
    c_nodes = np.array([[cgrid[x[y]] for y in range(NN)] for x in perms])
    ## split phase velocity vector to reduce size array an avoid memory errors
    c_nodes_sp = np.split(c_nodes, 100)
    print('bar')

    ## create a precomputed version of J0
    NJ=1001
    jargmin=0.5*wmin*r/chigh
    jargmax=2.0*wmax*r/clow
    LJ=jargmax-jargmin
    Dj=(NJ-1)/LJ
    jarg = np.linspace(jargmin, jargmax, NJ)
    nu = 0
    J0 = jv(nu, jarg)

    for n,x in enumerate(c_nodes_sp):
        print(n,len(c_nodes_sp))
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
            cnodes_best = c_nodes_sp[n][idx]
            c_best = ip(w)[idx]
        elif n > 0 and np.min(err) < err_best:
            # print('foo',np.min(err),cnodes_best)
            A_best = A_est[idx]
            rho_best = rho_pre[idx]
            err_best = err[idx]
            cnodes_best = c_nodes_sp[n][idx]
            c_best = ip(w)[idx]
    if create_figures:
        fig,ax=plt.subplots(2,1)
        #ax[0].plot(w, c_obs, 'k', linewidth=2, label='target phase vel')
        ax[0].plot(w, c_best, 'r', linewidth=2, label='est phase vel')
        ax[0].plot(wgrid, cnodes_best, 'go', linewidth=2)
        ax[0].set_xlabel('w')
        ax[0].set_ylabel('c(w)')
        ax[0].legend()

        ax[1].plot(w, rho_obs, 'k--', linewidth=2, label='target corr func')
        ax[1].plot(w, rho_best, 'r:', linewidth=2, label='estimated corr func')
        ax[1].legend()
        ax[1].set_xlabel('w')

    return A_best,rho_best,err_best,cnodes_best,c_best
    ## Set up
class TimeSeries():
    def __init__(self,paths,mp=False,create_data=True):
        if create_data:
            self.datas = []
            self.metadatas = []
            self.fs = 512
            self.name_syntax=1
            self.s=1

            tic = time.time()
            pool = multiprocessing.Pool()
            tic = time.time()
            print('empieza 1')
            R = pool.map(manage_data_new2, paths)
            pool.close()
            pool.join()
            self.datas=[R[x][0] for x in range(len(R))]
            self.metadatas=[R[x][1] for x in range(len(R))]
            toc = time.time()

            print('loaded data and metadata in {:.4f} seconds'.format(toc - tic))

            #datas_raw = [R[x][0] for x in range(len(R))]
            #metadatas_raw = [R[x][1] for x in range(len(R))]

            tic = time.time()
            N=len(self.datas)
            self.combs = self.get_combs(N)
            #self.datas_eff,self.t_eff,self.tv=self.get_eff_signal()
            datas_eff,t_eff,tv=self.get_eff_signal()

            toc = time.time()
            print('got the raw effective  data in {:.4f} seconds'.format(toc - tic))

            tic = time.time()

            self.tmaxes,self.problems,self.datas_eff,self.t_eff,self.tv=self.correlations(datas_eff,tv)

            toc = time.time()
            if len(self.problems)>0:
                print('fixed the incorrect data in {:.4f} seconds'.format(toc - tic))
            else:
                print('no problematic stations were found')

            pool = multiprocessing.Pool()
            tic = time.time()

            out=pool.map(partial(preprocess_data,self.datas_eff,self.fs,self.s,True),self.combs)
            self.C=np.array([x[0] for x in out])
            self.w=np.split(out[0][1],2)[0]
            self.coords,self.distances=self.get_coordinates(paths)
            pool.close()
            pool.join()
            toc = time.time()
            print(' got the correlation functions in in {:.4f} seconds'.format(toc - tic))

            tic = time.time()
            self.save_data_and_metadata(paths)
            self.create_NCF_data()
            self.save_NCF_data(paths)
            toc = time.time()
            print(' saved output in {:.4f} seconds'.format(toc - tic))

        if not create_data:
            #print('zz',paths)
            tic = time.time()

            self.datas_eff=np.loadtxt(paths+'_TS.txt')
            self.metadatas=read_csv(paths+'_MD.csv')
            self.w=np.loadtxt(paths+'_NCF_DATA.txt',max_rows=1)
            self.C=np.loadtxt(paths+'_NCF_DATA.txt',skiprows=1)
            self.md_ncf=read_csv(paths+'_NCF_MD.csv')
            toc = time.time()
            print(' loaded data in {:.4f} seconds'.format(toc - tic))


    def create_NCF_data(self):
        metadatas=[]
        for n,x in enumerate(self.combs):
            dic={}
            dic.update({'Combinations':x})
            dic.update({'Interstation_Distance':self.distances[n]})
            metadatas.append(dic)
        return metadatas
    def save_NCF_data(self,paths):
        meta_ncf=self.create_NCF_data()
        if 'SUR1A' in paths[0]:
            folder='/SUR1A'
        if 'SUR1B' in paths[0]:
            folder='/SUR1B'
        if 'SUR2A' in paths[0]:
            folder = '/SUR2A'
        if 'SUR2B' in paths[0]:
            folder='/SUR2B'

        med=paths[0].split('/')[-3]
        path = os.getcwd() + folder + '/Outputs/' +med + '_'
        #path=os.getcwd()+folder+'/Outputs/'+self.metadatas[0].get('Med_and_Node').split('_')[0]+'_'
        np.savetxt(path+'NCF_DATA.txt',np.vstack((self.w,self.C)))
        DataFrame(meta_ncf).to_csv(path+'NCF_MD.csv')
    def save_data_and_metadata(self,paths):
        if 'SUR1A' in paths[0]:
            folder='/SUR1A'
        if 'SUR1B' in paths[0]:
            folder='/SUR1B'
        if 'SUR2A' in paths[0]:
            folder = '/SUR2A'
        if 'SUR2B' in paths[0]:
            folder='/SUR2B'
        med=paths[0].split('/')[-3]
        path=os.getcwd()+folder+'/Outputs/'+med+'_'

        #path=os.getcwd()+folder+'/Outputs/'+self.metadatas[0].get('Med_and_Node').split('_')[0]+'_'
        np.savetxt(path+'TS.txt',self.datas_eff.T)
        DataFrame(self.metadatas).to_csv(path+'MD.csv')

    def get_coordinates(self,paths):
        print(paths[0])
        if 'SUR1A' in paths[0]:
            coords=np.loadtxt('coords_SUR1A.txt')
        if 'SUR1B' in paths[0]:
            coords=np.loadtxt('coords_SUR1B.txt')
        if 'SUR2A' in paths[0]:
            coords=np.loadtxt('coords_SUR2A.txt')
        if 'SUR2B' in paths[0]:
            coords=np.loadtxt('coords_SUR2B.txt')
        print('metadatas',self.metadatas)
        for x in self.metadatas:
                print('x',x)
                #print('a',x)

                #punto=int(x.get('Node')[1:])
                #punto=int(x.get('Node')[1:])
                punto=int(x.get('Node'))

                #print(punto)
                #print(coords[punto-1])
                x.update({'Local_Coordinates':coords[punto-1]})
        distances=[]
        for x in self.combs:
                 d=np.linalg.norm(self.metadatas[x[0]].get('Local_Coordinates')-self.metadatas[x[1]].get('Local_Coordinates'))
                 distances.append(d)

        return coords,distances


    def correlations(self,datas_eff,tv,figure=False):
        if figure:
            fig,ax=plt.subplots(len(self.combs)//2,2)
        rijs=[]
        tmaxes=[]
        ts=[]
        problematic_stations = []
        for n,x in enumerate(self.combs):
            print('n',n)
            rij=correlate(datas_eff[x[0]],datas_eff[x[1]])
            idx = correlation_lags(datas_eff[x[0]].size, datas_eff[x[1]].size)
            t = np.arange(-datas_eff[x[0]].size / 512, datas_eff[x[1]].size / 512, 1. / 512)[1:]
            #print('wow',t.shape,rij.shape)
            tmax=t[np.argmax(rij)]
            if tmax > 2.5:
                problematic_stations.append(self.combs[n][0])
            if tmax < -2.5:
                problematic_stations.append(self.combs[n][1])

            #rijs.append(rij)
            tmaxes.append(tmax)
            #ts.append(t)
            if figure:
                if n <= 4:
                    ax[n,0].plot(t,rij,'r')
                    ax[n,0].plot(tmax,np.max(rij),'go')
                if n>4:
                    ax[n-5, 1].plot(t, rij, 'r')
                    ax[n-5, 1].plot(tmax, np.max(rij), 'go')
        problematic_stations=list(set(problematic_stations))
        dataeff_new=[]
        for n in range(len(datas_eff)):
            #print('n',n)
            if n in problematic_stations:
                 dataeff_new.append(datas_eff[n,512*3:-1])

                 #self.datas_eff[n]=self.datas_eff[n][512 * 3: -1]
            if n not in problematic_stations:
                 dataeff_new.append(datas_eff[n,:-1-512*3])
        dataeff_new=np.array(dataeff_new)
        tv_new=tv[0:-1-512*3]
        teff_new=tv[-1]-tv[0]

        return tmaxes, problematic_stations,dataeff_new,tv_new,teff_new
        #return tmaxes
        #return rijs,tmaxes,ts
    def manage_data(self,paths):
            datas=[]
            metadatas=[]
            for path in paths:
                    metadata = []
                    with open(path, 'r', encoding='ISO-8859-1') as file:
                        c = 0
                        for line in file:
                            if line[0:6] != '[mm/s]':
                                tup = line.partition(':')

                                metadata.append([x.strip() for x in tup])

                                # print(line,c)
                                c += 1
                            else:

                                # print(line, c)
                                tup = line.partition(':')
                                metadata.append([x.strip() for x in tup])
                                c += 1
                                break
                        # print('c', c)
                    print('bub')
                    #da, mda = self.create_data(path, metadata, c)

                    da, mda = create_data(path,metadata,c)
                    print('bob')
                    self.datas.append(da)
                    self.metadatas.append(mda)
            print('du')
            self.datas_eff, self.t_eff, self.tt = self.get_eff_signal()
            print('da')

    def plot_traces(self,name,closefig=True):
        rows = len(self.datas_eff)
        fig, ax = plt.subplots(rows,figsize=(14, 10))
        print('lol',ax.shape)
        for i, data in enumerate(self.datas_eff):
            ax[i].plot(self.t_eff,data)
                #ax[i, 0].text(60, 0.4, 'ID:' + str(np.round(self.md_ncf['Interstation_Distance'][i], 2)) + '[m]')
                # ax[i,0].set_xlabel('ID:'+str(np.round(self.md_ncf['Interstation_Distance'][i],2)))
            #if i >= rows:
                #ax[i - rows, 1].plot(self.t_eff, data)
                #ax[i - rows, 1].text(60, 0.4, 'ID:' + str(np.round(self.md_ncf['Interstation_Distance'][i], 2)) + '[m]')

        ax[-1].set_xlabel('Time[s]')
        #ax[-1, 1].set_xlabel('Frequency (Hz)')
        # path = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/SUR2B/'
        path = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/' + name
        print('q',path)
        plt.savefig(path, format='png', dpi=300, bbox_inches='tight')
        if closefig:
            plt.close()
    def get_eff_signal(self):
        stime = [x.get('Start_Time') for x in self.metadatas]
        etime = [x.get('End_Time') for x in self.metadatas]
        teff_s = np.max(stime)
        teff_e = np.min(etime)
        ts_i = [(teff_s - x).seconds for x in stime]
        tf_i = [(x - teff_e).seconds for x in etime]
        t_eff = (teff_e - teff_s).seconds
        data_eff = [x[self.fs * y:self.fs * (y + t_eff)] for x, y in zip(self.datas, ts_i)]
        data_eff = np.array(data_eff)
        print(t_eff)
        dt=1/self.fs
        tv=np.arange(0,t_eff,dt)
        return data_eff, t_eff, tv
    def get_combs(self,N):
        ND = [x for x in range(N)]
        combs = list(itertools.combinations(ND, 2))
        print('ND', ND)
        print('combinaciones posibles', combs, len(combs))
        return combs

    def get_NCF(self,Cxy,w,name,closefig=True):
        wp = w
        ## analyze frequencies lower than 50 Hz
        idx = np.where((wp>0.5)&(wp < 75))[0]
        wp=wp[idx]
        Cxy=Cxy[:,idx]
        rows=len(Cxy)//2
        if rows>1:
            print('d',rows,len(Cxy))
            fig,ax=plt.subplots(rows,2,figsize=(14,10))
            for i,coh in enumerate(Cxy):
                if i < rows:
                    ax[i,0].plot(wp,coh)
                    ax[i,0].text(60,0.4,'ID:'+str(np.round(self.md_ncf['Interstation_Distance'][i],2))+'[m]')
                    #ax[i,0].set_xlabel('ID:'+str(np.round(self.md_ncf['Interstation_Distance'][i],2)))
                if i >= rows:
                    ax[i-rows,1].plot(wp,coh)
                    ax[i-rows,1].text(60,0.4,'ID:'+str(np.round(self.md_ncf['Interstation_Distance'][i],2))+'[m]')

            ax[-1,0].set_xlabel('Frequency (Hz)')
            ax[-1,1].set_xlabel('Frequency (Hz)')
        else:
            fig,ax=plt.subplots(len(Cxy),figsize=(14,10))
            for i,coh in enumerate(Cxy):
                ax[i].plot(wp, coh)
                ax[i].text(60, 0.4, 'ID:' + str(np.round(self.md_ncf['Interstation_Distance'][i], 2)) + '[m]')
        #path = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/SUR2B/'
        path = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/'+name
        print('k',path)

        plt.savefig(path,format='png', dpi=300, bbox_inches='tight')
        if closefig:
            plt.close()
        # med = 'M8B'
        # plt.savefig(path + '/NCF_' + med + '.png', format='png', dpi=300, bbox_inches='tight')

        #plt.savefig(path + '/NCF_' + med + '.png', format='png', dpi=300, bbox_inches='tight')
        #plt.plot(wp[idx], Cxy[0][idx])

