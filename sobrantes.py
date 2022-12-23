 def plot_signals(self):
        pass

    def check_convergence(self, Cxy_acum, lags):
        pass

    def gaussian_filter(self, Cw,wc,w):
        gamma = 1.0
        alpha = gamma ** 2 * wc
        for i in range(len(wc[i])):
            F = np.exp(-alpha * ((w / wc)-1)**2)
            CWF=F*W

    def plot_spectrogram_pilz(self, Cxy_acum, lags, N=0):

        fig, ax = plt.subplots(6, 1, figsize=(16, 14))
        ax[0].plot(lags, Cxy_acum[N, -1], 'k')
        CC = 10 * np.log10(2 * np.abs(fft(Cxy_acum[N, -1]) ** 2) / (len(lags)))
        f = fftfreq(CC.size, d=1 / self.fs)
        ax[1].plot(fftshift(f), fftshift(CC), 'k')
        #ax[1].plot(fftshift(f), fftshift(RR), 'r')
        fw1, CCw = welch(Cxy_acum[N, -1], fs=self.fs, nperseg=512, detrend=False, return_onesided=False)

        ax[2].plot(fftshift(fw1), fftshift(CCw), 'k')
        print('cheesasse')

        ax[3].specgram(Cxy_acum[N, -1], Fs=512, NFFT=256, noverlap=0, sides='twosided')
        ax[3].set_ylim([-30, 30])

        sosbp = butter(10, (0.8, 40), btype='bandpass', output='sos', fs=self.fs)
        filt_bpc = sosfilt(sosbp, Cxy_acum[N, -1])
        ax[4].plot(lags, filt_bpc, 'k')

        ax[5].specgram(Cxy_acum[N, -1], Fs=512, NFFT=128, noverlap=64, sides='twosided')
        ax[5].set_ylim([-30, 30])

        path_new='/home/doctor/Doctor/Magister/Tesis/databases/process_data/corr_results_new/'
        if self.name_syntax == 1:
            drc=self.metadatas[0]['Med_and_Node'].split('_')[0]
        if self.name_syntax == 2:
            drc=self.metadatas[0]['Med_and_Node'].split('_')[0]+self.metadatas[0]['Med_and_Node'].split('_')[2]
        #plt.savefig(path_new+'win_60sc.png', format='png', dpi=300, bbox_inches='tight')
        plt.savefig(path_new+drc+'/result_'+str(N)+'.png', format='png', dpi=300, bbox_inches='tight')
        plt.close(fig)

    def plot_spectrogram(self, Cxy_acum, Rxy_acum, lags, N=0):

        fig, ax = plt.subplots(6, 2, figsize=(16, 14))
        ax[0, 0].plot(lags, Cxy_acum[N, -1], 'k')
        ax[0, 1].plot(lags, Rxy_acum[N, -1], 'r')
        CC = 10 * np.log10(2 * np.abs(fft(Cxy_acum[N, -1]) ** 2) / (len(lags)))
        RR = 10 * np.log10(2 * np.abs(fft(Rxy_acum[N, -1]) ** 2) / (len(lags)))
        f = fftfreq(CC.size, d=1 / self.fs)
        ax[1, 0].plot(fftshift(f), fftshift(CC), 'k')
        ax[1, 1].plot(fftshift(f), fftshift(RR), 'r')
        fw1, CCw = welch(Cxy_acum[N, -1], fs=self.fs, nperseg=512, detrend=False, return_onesided=False)
        fw2, RRw = welch(Rxy_acum[N, -1], fs=self.fs, nperseg=512, detrend=False, return_onesided=False)

        ax[2, 0].plot(fftshift(fw1), fftshift(CCw), 'k')
        ax[2, 1].plot(fftshift(fw2), fftshift(RRw), 'r')
        print('cheesasse')

        ax[3, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=256, noverlap=128, sides='twosided')
        ax[3, 0].set_ylim([-30, 30])
        ax[3, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=256, noverlap=128, sides='twosided')
        ax[3, 1].set_ylim([-30, 30])

        sosbp = butter(10, (0.8, 40), btype='bandpass', output='sos', fs=self.fs)
        filt_bpc = sosfilt(sosbp, Cxy_acum[N, -1])
        filt_bpr = sosfilt(sosbp, Rxy_acum[N, -1])
        ax[4, 0].plot(lags, filt_bpc, 'k')
        ax[4, 1].plot(lags, filt_bpr, 'r')

        ax[5, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=128, noverlap=64, sides='twosided')
        ax[5, 0].set_ylim([-30, 30])
        ax[5, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=128, noverlap=64, sides='twosided')
        ax[5, 1].set_ylim([-30, 30])

        path_new='/home/doctor/Doctor/Magister/Tesis/databases/process_data/corr_results/'
        drc=self.metadatas[0]['Med_and_Node'].split('_')[0]+self.metadatas[0]['Med_and_Node'].split('_')[2]
        #plt.savefig(path_new+'win_60sc.png', format='png', dpi=300, bbox_inches='tight')
        plt.savefig(path_new+drc+'/result_'+str(N)+'.png', format='png', dpi=300, bbox_inches='tight')
        plt.close(fig)


        # fig, ax = plt.subplots(6, 2, figsize=(16, 14))
        # ax[0, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=8, noverlap=4, sides='twosided')
        # ax[0, 0].set_ylim([-30, 30])
        # ax[0, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=8, noverlap=4, sides='twosided')
        # ax[0, 1].set_ylim([-30, 30])
        #
        # ax[1, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=32, noverlap=16, sides='twosided')
        # ax[1, 0].set_ylim([-30, 30])
        # ax[1, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=32, noverlap=16, sides='twosided')
        # ax[1, 1].set_ylim([-30, 30])
        #
        # ax[2, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=16, noverlap=8, sides='twosided')
        # ax[2, 0].set_ylim([-30, 30])
        # ax[2, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=16, noverlap=8, sides='twosided')
        # ax[2, 1].set_ylim([-30, 30])
        #
        # ax[3, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=64, noverlap=32, sides='twosided')
        # ax[3, 0].set_ylim([-30, 30])
        # ax[3, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=64, noverlap=32, sides='twosided')
        # ax[3, 1].set_ylim([-30, 30])
        #
        # ax[4, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=4, noverlap=0, sides='twosided')
        # ax[4, 0].set_ylim([-30, 30])
        # ax[4, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=4, noverlap=0, sides='twosided')
        # ax[4, 1].set_ylim([-30, 30])
        #
        # ax[5, 0].specgram(Cxy_acum[N, -1], Fs=512, NFFT=16, noverlap=0, sides='twosided')
        # ax[5, 0].set_ylim([-30, 30])
        # ax[5, 1].specgram(Rxy_acum[N, -1], Fs=512, NFFT=16, noverlap=0, sides='twosided')
        # ax[5, 1].set_ylim([-30, 30])
        #
        # plt.savefig(path_new+'win_8s_specgrams_nolap.png', format='png', dpi=300, bbox_inches='tight')


        # f1,t1,Cxx=spectrogram(Cxy_acum[N,-1],fs=self.fs,window='hamming',detrend=False, return_onesided=False)
        # f2,t2,Rxx=spectrogram(Rxy_acum[N,-1],fs=self.fs,window='hamming',detrend=False, return_onesided=False)
        #
        # fig,ax=plt.subplots()
        # plt.figure(1)
        # plt.pcolormesh(t1,fftshift(f1),fftshift(Cxx),shading='gouraud')
        #
        # plt.figure(2)
        # plt.pcolormesh(t2, fftshift(f2), fftshift(Rxx), shading='gouraud')
    def preprocess_data_pilz2(self,s=60,filter=False):
        if filter:
            soshp = butter(10, 0.8, btype='highpass', output='sos', fs=512)
            self.datas_eff = sosfilt(soshp, self.datas_eff)

        size = int(self.fs * s)
        step = int(size / 2)
        W = tukey(size, alpha=0.1)
        ND = [x for x in range(len(self.datas_eff))]

        a1 = self.datas_eff[0]
        d1 = [a1[i:i + size] for i in range(0, len(a1), step) if len(a1[i:i + size]) == size]
        nwins = len(d1)

        combs = list(itertools.combinations(ND, 2))
        print('ND', ND)
        print('combinaciones posibles', combs, len(combs))

        Cxy_complex = np.zeros((len(combs), nwins, size), dtype='complex128')
        Cxy = np.zeros((len(combs), size//2))
        for n,i in enumerate(combs):
            print(combs[0],combs[1])
            Pxy=csd(self.datas_eff[i[0]],self.datas_eff[i[1]],fs=self.fs,window=W,noverlap=step,detrend='linear',return_onesided=False)


            Pxx = welch(self.datas_eff[i[0]],fs=self.fs, window=W, noverlap=step, detrend='linear',return_onesided=False)


            Pyy = welch(self.datas_eff[i[1]], fs=self.fs, window=W, noverlap=step, detrend='linear', return_onesided=False)


            Cxyp = (1 / len(Pxy[0])) * np.real(Pxy[1]) / np.sqrt(((1 / len(Pxy[0])) ** 2) * np.real(Pxx[1]) * np.real(Pyy[1]))
            data_sp=np.split(Cxyp,2)
            Cxy[n]=data_sp[0] + data_sp[1][::-1] / 2
        return Cxy,Pxy[0],combs
## highpass filter to each segment
        #soshp = butter(10, 0.8, btype='highpass', output='sos', fs=self.fs)

    def preprocess_data_pilz(self,s=60,correl=False):
        size = int(self.fs * s)
        step = int(size / 2)
        dd = []
        fftdd = []
        ## apply 10% cosine taper
        W = tukey(size, alpha=0.1)
        ## highpass filter to each segment
        soshp = butter(10, 0.8, btype='highpass', output='sos', fs=self.fs)
        for deff in self.datas_eff:
            d = np.stack(
                [deff[i:i + size] for i in range(0, len(deff), step) if len(self.datas_eff[0][i:i + size]) == size])
            ## detrended and tapered data
            ddet = detrend(d)
            dtap=W*ddet
            ## apply the filter
            dhp=sosfilt(soshp,dtap)
            ## apply one bit normalization
            dsig=np.sign(dhp)
            #dsig=dhp
            print('dsig',dsig)
            ## a
            fftd = fft(dsig, norm='ortho')

            dd.append(dsig)
            ddd=np.stack(dd)
            fftdd.append(fftd)
            fftddd = np.stack(fftdd)


        nwins = len(d)
        print('NW', nwins)

        ND = [x for x in range(len(self.datas_eff))]

        combs = list(itertools.combinations(ND, 2))
        print('ND', ND)
        print('combinaciones posibles', combs, len(combs))
        c = 0
        #Cxy = np.zeros((len(combs), nwins, size))
        Cxy_complex = np.zeros((len(combs), nwins, size), dtype='complex128')
        Cxy = np.zeros((len(combs), nwins, size))
        Rxy = np.zeros((len(combs), nwins, size))

        print('CXY', Cxy.shape)
        #correl = True
        if not correl:
            for x in combs:
                Cxy[c] = np.real(fftddd[x[0]] * np.conj(fftddd[x[1]]))  / (np.abs(fftddd[x[0]]) * np.abs(fftddd[x[1]]))


                c += 1
        if correl:
            for x in combs:
                Cxy[c] = np.real(fftddd[x[0]] * np.conj(fftddd[x[1]]))  / (np.abs(fftddd[x[0]]) * np.abs(fftddd[x[1]]))
                Rxy[c]= correlate(ddd[x[0]],ddd[x[1]],mode='same')
                c += 1
        Cxy_acum = np.zeros((len(combs), nwins, size))
        Rxy_acum = np.zeros((len(combs), nwins, size))

        lags = correlation_lags(Cxy[0, 0].size, Cxy[0, 0].size, mode='same')
        for j in range(len(combs)):
            for i in range(1, 61):
                Cxy_acum[j] = np.cumsum(Cxy[j], axis=0) / i
                Rxy_acum[j] = np.cumsum(Rxy[j], axis=0) / i

        lags = correlation_lags(Cxy[0, 0].size, Cxy[0, 0].size, mode='same')
        if not correl:
            return Cxy,Cxy_acum,lags
        if correl:
            return Cxy,Cxy_acum,Rxy,Rxy_acum,lags
    def get_NCF_tspace(self,Rxyf,t):

        fig,axs =plt.subplots(5,2,figsize=(14,10))
        idx = np.where(t < 25.0)[0]
        tt= t[idx]

        for i in range(len(Rxyf)):
            data_sp = np.split(Rxyf[i], 2)
            # ## obtener la media del causar and acausal branches
            data_stack = (data_sp[0] + data_sp[1][::-1]) / 2
            print(data_stack.shape)
            #data_stack = data_stack[idx]

            if i<=4:
                axs[i,0].plot(data_stack)
                #axs[i,0].plot(ww[zcs],ynew[i][zcs],'ro')
            if i>4:
                axs[i-5,1].plot(data_stack)
                #axs[i-5,1].plot(ww,ynew[i][zcs],'ro')
            if i==4:
                axs[i,0].set_xlabel('Frequency [Hz]')
            if i == 9:
                axs[i-5, 1].set_xlabel('Frequency [Hz]')
        path_new = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/corr_results_new/'

        #plt.savefig(path_new+'/result_tspace_1.png', format='png', dpi=300, bbox_inches='tight')

        return tt, data_stack

        # return Cxy,Cxy_acum,Rxy,Rxy_acum,lags,ddd

        # return ddd,fftddd,Cxy

        # data_cut.append(dd)
        #
        # dd=[deff[i:i + size] for i in range(0, len(deff), step)]
        # data_cut.append(dd)
        # if len(self.metadatas)==5:
        #     d1=[self.datas_eff[0][i:i + size] for i in range(0, len(self.datas_eff[0]), step) if len(self.datas_eff[0][i:i + size])==size]
        #     d2=[self.datas_eff[1][i:i + size] for i in range(0, len(self.datas_eff[1]), step) if len(self.datas_eff[0][i:i + size])==size]
        #     d3=[self.datas_eff[2][i:i + size] for i in range(0, len(self.datas_eff[2]), step) if len(self.datas_eff[0][i:i + size])==size]
        #     d4=[self.datas_eff[3][i:i + size] for i in range(0, len(self.datas_eff[3]), step) if len(self.datas_eff[0][i:i + size])==size]
        #     d5=[self.datas_eff[4][i:i + size] for i in range(0, len(self.datas_eff[4]), step) if len(self.datas_eff[0][i:i + size])==size]
        #
        #
        #
        # if len(self.metadatas) == 4:
        #     d1 = [self.datas_eff[0][i:i + size] for i in range(0, len(self.datas_eff[0]), step) if len(self.datas_eff[0][i:i + size])==size]
        #     d2 = [self.datas_eff[1][i:i + size] for i in range(0, len(self.datas_eff[1]), step) if len(self.datas_eff[0][i:i + size])==size]
        #     d3 = [self.datas_eff[2][i:i + size] for i in range(0, len(self.datas_eff[2]), step) if len(self.datas_eff[0][i:i + size])==size]
        #     d4 = [self.datas_eff[3][i:i + size] for i in range(0, len(self.datas_eff[3]), step) if len(self.datas_eff[0][i:i + size])==size]
        #
        # #for deff in self.datas_eff:
        #     #dd=[deff[i:i + size] for i in range(0, len(deff), step)]
        #     #data_cut.append(dd)
        #
        #     #dd=[deff[i:i + size] for i in range(0, len(deff), step)]
        #     #data_cut.append(dd)
        # return d1,d2,d3,d4,d5

# ## hacer la timeseries con toda la medicion altiro y obtener la señal efectiva
# path1=path+'M10A_001part1/EqualizedFile.dat'
# path2=path+'M10A_002part1/EqualizedFile.dat'
# path3=path+'M10A_003part1/EqualizedFile.dat'
# path4=path+'M10A_004part1/EqualizedFile.dat'
# path5=path+'M10A_005part1/EqualizedFile.dat'
#
# TS1=TimeSeries(path1)
# TS2=TimeSeries(path2)
# TS3=TimeSeries(path3)
# TS4=TimeSeries(path4)
# TS5=TimeSeries(path5)
# class TimeSeries():
#     def __init__(self,path):
#         metadata=[]
#         with open(path, 'r', encoding='ISO-8859-1') as file:
#             c = 0
#             for line in file:
#                 if line[0:6] != '[mm/s]':
#                     tup = line.partition(':')
#
#                     metadata.append([x.strip() for x in tup])
#
#                     # print(line,c)
#                     c += 1
#                 else:
#                     #print(line, c)
#                     tup = line.partition(':')
#                     metadata.append([x.strip() for x in tup])
#                     c += 1
#                     break
#         dic = {}
#         dic.update({'Site_ID': metadata[0][2].split(',')[0]})
#         dic.update({'Med_and_Node': metadata[0][2].split(',')[1].strip()[:-3]})
#         dic.update({'Instrument': metadata[0][2].split(',')[1].strip()[-2:]})
#         dic.update({'Instrument_Tag': metadata[2][2]})
#         dic.update({'Data_Format_Bytes': int(metadata[4][2][:2])})
#         dic.update({'Full_scale_mV': int(metadata[5][2][:3])})
#         dic.update({'N_Channels': int(metadata[6][2])})
#         dic.update({'Sampling_Rate': int(metadata[7][2][:3])})
#         dic.update({'Start_Time': datetime.strptime(metadata[9][2], "%d/%m/%y %H:%M:%S")})
#         dic.update({'End_Time': datetime.strptime(metadata[10][2], "%d/%m/%y %H:%M:%S")})
#         dic.update({'Trace_length': metadata[11][2]})
#         if dic['N_Channels']==4:
#             print('dudabi')
#             dic.update({'Start_Recording_UTC': metadata[18][2].split('\t')[1]})  ## will check this later
#             dic.update({'Latitude': metadata[20][2]})
#             dic.update({'Longitude': metadata[21][2]})
#             dic.update({'Horizontal_Diluition': metadata[23][2]})
#             dic.update({'Geoid_Altitude': metadata[24][2]})
#             self.data=np.loadtxt(path, skiprows=c, encoding='ISO-8859-1',usecols=(0,1,2))
#         if dic['N_Channels']==8:
#             print('dubido')
#             dic.update({'Start_Recording_UTC': metadata[22][2].split('\t')[1]})  ## will check this later
#             dic.update({'Latitude': metadata[24][2]})
#             dic.update({'Longitude': metadata[25][2]})
#             dic.update({'Horizontal_Diluition': metadata[27][2]})
#             dic.update({'Geoid_Altitude': metadata[28][2]})
#             self.data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1',usecols=(0,1,2))
#         self.metadata=dic

    def get_NCFF(self,Cxyf):
        ff=fftfreq(Cxyf.shape[1],d=1./512)
        w = ff[:len(ff) // 2]
        ynew=np.zeros((len(Cxyf),20000))
        zcs=[]
        fig,axs =plt.subplots(5,2,figsize=(14,10))
        idx = np.where(w < 250.0)[0]
        w = w[idx]
        ww = np.linspace(w[0], w[-1], 20000)

        for i in range(len(Cxyf)):
            data_sp = np.split(Cxyf[i], 2)
            # ## obtener la media del causar and acausal branches
            data_stack = (data_sp[0] + data_sp[1][::-1]) / 2

            data_stack = data_stack[idx]
            ## interpolate to get a more precise zero crossing
            func = interpolate.interp1d(w, data_stack, kind='linear')

            ynew[i] = func(ww)

            zcs.append(np.where(np.diff(np.sign(ynew)))[0])
            if i<=4:
                axs[i,0].plot(ww,ynew[i])
                #axs[i,0].plot(ww[zcs],ynew[i][zcs],'ro')
            if i>4:
                axs[i-5,1].plot(ww,ynew[i])
                #axs[i-5,1].plot(ww,ynew[i][zcs],'ro')
            if i==4:
                axs[i,0].set_xlabel('Frequency [Hz]')
            if i == 9:
                axs[i-5, 1].set_xlabel('Frequency [Hz]')
        path_new = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/corr_results_new/'

        #plt.savefig(path_new+'/result_1.png', format='png', dpi=300, bbox_inches='tight')
        print( ww, ynew)
        return ww, ynew

    def preprocess_data(self, s=60, correl=False):
        ## hacer el procesamiento de los datos aqui, leer las trazas, detrendearlas, calcular fft, y crosscoherencia, hacer un diccionario con las crosscoherencias
        ## predefined size=60seconds and 50% overlap (30s step)
        data_cut = []
        print('q')
        size = int(self.fs * s)
        step = int(size / 2)
        dd = []
        fftdd = []
        W = tukey(size, alpha=0.25)

        print(self.datas_eff[0])
        for deff in self.datas_eff:
            d = np.stack(
                [deff[i:i + size] for i in range(0, len(deff), step) if len(self.datas_eff[0][i:i + size]) == size])
            print(len(d), 'd')
            dp = W * detrend(d)
            fftd = fft(dp, norm='ortho')
            dd.append(dp)
            ddd = np.stack(dd)
            fftdd.append(fftd)
            fftddd = np.stack(fftdd)
        nwins = len(d)
        print('NW', nwins)

        ND = [x for x in range(len(self.datas_eff))]

        combs = list(itertools.combinations(ND, 2))
        print('ND', ND)
        print('combinaciones posibles', combs, len(combs))
        Rxy = np.zeros((len(combs), nwins, size))
        Rxy2 = np.zeros((len(combs), nwins, size))

        c = 0
        Cxy = np.zeros((len(combs), nwins, size))
        Cxy_complex = np.zeros((len(combs), nwins, size),dtype='complex128')
        print('CXY', Cxy.shape)
        for x in combs:
            Cxy_complex[c] = fftddd[x[0]] * np.conj(fftddd[x[1]]) / (np.abs(fftddd[x[0]]) * np.abs(fftddd[x[1]]))
            Cxy[c] = fftshift(np.real(ifft(Cxy_complex[c], norm='ortho')))
            print(Cxy.shape, c, self.metadatas[x[0]].get('Node_Instrument'),
                  self.metadatas[x[1]].get('Node_Instrument'))
            if correl:
                Rxy[c] = correlate(ddd[x[0]], ddd[x[1]], mode='same')
                # Rxy2[c]=correlate(ddd[x[0]],ddd[x[1]],mode='full')

            c += 1
        Cxy_acum = np.zeros((len(combs), nwins, size))
        Rxy_acum = np.zeros((len(combs), nwins, size))
        lags = correlation_lags(Cxy[0, 0].size, Cxy[0, 0].size, mode='same')
        for j in range(len(combs)):
            for i in range(1, 61):
                Cxy_acum[j] = np.cumsum(Cxy[j], axis=0) / i
                Rxy_acum[j] = np.cumsum(Rxy[j], axis=0) / i

        return Cxy_acum, Rxy_acum, lags,Cxy_complex
        
    def correlations2(self, figure=False):
        if figure:
            fig, ax = plt.subplots(len(self.combs) // 2, 2)
        rijs = []
        tmaxes = []
        ts = []
        for n, x in enumerate(self.combs):
            rij = correlate(self.datas[x[0]], self.datas[x[1]])
            idx = correlation_lags(self.datas[0].size, self.datas[1].size)
            t = np.arange(-self.datas[0].size / 512, self.datas[1].size / 512, 1. / 512)[1:]
            print('wow', t.shape, rij.shape)
            tmax = t[np.argmax(rij)]
            # rijs.append(rij)
            tmaxes.append(tmax)
            # ts.append(t)
            if figure:
                if n <= 4:
                    ax[n, 0].plot(t, rij, 'r')
                    ax[n, 0].plot(tmax, np.max(rij), 'go')
                if n > 4:
                    ax[n - 5, 1].plot(t, rij, 'r')
                    ax[n - 5, 1].plot(tmax, np.max(rij), 'go')
        problematic_stations=[]
        for n,x in enumerate(tmaxes):
            ## when x>2.5, the problematic station is the first one that enters the correlation function
            ## when x<-2.5, the problematic station is the second one that enters the correlation function
            if x>2.5:
                problematic_stations.append(self.combs[n][0])
                #data_new=self.datas_eff[n][512 * 3: -1]

            if x <-2.5:
                problematic_stations.append(self.combs[n][1])

        return tmaxes,problematic_stations

    def manage_data(paths):
    datas = []
    metadatas = []
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
        # da, mda = self.create_data(path, metadata, c)

        da, mda = create_data(path, metadata, c)
        print('bob')
        datas.append(da)
        metadatas.append(mda)
    return datas,metadatas

def manage_data_new(paths):
    datas = []
    metadatas = []
    for path in paths:
        metadata = []
        with open(path, 'r', encoding='ISO-8859-1') as file:
            print('lol')
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
        #print('bub')
        # da, mda = self.create_data(path, metadata, c)

        da, mda = create_data_new(path, metadata, c)
        datas.append(da)
        metadatas.append(mda)
    return datas,metadatas

# def create_data(self, path, metadata, c, max=False):
    #     dic = {}
    #     dic.update({'Site_ID': metadata[0][2].split(',')[0]})
    #     dic.update({'Med_and_Node': metadata[0][2].split(',')[1].strip()[:-3]})
    #     dic.update({'Instrument': metadata[0][2].split(',')[1].strip()[-2:]})
    #     dic.update({'Instrument_Tag': metadata[2][2]})
    #     dic.update({'Node_Instrument': metadata[0][2].split(',')[1].split('_')[1] + '_' + metadata[0][2].split(',')[
    #                                                                                           1].strip()[-2:]})
    #     dic.update({'Data_Format_Bytes': int(metadata[4][2][:2])})
    #     dic.update({'Full_scale_mV': int(metadata[5][2][:3])})
    #     dic.update({'N_Channels': int(metadata[6][2])})
    #     dic.update({'Sampling_Rate': int(metadata[7][2][:3])})
    #     dic.update({'Start_Time': datetime.strptime(metadata[9][2], "%d/%m/%y %H:%M:%S")})
    #     dic.update({'End_Time': datetime.strptime(metadata[10][2], "%d/%m/%y %H:%M:%S")})
    #     dic.update({'Trace_length': metadata[11][2]})
    #     if dic['N_Channels'] == 4:
    #         print('channel 4')
    #         dic.update({'Start_Recording_UTC': metadata[18][2].split('\t')[1]})  ## will check this later
    #         dic.update({'Latitude': metadata[20][2]})
    #         dic.update({'Longitude': metadata[21][2]})
    #         dic.update({'Horizontal_Diluition': metadata[23][2]})
    #         dic.update({'Geoid_Altitude': metadata[24][2]})
    #         if not max:
    #             #print('q')
    #             data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
    #         if max:
    #             #print('qa')
    #             data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)
    #     if dic['N_Channels'] == 7:
    #         print('channel 7')
    #         dic.update({'Start_Recording_UTC': metadata[21][2].split('\t')[1]})  ## will check this later
    #         dic.update({'Latitude': metadata[23][2]})
    #         dic.update({'Longitude': metadata[24][2]})
    #         dic.update({'Horizontal_Diluition': metadata[26][2]})
    #         dic.update({'Geoid_Altitude': metadata[27][2]})
    #         if not max:
    #             data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
    #         if max:
    #             data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)
    #     if dic['N_Channels'] == 8:
    #         print('channel 8')
    #         dic.update({'Start_Recording_UTC': metadata[22][2].split('\t')[1]})  ## will check this later
    #         dic.update({'Latitude': metadata[24][2]})
    #         dic.update({'Longitude': metadata[25][2]})
    #         dic.update({'Horizontal_Diluition': metadata[27][2]})
    #         dic.update({'Geoid_Altitude': metadata[28][2]})
    #         if not max:
    #             data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
    #         if max:
    #             data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)
    #
    #
    #     mdata_dic = dic
    #     return data, mdata_dic



    Hello Matt, I have a question that has been bugging me since a while.

    I attach some noise correlation functions of the same data obtained by averaging segments of the complete time series.

    It is clear to me that smaller windows leads to more averaging => less data and less noise, as we can see in the pictures. Therefore using small windows (e.g 1s is the smaller I show here) would lead to better fits of the bessel function. 

    However,  if I use so small windows I would lose information in the lower frequencies, which iseems to not change the shape of the NCF at all, and I'm not sure 