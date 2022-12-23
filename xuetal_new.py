import numpy as np
from scipy.signal import detrend,hann,coherence,correlate,periodogram,correlation_lags,tukey,butter,sosfilt
from scipy.fft import fft,ifft,fftshift,ifftshift,fftfreq,rfftfreq,rfft
import matplotlib.pyplot as plt
import os


path='/home/doctor/Doctor/Magister/Tesis/databases/Los_Presidentes_DA/'
#data=np.loadtxt(path+'e1e2.txt')
data=np.loadtxt(path+'SynchronizedZ.dat')
fs=256
## get only 25 min 30s
#end=fs*60*25+fs*30
## cantidad de minutos que se graba
mins=np.floor(len(data)/fs/60)
#T=int(len(data)*fs**(-1))
e1=data[:int(mins*fs*60),0]
e2=data[:int(mins*fs*60),1]
#e3=data[:int(mins*fs*60),2]
#e4=data[:int(mins*fs*60),3]
#e5=data[:int(mins*fs*60),4]
#plt.plot(e1)
#plt.plot(e4)
size=int(fs*60)
step=int(fs*60/2)
A=[e1[i:i+size] for i in range(0,len(e1),step)]
B=[e2[i:i+size] for i in range(0,len(e2),step)]
A=A[:-1]
B=B[:-1]
#W=hann(fs*60)
W=tukey(256*60,alpha=0.25)
#arr=np.zeros((50,int(fs*60)))
arr=np.empty([len(A),int(fs*60)]).astype(complex)
arrx=np.empty([len(A),int(fs*60)])
Csuc=[]
Rsuc=[]
c=0
for x,y in zip(A,B):
    xx = W * x
    yy = W * y
    xxx=detrend(xx,overwrite_data=True)
    yyy=detrend(yy,overwrite_data=True)
    #f,Sxx=periodogram(x,fs=fs,window=W,return_onesided=False)
    #f,Syy=periodogram(x,fs=fs,window=W,return_onesided=False)
    #Cxy=Sxx*np.conj(Syy)/(np.abs(Sxx)*np.abs(Syy))

    fftx = fft(xxx)[1:]
    ffty = fft(yyy)[1:]
    #Rxy=correlate(x,y,mode='same')
    Cxy=fftx*np.conj(ffty)/(np.abs(fftx)*np.abs(ffty))
    Cxy=np.hstack(([0.],Cxy))
    #ang=np.angle(Cxy)
    #fx,Cxy=coherence(xx,yy,fs=fs,window=W,noverlap=0,detrend=False)
    #print(Cxy.shape,fx.shape)
    #arr[c,:]=Cxy/(np.max(abs(Cxy)))
    arr[c,:]=Cxy
    #arrx[c,:]=Rxy
    Csuc.append(np.mean(arr[:c+1,:],axis=0))
    #Rsuc.append(np.mean(arrx[:c+1,:],axis=0))

    c+=1

t=np.linspace(-25.5,25.5,len(xx))
Cxyf=np.mean(arr,axis=0)
Cxyt=np.real(ifft(Cxyf))
Cxytsuc=[np.real(ifft(x)) for x in Csuc]
lags=correlation_lags(xx.size,yy.size,mode='same')
n=360
soslp=butter(10,0.10,btype='lowpass',output='sos',fs=fs)
soshp=butter(10,80,btype='highpass',output='sos',fs=fs)
sosbp=butter(10,(0.1,15),btype='bandpass',output='sos',fs=fs)

filt_lp=sosfilt(soslp,Cxytsuc[-1])
filt_hp=sosfilt(soshp,Cxytsuc[-1])
filt_bp=sosfilt(sosbp,Cxytsuc[-1])

fig,ax=plt.subplots(4,1,figsize=(10,8))
ns=[0,50,150,360]
ax[0].plot(lags/(256),fftshift(Cxytsuc[ns[0]]/(np.max(Cxytsuc[ns[0]]))),label='1 window')
ax[0].plot(lags/(256),fftshift(Cxytsuc[ns[1]]/(np.max(Cxytsuc[ns[1]]))),label='50 windows')
ax[0].plot(lags/(256),fftshift(Cxytsuc[ns[2]]/(np.max(Cxytsuc[ns[2]]))),label='150 windows')
ax[0].plot(lags/(256),fftshift(Cxytsuc[ns[3]]/(np.max(Cxytsuc[ns[3]]))),label='360 windows')
ax[1].plot(lags/(256),fftshift(filt_lp)/(np.max(filt_lp)),color='k',label='lowpass 0.1 Hz')
ax[2].plot(lags/(256),fftshift(filt_hp)/(np.max(filt_hp)),color='r',label='highpass 5 Hz')
ax[3].plot(lags/(256),fftshift(filt_bp)/(np.max(filt_bp)),color='m',label='bandpass 0.1-15 Hz')

#plt.plot(lags/(256*60),Rxyt/np.max(Rxyt),label='correl',zorder=1)
plt.xlabel('lag time [seconds]')
ax[0].legend()
ax[1].legend()
ax[2].legend()
ax[3].legend()
plt.savefig('result1.png',format='png',dpi=300,bbox_inches='tight')

#fsr=rfftfreq(len(Cxytsuc[n]),d=1/fs)

#Rxyt=np.mean(arrx,axis=0)


# coh2=coherence(e2,e3,fs=512,nperseg=len(x),noverlap=step,detrend='linear')
# plt.plot(t,fftshift(Cxyt)/np.max(Cxyt),label='coherence')
# plt.plot(t,Rxyt/np.max(Rxyt),label='correl')
# plt.xlabel('lag time [min]')
# plt.legend()
# plt.savefig('ex_.png',format='png',dpi='figure')

#Cxy=Cxy/np.max(Cxy)
#arcor=correlate(Cxy,Cxy)
#e1sp=np.split(e1,25)
#e2sp=np.split(e2,25)

#t=np.linspace(0,25,len(e1))
#plt.plot(t,e1)
#plt.plot(t,e1)