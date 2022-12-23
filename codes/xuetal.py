import numpy as np
from scipy.signal import detrend,hann,coherence,correlate,periodogram,correlation_lags
from scipy.fft import fft,ifft,fftshift,ifftshift
import matplotlib.pyplot as plt
import os


path='/home/doctor/Doctor/Magister/Tesis/databases/Los_Presidentes/Rect_1_27/'
#data=np.loadtxt(path+'e1e2.txt')
data=np.loadtxt(path+'SynchronizedZ.dat')
fs=512
## get only 25 min 30s
end=fs*60*25+fs*30
#T=int(len(data)*fs**(-1))
e1=data[:end,0]
e2=data[:end,0]
e3=data[:end,2]
e4=data[:end,3]
e5=data[:end,4]

size=int(fs*60)
step=int(fs*60/2)
A=[e2[i:i+size] for i in range(0,len(e2),step)]
B=[e3[i:i+size] for i in range(0,len(e3),step)]
A=A[:-1]
B=B[:-1]
W=hann(fs*60)
#arr=np.zeros((50,int(fs*60)))
arr=np.empty([50,int(fs*60)]).astype(complex)
arrx=np.empty([50,int(fs*60)])
Csuc=[]
Rsuc=[]
c=0
for x,y in zip(A,B):
    detrend(x,overwrite_data=True)
    detrend(y,overwrite_data=True)
    #f,Sxx=periodogram(x,fs=fs,window=W,return_onesided=False)
    #f,Syy=periodogram(x,fs=fs,window=W,return_onesided=False)
    #Cxy=Sxx*np.conj(Syy)/(np.abs(Sxx)*np.abs(Syy))
    xx=W*x
    yy=W*y
    fftx = fftshift(fft(xx, norm='ortho'))
    ffty = fftshift(fft(yy, norm='ortho'))
    Rxy=correlate(xx,yy,mode='same')
    Cxy=fftx*np.conj(ffty)/(np.abs(fftx)*np.abs(ffty))
    ang=np.angle(Cxy)
    #fx,Cxy=coherence(xx,yy,fs=fs,window=W,noverlap=0,detrend=False)
    #print(Cxy.shape,fx.shape)
    arr[c,:]=Cxy
    arrx[c,:]=Rxy
    Csuc.append(np.mean(arr[:c+1,:],axis=0))
    Rsuc.append(np.mean(arrx[:c+1,:],axis=0))

    c+=1

Cxyf=np.mean(arr,axis=0)
Cxyt=np.real(ifft(Cxyf))
Cxytsuc=[np.real(ifft(x)) for x in Csuc]
Rxyt=np.mean(arrx,axis=0)
lags=correlation_lags(xx.size,yy.size,mode='full')
coh2=coherence(e2,e3,fs=512,nperseg=len(x),noverlap=step,detrend='linear')
t=np.linspace(-25.5,25.5,len(xx))
plt.plot(t,fftshift(Cxyt)/np.max(Cxyt),label='coherence')
plt.plot(t,Rxyt/np.max(Rxyt),label='correl')
plt.xlabel('lag time [min]')
plt.legend()
plt.savefig('ex_.png',format='png',dpi='figure')

#Cxy=Cxy/np.max(Cxy)
#arcor=correlate(Cxy,Cxy)
#e1sp=np.split(e1,25)
#e2sp=np.split(e2,25)

#t=np.linspace(0,25,len(e1))
#plt.plot(t,e1)
#plt.plot(t,e1)