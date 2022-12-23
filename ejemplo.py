

import obspy
import numpy as np


javier=False
path='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'

#med='M'
k=9
#dataN=np.loadtxt(path+'SUR1A/'+med+'/SynchronizedNS.dat')
#dataN=np.loadtxt(path+'Los_Presidentes_HS/'+med+'/SynchronizedNS.dat')
#dataN=np.loadtxt(path+'SUR2A/'+med+'/SynchronizedNS.dat')
dataN=np.loadtxt(path+'ARR3D/SynchronizedNS.dat')

#dataE=np.loadtxt(path+'SUR1A/'+med+'/SynchronizedEW.dat')
#dataE=np.loadtxt(path+'Los_Presidentes_HS/'+med+'/SynchronizedEW.dat')
#dataE=np.loadtxt(path+'SUR2A/'+med+'/SynchronizedEW.dat')
dataE=np.loadtxt(path+'ARR3D/SynchronizedEW.dat')

#dataZ=np.loadtxt(path+'SUR1A/'+med+'/SynchronizedZ.dat')
#dataZ=np.loadtxt(path+'Los_Presidentes_HS/'+med+'/SynchronizedZ.dat')
#dataZ=np.loadtxt(path+'SUR2A/'+med+'/SynchronizedZ.dat')
dataZ=np.loadtxt(path+'ARR3D/SynchronizedZ.dat')

#orden_dat=np.loadtxt(path+'SUR1A/orden_instrumentos.txt').T
#orden_dat=np.loadtxt(path+'Los_Presidentes_HS/orden_instrumentos.txt').T
#orden_dat=np.loadtxt(path+'SUR2A/orden_instrumentos.txt').T
orden_dat=np.loadtxt(path+'orden_instrumentos_sept22.txt').T

nw='SD'
starttime1="2022/01/01 17:27:20"
endtime1= "2022/01/11 18:11:20"

fs=512
nodos=np.loadtxt(path+'MSEEDS/Dic21/nodos.txt').T
nodos=np.loadtxt(path+'MSEEDS/HS21/nodos.txt').T
nodos=np.loadtxt(path+'MSEEDS/Mar22/nodos.txt').T
nodos=np.loadtxt(path+'MSEEDS/Sept22/nodos.txt').T
stream=obspy.Stream()
for i in range(len(dataN.T)):
    nod=[str(int(x)) for x in nodos[k] if int(x) != 0]
    inst=[int(x) for x in orden_dat[k] if int(x) != 0]
    traceE= obspy.Trace(data=dataE[:,inst[i]-1], header={'sampling_rate': fs, 'network': nw, 'station': 'N'+nod[i],'channel':'HHE'})
    traceN= obspy.Trace(data=dataN[:,inst[i]-1], header={'sampling_rate': fs, 'network': nw, 'station': 'N'+nod[i],'channel':'HHN'})
    traceZ= obspy.Trace(data=dataZ[:,inst[i]-1], header={'sampling_rate': fs, 'network': nw, 'station': 'N'+nod[i],'channel':'HHZ'})
    stream.append(traceE)
    stream.append(traceN)
    stream.append(traceZ)
    stream.write(path+'MSEEDS/Sept22/ARR3D_all.mseed',format='MSEED')

if javier:
    jav = np.loadtxt(path + 'SUR2A/' + med + '/Seba4_001part1/EqualizedFile.dat', skiprows=42,
                         encoding='ISO-8859-1', usecols=(0, 1, 2))
    javN = jav[:, 0]
    javE = jav[:, 1]
    javZ = jav[:, 2]
    traceE= obspy.Trace(data=javE, header={'sampling_rate': fs, 'network': nw, 'station': 'N'+str(int(nodos[k][-1])),'channel':'HHE'})
    traceN= obspy.Trace(data=javN, header={'sampling_rate': fs, 'network': nw, 'station': 'N'+str(int(nodos[k][-1])),'channel':'HHN'})
    traceZ= obspy.Trace(data=javZ, header={'sampling_rate': fs, 'network': nw, 'station': 'N'+str(int(nodos[k][-1])),'channel':'HHZ'})
    stream.append(traceE)
    stream.append(traceN)
    stream.append(traceZ)
    stream.write(path+'MSEEDS/Sept22/ARR1A_all.mseed',format='MSEED')

        #'starttime':obspy.UTCDateTime(starttime1),'endtime':obspy.UTCDateTime(endtime1)})
#stream=obspy.read(path+'M1B.mseed')