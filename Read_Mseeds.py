import datetime
import itertools
import obspy
import funciones_linda as fl
import os
from datetime import datetime
import multiprocessing
from scipy.special import jv,jn_zeros
from scipy import interpolate
import matplotlib.pyplot as plt
from pandas import DataFrame,read_csv
import numpy as np
from scipy.fft import rfft, fft, ifft, fftshift, fftfreq,rfftfreq
from scipy.signal import detrend, correlate, correlation_lags, tukey, butter, sosfilt, \
    welch,csd,resample,coherence
import datetime
from functools import partial
from itertools import product
import time
from datetime import datetime
import scipy.sparse as scsp
from scipy.linalg import block_diag
from scipy.sparse.linalg import bicg,LinearOperator,aslinearoperator
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
Dic21=True
Hs21=False
Mar22=False
Sept22=False


path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'

if Dic21:
    path+='Dic21/'
    nw='SA'
if Hs21:
    path+='Hs21/'
if Mar22:
    path+='Mar22/'
if Sept22:
    path+='Sept22/'
files=os.listdir(path)

mseeds=sorted([x for x in files if '.mseed' in x])
coords=sorted([x for x in files if 'coords' in x])

stream_o=obspy.read(path+mseeds[0],format='MSEED')
coord=np.loadtxt(path+coords[0])
equipos=['E1','E2','E3','E4','E5']
nodos=['N1','N2','N3','N4','N5']
nodos=np.loadtxt(path+'nodos.txt',usecols=(0,1,3,4,9,10,11,12,13,14,15)).T
inv = Inventory(networks=[], source='Seba_LP')
net = Network(code=nw,
              stations=[],
              description="Dic 2021 Array.",
              start_date=obspy.UTCDateTime(2016, 1, 2))
stream=stream_o.copy()


for n,x in enumerate(stream):
    sta=Station(code=equipos[n],latitude=coord[n,0],longitude=coord[n,1],elevation=0.0,site=Site(name=nodos[n]),creation_date=obspy.UTCDateTime(2016, 1, 2))
    cha=Channel(code='HHZ',location_code="",latitude=coord[n,0],longitude=coord[n,1],elevation=0.0,depth=0.0,sample_rate=512)
    sta.channels.append(cha)
    net.stations.append(sta)

inv.networks.append(net)
inv.write("new_station.xml", format="stationxml", validate=True)
# for z in stream:
#     z.stats.network = nw
#     sta = Station(code=)

    # for x,y in zip(mseeds,coords):
#     stream = obspy.read(path + x, format='MSEED')
#     coord=np.loadtxt(path+y)


    #stream.write(path+x.split('.')[0]+'_N.'+x.split('.')[1],format='MSEED')


    #stream.write(path+x.split('.')[0]+'_N.'+x.split('.')[1],format='MSEED')




# if Dic21:s
#     for x in mseeds:
#         obspy.read()
# #coords=