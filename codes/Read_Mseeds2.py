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
from collections import Counter
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
Dic21=False
Hs21=False
Mar22=False
Sept22=True


path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'

if Dic21:
    path+='Dic21/'
    nw='SA'
if Hs21:
    path+='HS21/'
    nw='SB'
if Mar22:
    path+='Mar22/'
    nw='SC'
if Sept22:
    path+='Sept22/'
    nw='SD'
files=os.listdir(path)
fs=512
mseeds=sorted([x for x in files if '.mseed' in x])
#mseeds=mseeds[2:]+mseeds[:2]
coords=sorted([x for x in files if 'coords' in x])
#coords=coords[2:]+coords[:2]

stream_o=obspy.read(path+mseeds[0],format='MSEED')
coord=np.loadtxt(path+coords[0])
#equipos=['E1','E2','E3','E4','E5'] ## dic 21
equipos=['E1','E2','E3','E4','E5','E6'] # mar 22
equipos=['E1','E2','E3','E4','E5'] # HS 22
equipos=['E1','E2','E3','E4','E5'] # Sept 22

#nodos=['N1','N2','N3','N4','N5']
#meds=['M1B','M2A','M3B','M4','M6B','M6C','M7A','M8A','M9A','M10A','M11A'] ## dic 21
#meds=['M1P','M2P','M3P','M4P','M5P','M6P','M7P']
#meds=['Lin_1_20','Lin_2_20','Rect_1_27']
meds=['ARR1A','ARR1B','ARR1C','ARR2A','ARR2B','ARR2C','ARR3A','ARR3B','ARR3C','ARR3D'] # p1

#meds=['ARR2A','ARR2B','ARR2C'] # p2
#meds=['ARR3A','ARR3B','ARR3C','ARR3D'] #p3
#nodos=np.loadtxt(path+'nodos.txt',usecols=(0,1,3,4,9,10,11,12,13,14,15)).T ## dic 21
nodos=np.loadtxt(path+'nodos.txt').T ## mar 22, HS21, Sept 22
inv = Inventory(networks=[], source='Seba_LP')
net = Network(code=nw,
              stations=[],
              description="Sept 22 arrays.",
              start_date=obspy.UTCDateTime(2016, 1, 2))
stream=stream_o.copy()
c=0
combs=[]
meds_l=[]
for x,y in zip(mseeds,coords):
    stream = obspy.read(path + x, format='MSEED')

    coord=np.loadtxt(path+y)
    net = Network(code=nw,
                  stations=[],
                  description=meds[c],
                  start_date=obspy.UTCDateTime(2016, 1, 2))
    try:
        nods = ['N' + str(int(i)) for i in nodos[c] if i != 0]
    except ValueError:
        nods=['N' + i for i in nodos[c] if i != 0]
    print(nods)
    comb = list(itertools.combinations(nods, 2))
    combs.append(comb)
    meds_l.append(len(comb)*[meds[c]])
    for n,k in enumerate(stream):
        #comb = list(itertools.combinations(nods, 2))
        #combs.append(comb)
        print(coord[n,0],coord[n,1])
        sta=Station(code=equipos[n],latitude=coord[n,1],longitude=coord[n,0],elevation=0.0,site=Site(name=nods[n]),creation_date=obspy.UTCDateTime(2016, 1, 2))
        cha=Channel(code='HHZ',location_code="",latitude=coord[n,1],longitude=coord[n,0],elevation=0.0,depth=0.0,sample_rate=512)
        sta.channels.append(cha)
        net.stations.append(sta)
    inv.networks.append(net)

    c+=1
combs_flat=fl.flatten(combs)
meds_flat=fl.flatten(meds_l)
dick=fl.repeticiones(combs_flat,meds_flat)
#dick={('N1', 'N2'):['M1B','M2A']}
net_names = [x.description for x in inv.networks]
tr1=[]
tr2=[]
xtream=obspy.Stream()
k=0
stack1 = None
stack2 = None
keystack=[]
coords=np.zeros((len(set(combs_flat)),4))
for key,values in zip(dick.keys(),dick.values()):
    print(key,values)
    keystack.append(key)

    if len(values)>1:
        del stack1
        del stack2
        for x in values:
            idx_a=net_names.index(x)
            print('foo')
            #print('foo',inv[idx_a])
            nodes=[y.site.name for y in inv[idx_a].stations]
            idx_b1=nodes.index(key[0])
            idx_b2=nodes.index(key[1])
            lat1=inv[idx_a].stations[idx_b1].latitude
            lon1=inv[idx_a].stations[idx_b1].longitude
            lat2 = inv[idx_a].stations[idx_b2].latitude
            lon2 = inv[idx_a].stations[idx_b2].longitude
            coords[k,:]=np.array([lon1,lat1,lon2,lat2])
            print(coords[k, :])

            #print(lat1,lat2,lon1,lon2)
            #arr=np.array([lon1,lat1,lon2,lon2])
            #print('bar',inv[idx_a].stations[idx_b1],inv[idx_a].stations[idx_b2])
            trace1=obspy.read(path + x+'.mseed', format='MSEED')[idx_b1].data
            trace2=obspy.read(path + x+'.mseed', format='MSEED')[idx_b2].data
            #print('len',trace1.shape,trace2.shape)

            try:
                stack1 = np.hstack((stack1, trace1))
                stack2 = np.hstack((stack2, trace2))
                #coordinates= np.vstack((coordinates,arr))
            except NameError:
                print('exc')
                stack1 = trace1
                stack2 = trace2
                #coordinates = arr
        #print(coordinates)
        tr1 = obspy.Trace(data=stack1, header={'sampling_rate': fs, 'station': 'A'})
        tr2 = obspy.Trace(data=stack2, header={'sampling_rate': fs, 'station': 'B'})
        st=obspy.Stream(traces=[tr1,tr2])
        #xtream+=st
        st.write(path+'STACKS/par'+str(k).zfill(2)+'.mseed',format='MSEED')
        k+=1
    else:
        del stack1
        del stack2
        idx_a = net_names.index(values[0])
        print('koo')
        #print('koo', inv[idx_a])
        nodes = [y.site.name for y in inv[idx_a].stations]

        idx_b1 = nodes.index(key[0])
        idx_b2 = nodes.index(key[1])
        lat1 = inv[idx_a].stations[idx_b1].latitude
        lon1 = inv[idx_a].stations[idx_b1].longitude
        lat2 = inv[idx_a].stations[idx_b2].latitude
        lon2 = inv[idx_a].stations[idx_b2].longitude
        coords[k, :] = np.array([lon1, lat1, lon2, lat2])
        print(coords[k,:])
        #print('bar', inv[idx_a].stations[idx_b1], inv[idx_a].stations[idx_b2])
        trace1 = obspy.read(path + values[0] + '.mseed', format='MSEED')[idx_b1].data
        trace2 = obspy.read(path + values[0] + '.mseed', format='MSEED')[idx_b2].data

        try:
            stack1 = np.hstack((stack1, trace1))
            stack2 = np.hstack((stack2, trace2))
            #coordinates = np.vstack((coordinates, arr))

        except NameError:
            #print('exc')
            stack1 = trace1
            stack2 = trace2
            #coordinates = arr

        tr1 = obspy.Trace(data=stack1, header={'sampling_rate': fs, 'station': 'A'})
        tr2 = obspy.Trace(data=stack2, header={'sampling_rate': fs, 'station': 'B'})
        st = obspy.Stream(traces=[tr1, tr2])
        st.write(path+'STACKS/par'+str(k).zfill(2)+'.mseed',format='MSEED')

        #xtream += st
        k+=1
np.savetxt(path+'coords_all.txt',coords)
coords1=np.split(coords,2,axis=1)[0]
coords2=np.split(coords,2,axis=1)[1]

for x,y in zip(coords1,coords2):
    xd=np.vstack((x,y))
    print(xd)
    plt.plot(xd[:,0],xd[:,1],'ko-')
inv.write(path+"Dic21.xml", format="stationxml", validate=True)

#print(xtream.__str__(extended=True))
#xtream.write(path+'STACKS/stacked.mseed',format='MSEED')

#read=obspy.read(path+'STACKS/stacked.mseed',format='MSEED')
#print(xtream)
    # else:
    #     idx_a = net_names.index(x)
    #     #print('foo', inv[idx_a])
    #     nodes = [y.site.name for y in inv[idx_a].stations]
    #
    #     idx_b1 = nodes.index(key[0])
    #     idx_b2 = nodes.index(key[1])
    #     #print('bar', inv[idx_a].stations[idx_b1], inv[idx_a].stations[idx_b2])
    #     tr1 = obspy.read(path + x + '.mseed', format='MSEED')[idx_b1].data
    #     tr2 = obspy.read(path + x + '.mseed', format='MSEED')[idx_b2].data
    #     st=obspy.Stream(traces=[tr1,tr2])
    #     xtream+=st


            # try:
            #     stack1=np.vstack((stack1,trace1))
            #     stack2=np.vstack((stack2,trace2))
            #
            # except:
            #     stack1=trace1
            #     stack2=trace2


            #strea
            #trace1=obspy.read()
            #nod=inv[idx].stations.site.name
            #print(idx)
            #if z == net_names[1]:
                #print(x,z)

#counter=Counter(map(combs_flat,meds_flat))
# dicc = {}
# for x, y in zip(combs_flat, meds_flat):
#
#     dicc.setdefault(x, []).append(y)
#
# # for z in stream:
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