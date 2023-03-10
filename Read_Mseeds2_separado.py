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
Sept22=False
Dic22=True


path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'

if Dic21:
    path+='Dic21/'
    files = os.listdir(path)
    mseeds = sorted([x for x in files if '.mseed' in x])
    mseeds=mseeds[2:]+mseeds[:2] ## dic21
    coords = sorted([x for x in files if 'coords' in x])
    coords = coords[2:] + coords[:2]  ## dic21
    equipos = ['E1', 'E2', 'E3', 'E4', 'E5']
    nw='SA'
    meds = ['M1B', 'M2A', 'M3B', 'M4', 'M6B', 'M6C', 'M7A', 'M8A', 'M9A', 'M10A', 'M11A']  ## dic 21

if Hs21:
    path+='HS21/'
    files = os.listdir(path)
    nw='SB'
    equipos = ['E1', 'E2', 'E3', 'E4', 'E5']
    meds=['Lin_1_20','Lin_2_20','Rect_1_27'] ## hs21
    mseeds = sorted([x for x in files if '.mseed' in x])
    coords = sorted([x for x in files if 'coords' in x])

if Mar22:
    path+='Mar22/'
    files = os.listdir(path)
    nw='SC'
    equipos = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6']
    meds=['M1P','M2P','M3P','M4P','M5P','M6P','M7P'] ## mar22
    mseeds = sorted([x for x in files if '.mseed' in x])
    coords = sorted([x for x in files if 'coords' in x])


if Sept22:
    path+='Sept22/'
    files = os.listdir(path)
    nw='SD'
    equipos = ['E1', 'E2', 'E3', 'E4', 'E5']
    meds=['ARR1A','ARR1B','ARR1C','ARR2A','ARR2B','ARR2C','ARR3A','ARR3B','ARR3C','ARR3D']
    mseeds = sorted([x for x in files if '.mseed' in x])
    coords = sorted([x for x in files if 'coords' in x])

if Dic22:
    path+='Dic22/'
    files = os.listdir(path)
    nw='SE'
    equipos = ['E1', 'E2', 'E3', 'E4', 'E5']
    meds=['MED1A','MED1B'] ## dic22
    mseeds = sorted([x for x in files if '.mseed' in x])
    coords = sorted([x for x in files if 'coords' in x])

fs=512

coords=[x for x in coords if 'all' not in x]
stream_o=obspy.read(path+mseeds[0],format='MSEED')
coord=np.loadtxt(path+coords[0])



nodos=np.loadtxt(path+'nodos.txt').T
inv = Inventory(networks=[], source='Seba_LP')
net = Network(code=nw,
              stations=[],
              description="Dict 22 arrays.",
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
        number=int(k.stats.station[1])-1 ## for assigning correct coordinates!
        sta=Station(code=equipos[number],latitude=coord[number,1],longitude=coord[number,0],elevation=0.0,site=Site(name=nods[n]),creation_date=obspy.UTCDateTime(2016, 1, 2))
        cha=Channel(code='HHZ',location_code="",latitude=coord[number,1],longitude=coord[number,0],elevation=0.0,depth=0.0,sample_rate=512)
        sta.channels.append(cha)
        net.stations.append(sta)
    inv.networks.append(net)

    c+=1
combs_flat=fl.flatten(combs)
meds_flat=fl.flatten(meds_l)
dicc=fl.repeticiones(combs_flat,meds_flat)
net_names = [x.description for x in inv.networks]
tr1=[]
tr2=[]
xtream=obspy.Stream()
k=0
stack1 = None
stack2 = None
keystack=[]
valuestack=[]
#coords=np.zeros((len(set(combs_flat)),4))
coords=np.zeros((len(combs_flat),4))
valores=[]
for key,values in zip(dicc.keys(),dicc.values()):
    print(key,values)
    keystack.append(key)
    valuestack.append(values)
    if len(values)>1:
        #del stack1
        #del stack2
        for x in values:
            valores.append(x)
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
            tr1 = obspy.Trace(data=trace1, header={'sampling_rate': fs, 'station': 'A'})
            tr2 = obspy.Trace(data=trace2, header={'sampling_rate': fs, 'station': 'B'})
            st = obspy.Stream(traces=[tr1, tr2])
            # xtream+=st
            st.write(path + 'STACKS_2/par' + str(k).zfill(3) + '.mseed', format='MSEED')
            k += 1
            #print('len',trace1.shape,trace2.shape)

            # try:
            #     stack1 = np.hstack((stack1, trace1))
            #     stack2 = np.hstack((stack2, trace2))
            #     #coordinates= np.vstack((coordinates,arr))
            # except NameError:
            #     print('exc')
            #     stack1 = trace1
            #     stack2 = trace2
                #coordinates = arr
        #print(coordinates)
    else:
        #del stack1
        #del stack2
        valores.append(values[0])
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
        print('indices',idx_b1,idx_b2)
        print('nodes',nodes.index(key[0]))
        #print('bar', inv[idx_a].stations[idx_b1], inv[idx_a].stations[idx_b2])
        trace1 = obspy.read(path + values[0] + '.mseed', format='MSEED')[idx_b1].data
        trace2 = obspy.read(path + values[0] + '.mseed', format='MSEED')[idx_b2].data


        tr1 = obspy.Trace(data=trace1, header={'sampling_rate': fs, 'station': 'A'})
        tr2 = obspy.Trace(data=trace2, header={'sampling_rate': fs, 'station': 'B'})
        st = obspy.Stream(traces=[tr1, tr2])
        st.write(path+'STACKS_2/par'+str(k).zfill(3)+'.mseed',format='MSEED')

        #xtream += st
        k+=1
np.savetxt(path+'coords_all2.txt',coords)
with open(path+'keys_nodes.txt','w') as fileobject:
    for x in keystack:
        fileobject.write(x[0]+' '+x[1])
        fileobject.write('\n')

with open(path+'values_nodes.txt','w') as fileobject:
    for x in valuestack:
        for y in x:
            fileobject.write(y+' ')
        fileobject.write('\n')
with open(path+'values_separado.txt','w') as fileobject:
    for x in valores:
        fileobject.write(x)
        fileobject.write('\n')

coords1=np.split(coords,2,axis=1)[0]
coords2=np.split(coords,2,axis=1)[1]

for x,y in zip(coords1,coords2):
    xd=np.vstack((x,y))
    print(xd)
    plt.plot(xd[:,0],xd[:,1],'ko-')
inv.write(path+"HS21.xml", format="stationxml", validate=True)

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