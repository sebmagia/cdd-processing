import hvsrpy
import numpy as np
import obspy
import os
path='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/HV/Sept22/'
files=os.listdir(path)
ms=[x for x in files if '.mseed' in x]
#coords=[x for x in files if 'coords' in x]
arrys=[]
#coordenadas=[]
strings=[]
for x in ms:
    med=x.split('_')[0]
    #coords = np.loadtxt(path + 'coords_' + med + '.txt')

    try:
        coords=np.loadtxt(path+'coords_'+med+'.txt')
    except OSError:
        print('HS')
        med=x[:-10]
        coords=np.loadtxt(path+'coords_'+med+'.txt')
    stream=obspy.read(path+x,format='MSEED')
    nsta=int(len(stream)/3)
    for i in range(nsta):
        streamnew=stream[3*i:3*i+3]
        streamnew.write(path+'/SSHV/'+med+'_'+streamnew[0].stats.station+'.mseed',format='MSEED')
        strings.append(med+'_'+streamnew[0].stats.station)
        arry=[x.stats.channel for x in streamnew]
        arrys.append(arry)
        print('koo')
        try:
            coordenadas=np.vstack((coordenadas,coords[i]))
        except:
            coordenadas=coords[i]

with open(path+'SSHV/single_station_coords.txt','w') as fileobject:
    for x,y in zip(strings,coordenadas):
        fileobject.write(x+','+str(np.round(y[0],2))+','+str(np.round(y[1],2)))
        fileobject.write('\n')
        #streams.append(streamnew)
#stream=obspy.read('UT.STN11.A2_C50.miniseed',format='MSEED')

#ver=np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/veri.txt')
#hor=np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/hori.txt')

#horf=hor.flatten()
#verf=ver.flatten()

#div=horf/verf

