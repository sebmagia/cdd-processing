import os
import numpy as np
from datetime import datetime,time
path='/home/doctor/Doctor/Magister/Tesis/databases/LP_DENSEARRAY/'

path1=path+'M10A/M10A_002part1/EqualizedFile.dat'
print(os.listdir(path))
## idea: hacer una clase cuyos atributos sean la trace length y todas esas weas
metadata=[]
with open(path1,'r',encoding='ISO-8859-1') as file:
    c=0
    for line in file:
        if line[0:6]!='[mm/s]':
            tup=line.partition(':')

            metadata.append([x.strip() for x in tup])

            #print(line,c)
            c+=1
        else:
            print(line,c)
            tup=line.partition(':')
            metadata.append([x.strip() for x in tup])
            c+=1
            break
dic={}
dic.update({'Site_ID':metadata[0][2].split(',')[0]})
dic.update({'Med_and_Node':metadata[0][2].split(',')[1].strip()[:-3]})
dic.update({'Node_Instrument':metadata[0][2].split(',')[1].split('_')[1]+'_'+metadata[0][2].split(',')[1].strip()[-2:]})
dic.update({'Instrument':metadata[0][2].split(',')[1].strip()[-2:]})
dic.update({'Instrument_Tag':metadata[2][2]})
dic.update({'Data_Format_Bytes':int(metadata[4][2][:2])})
dic.update({'Full_scale_mV':int(metadata[5][2][:3])})
dic.update({'N_Channels':int(metadata[6][2])})
dic.update({'Sampling_Rate':int(metadata[7][2][:3])})
dic.update({'Start_Time':datetime.strptime(metadata[9][2],"%d/%m/%y %H:%M:%S")})
dic.update({'End_Time':datetime.strptime(metadata[10][2],"%d/%m/%y %H:%M:%S")})
dic.update({'Trace_length':metadata[11][2]})
dic.update({'Start_Recording_UTC':metadata[18][2].split('\t')[1]}) ## will check this later
dic.update({'Latitude':metadata[20][2]})
dic.update({'Longitude':metadata[21][2]})
dic.update({'Horizontal_Diluition':metadata[23][2]})
dic.update({'Geoid_Altitude':metadata[24][2]})

print(dic)
data=np.loadtxt(path1,skiprows=c,encoding='ISO-8859-1')

