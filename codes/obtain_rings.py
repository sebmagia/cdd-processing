import numpy as np
import obspy
from itertools import combinations
from collections import Counter
path='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'

folders = ['Dic21/', 'HS21/','Mar22/', 'Sept22/']
files_dic21=['coords_M1B.txt','coords_M2A.txt','coords_M3B.txt','coords_M4.txt','coords_M6B.txt','coords_M6C.txt','coords_M7A.txt','coords_M8A.txt',\
             'coords_M9A.txt','coords_M10A.txt','coords_M11A.txt']
files_hs21=['coords_Lin_1_20.txt','coords_Lin_2_20.txt','coords_Rect_1_27.txt']
files_mar22=['coords_M1P.txt','coords_M2P.txt','coords_M3P.txt','coords_M4P.txt','coords_M5P.txt','coords_M6P.txt','coords_M7P.txt']
files_sept22=['coords_ARR1A.txt','coords_ARR1B.txt','coords_ARR1C.txt','coords_ARR2A.txt','coords_ARR2B.txt','coords_ARR2C.txt','coords_ARR3A.txt',\
              'coords_ARR3B.txt','coords_ARR3C.txt','coords_ARR3D.txt']

files=[files_dic21,files_hs21,files_mar22,files_sept22]
#path2=path+folders[1]
#nodos=np.loadtxt(path2+'nodos.txt')
#st=obspy.read(path2+'M7P.mseed',format='MSEED')
header='# MinRadius	MaxRadius Red Green	Blue'
for n,fold in enumerate(folders):
    path2 = path + fold
    nodos = np.loadtxt(path2 + 'nodos.txt')
    c = 0
    for file in files[n]:
        if n==0 or n==2:
            med=file.split('_')[1][:3]
        elif n==1:
            med=file[7:].split('.')[0]
        elif n==3:
            med=file.split('_')[1][:-4]
        print(path2+file)
        data=np.loadtxt(path2+file)

        nodos_c=nodos[:,c]
        idxes=np.where(nodos_c>1e-6)
        stas=['SA_E'+str(int(x)+1) for x in idxes[0]]
        data=np.round(data[idxes],1)
        print(data.shape)
        c+=1
        combis=np.array(list(combinations(data,2)))
        distancias = []
        for comb in combis:
            dist=np.linalg.norm(comb[0,:]-comb[1,:])
            distancias.append(np.round(dist,1))
        eles,repeats=np.unique(distancias,return_counts=True)
        eles=np.round(eles,2)
        ring=np.vstack((eles-0.05,eles+0.05,255*np.ones(len(eles)),np.zeros(len(eles)),np.zeros(len(eles)),repeats)).T
        format=['%2.2f','%2.2f','%i','%i','%i','%i']
        np.savetxt(path2+med+'.rings',ring,header=header,fmt=format)
        with open(path2+'coords_'+med+'_new.txt','w') as fileobject:
            fileobject.write('# Station_name        X      Y       Z')
            fileobject.write('\n')
            fileobject.write('# (station_name cannot contain blanks)')
            fileobject.write('\n')
            for i in range(len(stas)):
                fileobject.write(stas[i]+' '+str(data[i,0])+' '+str(data[i,1])+' '+str(0.0))
                fileobject.write('\n')
        #np.savetxt(path2+'coords_'+med+'_new.txt',data,fmt=['%2.2f','%2.2f'])

