import obspy
import hvsrpy
import numpy as np
import os
path='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'
path_output='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/HV2/'

folders=os.listdir(path)
folders=[path+x for x in folders]
folders=filter(os.path.isdir,folders)


folders=[x for x in folders if 'MSEEDS' not in x and 'SUR1B' not in x and 'SUR2B' not in x]
#folders=[x for x in folders if '.txt' not in x]
#folders=[x for x in folders if 'SUR1B' not in x]
#folders=[x for x in folders if 'SUR2B' not in x]

folders=sorted(folders)
print(folders)
k=0
for folder in folders:
    if 'SUR1AX' in folder:
        k = 0
        subfolders=os.listdir(folder)
        subfolders=[folder + '/' + x for x in subfolders]
        subfolders = filter(os.path.isdir, subfolders)
        subfolders = sorted([ x for x in subfolders if 'COPYS' not in x and 'REJECTED' not in x])
        subfolders12=subfolders[:2]
        subfolders=subfolders[2:]
        subfolders.append(subfolders12[0])
        subfolders.append(subfolders12[1])

        nodos_file='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Dic21/nodos.txt'
        nodos=np.loadtxt(nodos_file).astype(int)-1
        coords_file='/home/doctor/Doctor/Magister/Tesis/databases/process_data/coords_SUR1A.txt'
        coords=np.loadtxt(coords_file)
        for i,subf in enumerate(subfolders):
            subsubfolders = os.listdir(subf)
            subsubfolders = [subf + '/' + x for x in subsubfolders]
            subsubfolders = sorted(filter(os.path.isdir, subsubfolders))
            for n,subsubf in enumerate(subsubfolders):
                nodos_i=nodos[:,i]
                nodos_i=[x for  x in nodos_i if x != -1 ]
                sta_xy=coords[nodos_i[n]]
                arr=np.genfromtxt(subsubf+'/EqualizedFile.dat', skip_header=34,encoding='ISO-8859-1',usecols=(0,1,2))
                trN=obspy.Trace(arr[:,0], header={'sampling_rate':512,'channel':'HHN'})
                trE=obspy.Trace(arr[:,1], header={'sampling_rate':512,'channel':'HHE'})
                trZ=obspy.Trace(arr[:,2], header={'sampling_rate':512,'channel':'HHZ'})
                stream=obspy.Stream(traces=[trN,trE,trZ])
                stream.write(path_output+'Dic21/data_3C_'+str(k).zfill(3)+'.mseed')
                with open(path_output+'Dic21/single_coords.txt','a') as fileobject:
                    fileobject.write(str(k).zfill(3)+' '+str(np.round(sta_xy[0],2))+' '+str(np.round(sta_xy[1],2)))
                    fileobject.write('\n')
                k+=1
                print(subsubf)
                print(sta_xy)
                print(arr[0])

    if 'SUR2AX' in folder:
        k = 0
        subfolders = os.listdir(folder)
        subfolders = [folder + '/' + x for x in subfolders]
        subfolders = filter(os.path.isdir, subfolders)
        subfolders = sorted([x for x in subfolders if 'Seba_Udec' not in x ])


        nodos_file = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Mar22/nodos.txt'
        nodos = np.loadtxt(nodos_file).astype(int) - 1
        coords_file = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/coords_SUR2A.txt'
        coords = np.loadtxt(coords_file)
        for i, subf in enumerate(subfolders):
            subsubfolders = os.listdir(subf)
            subsubfolders = [subf + '/' + x for x in subsubfolders]
            subsubfolders = sorted(filter(os.path.isdir, subsubfolders))
            for n, subsubf in enumerate(subsubfolders):
                nodos_i = nodos[:, i]
                nodos_i = [x for x in nodos_i if x != -1]
                sta_xy = coords[nodos_i[n]]
                arr = np.genfromtxt(subsubf + '/EqualizedFile.dat', skip_header=34, encoding='ISO-8859-1',usecols=(0, 1, 2))
                trN = obspy.Trace(arr[:, 0], header={'sampling_rate': 512, 'channel': 'HHN'})
                trE = obspy.Trace(arr[:, 1], header={'sampling_rate': 512, 'channel': 'HHE'})
                trZ = obspy.Trace(arr[:, 2], header={'sampling_rate': 512, 'channel': 'HHZ'})
                stream = obspy.Stream(traces=[trN, trE, trZ])
                stream.write(path_output + 'Mar22/data_3C_' + str(k).zfill(3) + '.mseed')
                with open(path_output + 'Mar22/single_coords.txt', 'a') as fileobject:
                     fileobject.write(str(k).zfill(3) + ' ' + str(np.round(sta_xy[0], 2)) + ' ' + str(np.round(sta_xy[1], 2)))
                     fileobject.write('\n')
                k += 1
                print(subsubf)
                print(sta_xy)
                print(arr[0])
    if 'Los_Presidentes_HSX' in folder:
        k = 0
        subfolders = os.listdir(folder)
        subfolders = [folder + '/' + x for x in subfolders]
        subfolders = filter(os.path.isdir, subfolders)
        subfolders = sorted([x for x in subfolders if 'Act' not in x ])
        subfolders = sorted([x for x in subfolders if 'Lin_3_20' not in x and 'Lin_4_20' not in x and 'Lin_5_20' not in x ])
        subfolders = sorted([x for x in subfolders if 'Rect_2_26' not in x and 'Rect_3_23' not in x])

        nodos_file = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/HS21/nodos.txt'
        nodos = np.loadtxt(nodos_file).astype(int) - 1

        coords_file = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/coords_HS21.txt'
        coords = np.loadtxt(coords_file)
        for i, subf in enumerate(subfolders):
            subsubfolders = os.listdir(subf)
            subsubfolders = [subf + '/' + x for x in subsubfolders]
            subsubfolders = sorted(filter(os.path.isdir, subsubfolders))
            nodos_i=nodos[:,i]
            for n, subsubf in enumerate(subsubfolders):
                sta_xy = coords[nodos_i[n]]
                arr = np.genfromtxt(subsubf + '/EqualizedFile.dat', skip_header=44, encoding='ISO-8859-1',usecols=(0, 1, 2))
                trN = obspy.Trace(arr[:, 0], header={'sampling_rate': 512, 'channel': 'HHN'})
                trE = obspy.Trace(arr[:, 1], header={'sampling_rate': 512, 'channel': 'HHE'})
                trZ = obspy.Trace(arr[:, 2], header={'sampling_rate': 512, 'channel': 'HHZ'})
                stream = obspy.Stream(traces=[trN, trE, trZ])
                stream.write(path_output + 'HS21/data_3C_' + str(k).zfill(3) + '.mseed')
                with open(path_output + 'HS21/single_coords.txt', 'a') as fileobject:
                     fileobject.write(str(k).zfill(3) + ' ' + str(np.round(sta_xy[0], 2)) + ' ' + str(np.round(sta_xy[1], 2)))
                     fileobject.write('\n')
                k += 1
                print(subsubf)
                print(sta_xy)
                print(arr[0])
    if 'ARRX' in folder:
        subfolders = os.listdir(folder)
        subfolders = [folder + '/' + x for x in subfolders]
        subfolders = sorted(filter(os.path.isdir, subfolders))        #break
        coords=np.loadtxt(folder+'/coords_'+folder.split('/')[-1]+'.txt')
        for n,subf in enumerate(subfolders):
            sta_xy = coords[n]
            arr = np.genfromtxt(subf + '/EqualizedFile.dat', skip_header=44, encoding='ISO-8859-1',usecols=(0, 1, 2))
            trN = obspy.Trace(arr[:, 0], header={'sampling_rate': 512, 'channel': 'HHN'})
            trE = obspy.Trace(arr[:, 1], header={'sampling_rate': 512, 'channel': 'HHE'})
            trZ = obspy.Trace(arr[:, 2], header={'sampling_rate': 512, 'channel': 'HHZ'})
            stream = obspy.Stream(traces=[trN, trE, trZ])
            stream.write(path_output + 'Sept22/data_3C_' + str(k).zfill(3) + '.mseed')
            with open(path_output + 'Sept22/single_coords.txt', 'a') as fileobject:
                fileobject.write(
                    str(k).zfill(3) + ' ' + str(np.round(sta_xy[0], 2)) + ' ' + str(np.round(sta_xy[1], 2)))
                fileobject.write('\n')
            k += 1
            print(subf)
            print(sta_xy)
            print(arr[0])
    if 'LP_DIC22' in folder:
        subfolders = os.listdir(folder)
        subfolders = [folder + '/' + x for x in subfolders]
        subfolders = sorted(filter(os.path.isdir, subfolders))
        subfolders = [x for x in subfolders if 'Act' not in x ]
        coordsA=np.loadtxt(subfolders[0]+'/coords_'+subfolders[0].split('/')[-1]+'.txt')
        coordsB=np.loadtxt(subfolders[1]+'/coords_'+subfolders[1].split('/')[-1]+'.txt')
        for subf in subfolders:
