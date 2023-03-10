import numpy as np
import matplotlib.pyplot as plt
import shutil
import swprepost
import os
import matplotlib


path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO/'
path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/INV/'

#cdd_names=['cdd00','cdd01','cdd02', 'cdd03', 'cdd04', 'cdd05', 'cdd06', 'cdd08', 'cdd09', 'cdd12', 'cdd13', 'cdd14', 'cdd15', 'cdd16' ]
# cdd_names=[0,1,3,4,5,6,8,9,12,13,14,15,16,17,19,21,22,23,24,30,31,32,35,37,
#            38,40,41,42,43,45,46,47,48,49,50,51,53,54,55,56,57,58,59,60,61,62,63,64,66,67,70,71,72,73,74,75,77,78,79,82,83,84,90,
#            91,93,95,96,98,99,100,101,102,103,104,107,108,109,110,111,112,116,119,121,123,124,125,127,128,133,134,136,137,138,
#            139,140,141,142,143,145,146,147,148,149]
cdd_names=os.listdir(path_cdd)
cdd_names=sorted([x for x in cdd_names if '.dat' in x and 'coordinates' not in x])
os.listdir('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/qinvs_feb15')
#estado=np.loadtxt(path_cdd+'estado_2.txt').astype(int)
#cdd_names=[x for x,y in zip(cdd_names,estado) if y != 0]
#cdd_names=cdd_names[0:13]
#cdd_names_read=['cdd'+str(x).zfill(2) for x in cdd_names]
#cdd_names=[str(x).zfill(2) for x in cdd_names]
#cdd_names=cdd_names[:10]
#idx=np.where(estado !=0)
coords=np.loadtxt(path_cdd+'cdd_coordinates.dat')
#coords=coords[idx]
subfolder='qinvs_feb15'
hv = False
gen_setup=True
# with open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + subfolder + '/script.sh',
#           'w') as fileobject:
#     fileobject.write('#!/bin/bash')
#     fileobject.write('\n')
#     fileobject.close()

nummer=0

for cdd_name in cdd_names:
    #cdd_name='cdd00'
    cdd=np.loadtxt(path_cdd+cdd_name)

    #cdd_name='cdd00'
    num = int(cdd_name.split('_')[1][:3])
    print(num)
    dist = ((coords[num, 0] - coords[num, 2]) ** 2 + (coords[num, 1] - coords[num, 3]) ** 2) ** 0.5
    #lambdar = 0.45
    #idx = np.where( cdd[:, 3] / dist >= lambdar)[0]
    #cdd = cdd[idx, :]
    np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/dispR',cdd)

    if hv:
        path_hv='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/HV2/Dic21/HV/'
        hv1=np.loadtxt(path_hv+'data_3C_000.hv',skiprows=24,delimiter=',')
        hv2=np.loadtxt(path_hv+'data_3C_002.hv',skiprows=24,delimiter=',')

        gen_setup=True

        amax1=np.argmax(hv1[:,1])
        fmin=1.7
        fmax=cdd[0,0]
        f0=0.67
        idx=np.where((hv1[:,0] <= fmax ) & (hv1[:,0] > fmin))
        freq=hv1[idx,0][0]
        amp=hv1[idx,1][0]
        np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/hvDushanbe',np.vstack((freq,amp)).T)



    if gen_setup:
        inv_for=0
        kind_inv=1
        app_fun=0
        stdrel=0
        file_disR='dispR'
        file_disL='dispL'
        file_hv='hvDushanbe'
        file_test='forw.dat_test'
        seeds=[3,5,7,9,11,13,15,17,19]
        seeds=[3,5,7]
        #seeds=[5]
        ga_params=[100,100,0.7,0.01]
        ga_params=[100,100,0.7,0.01]
        wmin=cdd[0,1]/cdd[0,0]
        wmax=cdd[-1,1]/(cdd[-1,0])
        vmin=cdd[0,1]
        vmax=cdd[-1,1]
        lns=[3,4,5,6,7]
        #lns=[6]
        #lrs=[2.5,1.8,1.2]
        depth_factor=2
        fac=0.15 ## factor por minmax velocities
        #depths_lr = swprepost.Parameter.depth_lr(wmin, wmax, lrs[3], depth_factor)
        depths_ln = swprepost.Parameter.depth_ln(wmin, wmax, 4, depth_factor)
        # try:
        #     os.mkdir('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' +  subfolder+ '/' + cdd_name)
        # except:
        #     pass
        # with open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+ subfolder + '/' +cdd_name + '/script.sh', 'w') as fileobject:
        #    fileobject.write('#!/bin/bash')
        #    fileobject.write('\n')
        #    fileobject.close()
        # with open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+ subfolder + '/' +cdd_name + '/script.sh', 'w') as fileobject:
        #    fileobject.write('#!/bin/bash')
        #    fileobject.write('\n')
        #    fileobject.close()
        lrs_trial=np.arange(1,10,0.05)[1:]
        trials3=[]
        trials4=[]
        trials5=[]
        trials6=[]
        layer_trials=[]

        for trial in lrs_trial:
            try:
                depths = swprepost.Parameter.depth_lr(wmin, wmax, trial, depth_factor)
            except IndexError:
                print('aq')
                break
            #print(depths)

            nl = len(depths[0])

            #if nl >= 3 and nl <= 6:
            if nl == 3:
                trials3.append(trial)
            if nl == 4:
                trials4.append(trial)
            if nl == 5:
                trials5.append(trial)
            if nl == 6:
                trials6.append(trial)
        lrs=np.array([np.median(trials3),np.median(trials4),np.median(trials5),np.median(trials6)])
        lrs=lrs[~np.isnan(lrs)]
        lrs=np.round(lrs,2)

        #lrs=[4,5,6,7]
        for seed in seeds:
            for lr in lrs:
                depths = swprepost.Parameter.depth_lr(wmin, wmax, lr, depth_factor)
                #depths = swprepost.Parameter.depth_ln(wmin, wmax, lr, depth_factor)

                #print(lr, depths)
                hmin = np.cumsum(depths[0])[:-1]
                hmax = (hmin + depths[1][0])
                vsmin = ((1 - fac) * vmin) * np.ones(len(hmin))
                # vsmin=np.array([100,100,100,225,225])
                # vsmax=np.arra
                vsmax = (1 + fac) * vmax * np.ones(len(hmax))
                # vsmax=np.array([250,250,250,350,350])
                # vsmax=450*np.ones(len(hmax))
                ro = 1.9 * np.ones(len(vsmin))
                sigma = 0.40 * np.ones(len(vsmin))
                vp = 1350 * np.ones(len(vsmin))
                bedrock = np.array([vmax*(1-fac), 600, 2.0, 0.45, 1500, 0.0, 0.0]).reshape((7, 1))
                nl = len(hmin) + 1
                nsl = len(hmin)
                a, b = [1290., 1.11]
                fixhv = 9999
                herman = 9999
                weights = [0.10, 0.40]
                vsrep = (vsmin + vsmax) / 2
                width = hmax - hmin
                travel_time = width / (vsrep)
                vsavg = np.sum(width) / (np.sum(travel_time))
                maxdepth = vsavg / (4 * 0.67)
                #print(vsavg / (4 * 0.67))

                header = str(inv_for) + '\n'
                header += str(kind_inv) + '\n'
                header += str(app_fun) + '\n'
                header += str(stdrel) + '\n'
                header += file_disR + '\n'
                header += file_disL + '\n'
                header += file_hv + '\n'
                header += file_test + '\n'
                header += str(seed) + '\n'
                header += str(ga_params[0]) + ' ' + str(ga_params[1]) + ' ' + str(ga_params[2]) + ' ' + str(ga_params[3]) + '\n'
                header += str(nl) + '\n'
                header += str(nsl) + '\n'
                header += "#vsmin(m/s),vsmax(m/s), hmin(m),hmax(m),ro,sigma,vp"

                footer = "#" + '\n'
                footer += str(a) + ' ' + str(b) + '\n'
                footer += str(fixhv) + '\n'
                footer += str(herman) + '\n'
                footer += str(weights[0]) + ' ' + str(weights[1]) + '\n'

                table = np.vstack((vsmin, vsmax, hmin, hmax, ro, sigma, vp))
                table = np.hstack((table, bedrock)).T
                np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/set_up.dat', table, header=header,
                           comments='', footer=footer, fmt='%1.3f')
                # try:
                #     shutil.rmtree('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+qinv_name)
                # except:
                #     print('could not remove tree, check later')
                #qinv_name = cdd_name[:-4] + f'_lr{int(lr * 10)}_' + 'seed'+ str(seed).zfill(2)
                qinv_name = cdd_name[:-4] + f'_ln{int(len(depths[0]))}_' + 'seed'+ str(seed).zfill(2)
                # try:
                #     os.mkdir('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + subfolder + '/' + cdd_name)
                # except:
                #     pass
                # shutil.copytree('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/',
                #                 '/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + subfolder +'/' + cdd_name + '/' + qinv_name)
                shutil.copytree('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/',
                                '/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + subfolder + '/' + qinv_name)
                #print(num)
                # if nummer == 3:
                #     with open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + subfolder + '/script.sh',
                #               'a') as fileobject:
                #         fileobject.write('wait')
                #         fileobject.write('\n')
                #     nummer = 0

                # with open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+ subfolder + '/script.sh', 'a') as fileobject:
                #     fileobject.write('cd ~/inversion/qinvs_tuesday_24/' + cdd_name+ '/' + qinv_name + ' && ./compila.sh')
                #     fileobject.write('\n')
                # if seed == seeds[-1] and lr == lrs[-1]:
                #     nummer+=1


                    #fileobject.close()
                # with open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+ subfolder + '/'+ cdd_name+'/script.sh', 'a') as fileobject:
                #     fileobject.write('cd ' +  qinv_name + ' && ./compila.sh &')
                #     fileobject.write('\n')
                #     fileobject.close()
    # for seed in seeds:
    #     for ln in lns:
    #         depths = swprepost.Parameter.depth_ln(wmin, wmax, ln, depth_factor)
    #         hmin = np.cumsum(depths[0])[:-1]
    #         hmax = (hmin + depths[1][0])
    #         vsmin = ((1-fac)*vmin)*np.ones(len(hmin))
    #         #vsmin=np.array([100,100,100,225,225])
    #         #vsmax=np.arra
    #         vsmax = (1+fac)*vmax*np.ones(len(hmax))
    #         #vsmax=np.array([250,250,250,350,350])
    #         #vsmax=450*np.ones(len(hmax))
    #         ro = 1.9*np.ones(len(vsmin))
    #         sigma=0.40*np.ones(len(vsmin))
    #         vp=1350*np.ones(len(vsmin))
    #         bedrock = np.array([300, 1000, 2.0, 0.45, 1500, 0.0, 0.0]).reshape((7, 1))
    #         nl = len(hmin) + 1
    #         nsl = len(hmin)
    #         a, b = [1290., 1.11]
    #         fixhv = 9999
    #         herman = 9999
    #         weights = [0.10, 0.40]
    #         vsrep = (vsmin + vsmax) / 2
    #         width = hmax - hmin
    #         travel_time = width / (vsrep)
    #         vsavg = np.sum(width) / (np.sum(travel_time))
    #         maxdepth = vsavg / (4 * 0.67)
    #         print(vsavg / (4 * 0.67))
    #
    #         header = str(inv_for) + '\n'
    #         header += str(kind_inv) + '\n'
    #         header += str(app_fun) + '\n'
    #         header += str(stdrel) + '\n'
    #         header += file_disR + '\n'
    #         header += file_disL + '\n'
    #         header += file_hv + '\n'
    #         header += file_test + '\n'
    #         header += str(seed) + '\n'
    #         header += str(ga_params[0]) + ' ' + str(ga_params[1]) + ' ' + str(ga_params[2]) + ' ' + str(ga_params[3]) + '\n'
    #         header += str(nl) + '\n'
    #         header += str(nsl) + '\n'
    #         header += "#vsmin(m/s),vsmax(m/s), hmin(m),hmax(m),ro,sigma,vp"
    #
    #         footer = "#" + '\n'
    #         footer += str(a) + ' ' + str(b) + '\n'
    #         footer += str(fixhv) + '\n'
    #         footer += str(herman) + '\n'
    #         footer += str(weights[0]) + ' ' + str(weights[1]) + '\n'
    #
    #         table = np.vstack((vsmin, vsmax, hmin, hmax, ro, sigma, vp))
    #         table = np.hstack((table, bedrock)).T
    #         np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/set_up.dat', table, header=header,
    #                    comments='', footer=footer, fmt='%1.3f')
    #         # try:
    #         #     shutil.rmtree('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+qinv_name)
    #         # except:
    #         #     print('could not remove tree, check later')
    #         qinv_name = cdd_name+f'_ln{ln}_'+f'seed{seed}'
    #
    #
    #         shutil.copytree('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/',
    #                         '/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' +cdd_name+ '/'+ qinv_name)
    #
    #         with open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + cdd_name+'/script.sh', 'a')  as fileobject:
    #             fileobject.write('cd '+qinv_name+ ' && ./compila.sh &')
    #             fileobject.write('\n')
    #             fileobject.close()
