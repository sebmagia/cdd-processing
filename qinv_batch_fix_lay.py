import numpy as np
import matplotlib.pyplot as plt
import shutil
import swprepost
import os
import matplotlib


path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO/'
path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/INV/'

cdd_names=os.listdir(path_cdd)
cdd_names=sorted([x for x in cdd_names if '.dat' in x and 'coordinates' not in x])
os.listdir('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/qinvs_feb15')
coords=np.loadtxt(path_cdd+'cdd_coordinates.dat')
subfolder='qinvs_feb15_fl'
hv = False
gen_setup=True
## read number of layers file
lrs=np.loadtxt(path_cdd+'numberoflayers.txt',dtype=int)
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
        seeds=[3,5,7,9,11,13,15,17,19,21,23,25]
        #seeds=[3,5,7]
        #seeds=[5]
        ga_params=[100,100,0.7,0.01]
        ga_params=[100,100,0.7,0.01]
        wmin=cdd[0,1]/cdd[0,0]
        wmax=cdd[-1,1]/(cdd[-1,0])
        vmin=cdd[0,1]
        vmax=cdd[-1,1]

        depth_factor=2
        fac=0.15 ## factor por minmax velocities
        #depths_lr = swprepost.Parameter.depth_lr(wmin, wmax, lrs[3], depth_factor)

        lrs_trial=np.arange(1,10,0.05)[1:]
        trials=[]


        for trial in lrs_trial:
            try:
                depths = swprepost.Parameter.depth_lr(wmin, wmax, trial, depth_factor)
            except IndexError:
                #print('aq')
                break

            nl = len(depths[0])
            #print('nl',nl)
            if nl == lrs[num]:
                trials.append(trial)
            # if nl > lrs[num]:
            #     break
        lr=np.median(trials)

        lr=np.round(lr,2)
        print('cdd',num)
        for seed in seeds:
            depths = swprepost.Parameter.depth_lr(wmin, wmax, lr, depth_factor)
            hmin = np.cumsum(depths[0])[:-1]
            hmax = (hmin + depths[1][0])
            vsmin = ((1 - fac) * vmin) * np.ones(len(hmin))
            vsmax = (1 + fac) * vmax * np.ones(len(hmax))
            ro = 1.9 * np.ones(len(vsmin))
            sigma = 0.40 * np.ones(len(vsmin))
            vp = 1350 * np.ones(len(vsmin))
            bedrock = np.array([vmax*(1-fac), 450, 2.0, 0.45, 1500, 0.0, 0.0]).reshape((7, 1))
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

            qinv_name = cdd_name[:-4] + f'_ln{int(len(depths[0]))}_' + 'seed'+ str(seed).zfill(2)

            shutil.copytree('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/',
                            '/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + subfolder + '/' + qinv_name)
