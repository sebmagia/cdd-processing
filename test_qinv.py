import numpy as np
import matplotlib.pyplot as plt
import shutil
import swprepost

path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO/'
cdd=np.loadtxt(path_cdd+'cdd01.dat')
path_hv='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/HV2/Dic21/HV/'
hv1=np.loadtxt(path_hv+'data_3C_000.hv',skiprows=24,delimiter=',')
hv2=np.loadtxt(path_hv+'data_3C_002.hv',skiprows=24,delimiter=',')

gen_setup=True

amax1=np.argmax(hv1[:,1])
fmin=1.7
fmax=15
f0=0.67
idx=np.where((hv1[:,0] <= cdd[0,0] ) & (hv1[:,0] > fmin))
freq=hv1[idx,0][0]
amp=hv1[idx,1][0]
#plt.figure(1)
# plt.plot(hv1[:,0],hv1[:,1],'r',linewidth=2)
# plt.plot(hv2[:,0],hv2[:,1],'k',linewidth=2)
# plt.plot(freq,amp,'b',linewidth=2)
# plt.xlabel('Frequency')
# plt.ylabel('HVSR')
#
# plt.figure(2)
# plt.plot(cdd[:,0],cdd[:,1])
# plt.figure(3)
# plt.plot(cdd[:,3],cdd[:,1])
#
# plt.xlabel('Frequency')
# plt.ylabel('Velocity')
#np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/dispR',cdd[:,:2])
#np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/hvDushanbe',np.vstack((freq,amp)).T)

np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/dispR',cdd[:,:2])
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
    seed=3
    ga_params=[50,50,0.7,0.01]
    wmin=cdd[0,3]
    wmax=cdd[-1,3]
    vmin=cdd[0,1]
    vmax=cdd[-1,1]
    lns=[3,4,5,6,7]
    lrs=[3.0,2.0,1.5,1.2]
    depth_factor=3

    for ln in lns:
        depths = swprepost.Parameter.depth_ln(wmin, wmax, ln, depth_factor)
        hmin = np.cumsum(depths[0])[:-1]
        hmax = (hmin + depths[1][0])
        vsmin = np.array([100, 100, 100, 100])
        vsmax = np.array([350, 350, 350, 350])
        ro = 1.9*np.ones(len(vsmin))
        sigma=0.40*np.ones(len(vsmin))
        vp=1350*np.ones(len(vsmin))
        bedrock = np.array([300, 760, 2.0, 0.45, 1500, 0.0, 0.0]).reshape((7, 1))
        nl = len(hmin) + 1
        nsl = len(hmin)




    depths = swprepost.Parameter.depth_ln(wmin, wmax, ln, depth_factor)
    hmin = np.cumsum(depths[0])[:-1]
    hmax = (hmin + depths[1][0])
    vsmin = np.array([100, 100, 100,100])
    vsmax = np.array([350, 350, 350,350])
    ro = np.array([1.9, 1.9, 1.9,1.9])
    sigma = np.array([0.35, 0.35, 0.35,0.35])
    vp = np.array([1350, 1350, 1350,1350])
    bedrock=np.array([300,1500,2.0,0.45,1500,0.0,0.0]).reshape((7,1))
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

    # vsbedrock=np.array([400,1700])
    # meanvs = (np.sum(vsmin[:-1]) + np.sum(vsmax[:-1])) / (len(vsmin) + len(vsmax)-2)

    print(vsavg / (4 * 0.67))

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
    qinv_name = 'qinv_48'
    np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/set_up.dat', table, header=header,
               comments='', footer=footer, fmt='%1.3f')
    # try:
    #     shutil.rmtree('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+qinv_name)
    # except:
    #     print('could not remove tree, check later')
    shutil.copytree('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/',
                    '/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + qinv_name)

    # if nlay ==3:
    #     ## param1

    #nlay=3

    #vsmin = np.array([120, 140, 180, 200])
    #vsmax = np.array([150, 200, 250, 350])
    #vsmin = np.array([120, 160, 210, 240,300])

    #vsmax = np.array([180, 230, 270, 320,450])
    # hmin = np.array([3, 6, 10, 20, 30, 50])
    # vsmin = np.array([120, 150, 180, 200, 240,300])
    #
    # hmax = np.array([8, 15, 25, 40, 60, 100])
    # vsmax = np.array([160, 200, 230, 260, 300,400])
#vsmax = np.array([160, 225, 250, 300])
    #
    # ro = np.array([1.9,1.9,1.9,1.9])
    # sigma = np.array([0.35, 0.35,0.35,0.35,0.35])
    # vp = np.array([1350, 1350,1350,1350,1350])
    #
    # bedrock=np.array([300,1500,2.0,0.45,1500,0.0,0.0]).reshape((7,1))
    # nl = len(hmin) + 1
    # nsl = len(hmin)
    #
    # a, b = [1290., 1.11]
    # fixhv = 9999
    # herman = 9999
    # weights = [0.10, 0.40]
    # vsrep=(vsmin+vsmax)/2
    # width = hmax - hmin
    # travel_time = width / (vsrep)
    # vsavg=np.sum(width)/(np.sum(travel_time))
    # maxdepth=vsavg/(4*0.67)
    #
    # #vsbedrock=np.array([400,1700])
    # #meanvs = (np.sum(vsmin[:-1]) + np.sum(vsmax[:-1])) / (len(vsmin) + len(vsmax)-2)
    #
    # print(vsavg / (4 * 0.67))
    #
    # header=str(inv_for)+'\n'
    # header+=str(kind_inv)+'\n'
    # header+=str(app_fun)+'\n'
    # header+=str(stdrel)+'\n'
    # header+=file_disR+'\n'
    # header+=file_disL+'\n'
    # header+=file_hv+'\n'
    # header+=file_test+'\n'
    # header+=str(seed)+'\n'
    # header+=str(ga_params[0])+' '+str(ga_params[1])+' '+str(ga_params[2])+' '+str(ga_params[3])+'\n'
    # header+=str(nl)+'\n'
    # header+=str(nsl)+'\n'
    # header+="#vsmin(m/s),vsmax(m/s), hmin(m),hmax(m),ro,sigma,vp"
    #
    # footer="#"+'\n'
    # footer+=str(a)+' '+str(b)+'\n'
    # footer+=str(fixhv)+'\n'
    # footer+=str(herman)+'\n'
    # footer+=str(weights[0])+' '+str(weights[1])+'\n'
    #
    #
    #
    # table = np.vstack((vsmin, vsmax, hmin, hmax, ro, sigma, vp))
    # table = np.hstack((table, bedrock)).T
    # qinv_name='qinv_46'
    # np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/set_up.dat',table,header=header,comments='',footer=footer,fmt='%1.3f')
    # # try:
    # #     shutil.rmtree('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+qinv_name)
    # # except:
    # #     print('could not remove tree, check later')
    # shutil.copytree('/home/doctor/Doctor/Magister/Tesis/databases/process_data/qinv/','/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+qinv_name)

    # with  open('/home/doctor/Doctor/Magister/Tesis/databases/qinvs/test.csv', "ab") as f:
    #     f.write(b"\n")
    #     np.savetxt(f,table, fmt='%1.3f',delimiter=',',header='vsmin(m/s),vsmax(m/s), hmin(m),hmax(m),ro,sigma,vp')

    #np.savetxt('/home/doctor/Doctor/Magister/Tesis/databases/qinv/set_up.dat',table,header=header,comments='',footer=footer,fmt='%1.3f')
    # with open('set_up.dat', 'w') as fileobject:
    #     fileobject.write(str(inv_for))
    #     fileobject.write('\n')
    #     fileobject.write(str(kind_inv))
    #     fileobject.write('\n')
    #     fileobject.write(str(app_fun))
    #     fileobject.write('\n')
    #     fileobject.write(str(stdrel))
    #     fileobject.write('\n')
    #     fileobject.write(file_disR)
    #     fileobject.write('\n')
    #     fileobject.write(file_disL)
    #     fileobject.write('\n')
    #     fileobject.write(file_hv)
    #     fileobject.write('\n')
    #     fileobject.write(file_test)
    #     fileobject.write('\n')
    #     fileobject.write(str(seed))
    #     fileobject.write('\n')
    #     fileobject.write(str(ga_params[0])+' '+str(ga_params[1])+' '+str(ga_params[2])+' '+str(ga_params[3]))
    #     fileobject.write('\n')
    #     fileobject.write(str(nl))
    #     fileobject.write('\n')
    #     fileobject.write(str(nsl))
    #     fileobject.write('\n')
    #     i=0
    #     while i<nl-1:
    #         fileobject.write(str(vsmin[i])+' '+str(vsmax[i])+' '+str(hmin[i])+' '+str(hmin[i])+' '+str(hmax[i])+' '+str(ro[i])+' '+str(sigma[i])+' '+str(vp[i])+' ')
    #         fileobject.write('\n')
    #         i+=1
    #     fileobject.write(str(vsmin[i])+' '+str(vsmax[i])+' '+'000'+' '+'000'+' '+str(ro[i])+' '+str(sigma[i])+' '+str(vp[i]))





#fmax1=hv1[amax1,0]

#amax2=np.argmax(hv2[:,1])
#fmax2=hv2[amax2,0]

#idx=np.where((hv1[:,0] < fmax1 +0.25) & (hv1[:,0] > fmax1 -0.25))

freq=hv1[idx,0][0]
amp=hv1[idx,1][0]
plt.figure(1)
plt.plot(hv1[:,0],hv1[:,1],'r',linewidth=2)
plt.plot(hv2[:,0],hv2[:,1],'k',linewidth=2)
plt.plot(freq,amp,'b',linewidth=2)
plt.xlabel('Frequency')
plt.ylabel('HVSR')

plt.figure(2)
plt.plot(cdd[:,0],cdd[:,1])
plt.xlabel('Frequency')
plt.ylabel('Velocity')
np.savetxt('dispR',cdd[:,:2])
np.savetxt('hvDushanbe',np.vstack((freq,amp)).T)
#
# vsmin=np.array([120,150,200,270])
# vsmax=np.array([200,250,350,450])
# hmin=np.array([5,5,15,30])
# hmax=np.array([15,30,50,80])
#
# meanvs=(np.sum(vsmin)+np.sum(vsmax))/(len(vsmin)+len(vsmax))
# print(meanvs/(4*0.67))