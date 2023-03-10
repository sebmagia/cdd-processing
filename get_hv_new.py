import simple_hvsrpy_interface as shv
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# path='Sept22/SSHV/'
path='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/HV2/'

step1=True
step2=False
# files=os.listdir(path)
# files=[x for x in files if '.mseed' in x]
# freqs=[]
# amps=[]
if step1:
    #folders=['Dic21/SSHV/','Mar22/SSHV/','Sept22/SSHV/']
    folders=['Dic21/']
    for fold in folders:
        path=path+fold
        files = os.listdir(path)
        files = sorted([x for x in files if '.mseed' in x])
        for x in files:
            print(files)
            hv,rel,clar=shv.do_hvsr(path+x)
        #freqs.append(hv.mc_peak_frq('lognormal'))
        #amps.append(hv.mc_peak_amp('lognormal'))
if step2:
    folders=['Dic21/SSHV/','Mar22/SSHV/','Sept22/SSHV/']
    for fold in folders:
        path=fold
        names=[]
        xx=[]
        yy=[]
        hv_okf=sorted(os.listdir(path+'Results'))
        with open(path+'single_station_coords.txt','r') as fileobject:
            for row in fileobject.readlines():
                names.append(row.split(',')[0])
                xx.append(float(row.split(',')[1]))
                yy.append(float(row.split(',')[2]))
        #fig,ax=plt.subplots()
        medianf0=[]
        stdf0=[]
        cxs=[]
        cys=[]
        for hv in hv_okf:
            data=np.loadtxt(path+'Results/'+hv,skiprows=23,delimiter=',')
            with open(path+'Results/'+hv,'r') as fileobject:
                for n,row in enumerate(fileobject.readlines()):
                    if n==13:
                        mf0=float(row.split(',')[1])
                        medianf0.append(mf0)
                    if n==14:
                        st0=float(row.split(',')[1])
                        stdf0.append(st0)
            idx=[n for n, x in enumerate(names) if x == hv.split('.')[0]][0]
            cxs.append(xx[idx])
            cys.append(yy[idx])
            try:
                median=np.vstack((median,data[:,1]))
                std_low=np.vstack((std_low,data[:,2]))
                std_upp=np.vstack((std_upp,data[:,3]))
            except:
                print('bru')
                median=data[:,1]
                std_low=data[:,2]
                std_upp=data[:,3]
            #ax.plot(data[:,0],data[:,1],'r')
            #ax.plot(data[:,0],data[:,2],'k--')
            #ax.plot(data[:,0],data[:,3],'k--')
        print(len(cxs))
        tocsv=np.vstack((cxs,cys,medianf0,stdf0)).T
        print(tocsv.shape)
        freqs=data[:,0]
        median=median.T
        argmax=np.argmax(median,axis=0)
        fmax=freqs[argmax]
        #ampmax=
        std_low=std_low.T
        std_upp=std_upp.T
        try:
            medians=np.hstack((medians,median))
            stdups=np.hstack((stdups,std_upp))
            stdlows=np.hstack((stdlows,std_low))

            csv_data=np.vstack((csv_data,tocsv))
        except:
            csv_data=tocsv
            medians=median
            stdups=std_upp
            stdlows=std_low


    fig,ax=plt.subplots()
    for i in range(medians.shape[1]):
        ax.semilogx(freqs,medians[:,i],'r',linewidth=1,label='median curve')
        ax.semilogx(freqs,stdups[:,i],'k--',linewidth=1,label='median+1std')
        ax.semilogx(freqs,stdlows[:,i],'k--',linewidth=1,label='median-1std')
    meanmed=np.mean(medians,axis=1)
    hv_output=np.vstack((freqs,meanmed)).T
    f0=hv_output[np.argmax(hv_output[:,1]),0]

    data_inver=hv_output[np.where((hv_output[:,0]>= f0-0.2) & (hv_output[:,0]<= f0+0.2) )]
    pathf='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/HV/'

    np.savetxt(pathf+'hvDushanbe.txt',data_inver)
    ax.set_xlabel('Frequency')
    ax.set_ylabel('HVSR Amplitude')


    plt.savefig('HVSR_ALL.png', format='png', dpi=300, bbox_inches='tight')

    df=pd.DataFrame(csv_data, columns=['x', 'y', 'mean', 'stddev'])
    df.to_csv('LP_hv.csv',index=False)
    #c1=[-5,25,25,-5]
    c1=[-5,-5,25,25]
    c2=[25,-15,-15,25]
    cc=np.vstack((c1,c2)).T
    plt.tricontourf(csv_data[:, 0], csv_data[:, 1], csv_data[:, 2])
    # dflims=pd.DataFrame(cc, columns=['x', 'y'])
    # dflims.to_csv('LP_lims.csv',index=False)
    # points=csv_data[:,0:2]
    # points_g=pd.read_csv('sp_centreport.csv')
    # points_g=np.vstack((np.array(points_g.x),np.array(points_g.y))).T
    # from scipy.spatial import Voronoi,voronoi_plot_2d
    # vor = Voronoi(points)
    # fig = voronoi_plot_2d(vor)
    #
    # vor_g=Voronoi(points_g)
    # fig=voronoi_plot_2d(vor_g)
        #coords=np.load(path+'single_station_coords.txt')