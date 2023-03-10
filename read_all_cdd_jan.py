import numpy as np
import swprepost
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

#folders=['AJUSTES_DIC21/','AJUSTES_HS21/','AJUSTES_MAR22/','AJUSTES_SEPT22/']
folder='AJUSTES_NUEVOS_GS_FINAL/'
folders=['dic21/','hs21/','mar22/','sept22/','dic22/']
folders=[folder + x for x in folders]
path_tomo='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/'
path_estado = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/'
estado=np.loadtxt(path_estado+'estado.txt')
estado=np.loadtxt(path_estado+'estados_nuevos.txt')

distancias=np.loadtxt(path_estado+'distancias.txt')

#path_estado_2 = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_JAN/estado_2.txt'


def plot_target(target):
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6, 3), dpi=150)
    target.plot(x="frequency", y="velocity", ax=axs[0])
    target.plot(x="wavelength", y="velocity", ax=axs[1])
    axs[1].set_ylabel("")
    axs[1].legend()
    return (fig, axs)

resample=False
check_reliable=True
merge_clones=True
check_by_single=True
todas=True
check_by_location=False
check_by_clone=False
n=0

obtain_coords=False
if not obtain_coords:
    coords = np.loadtxt(path_tomo + '/cdd_coordinates.dat')
    coords_sp=np.split(coords,2,axis=1)
    coords_uni=np.vstack((coords_sp[0],coords_sp[1]))

if resample:
    for folder in folders:
        # Approach 1: Import from comma seperated text file (see swprepost documentation for details).
        path='/home/doctor/Doctor/Magister/Tesis/databases/process_data/'+folder
        cdds=os.listdir(path)
        cdds=[x for x in cdds if '.csv' in x]
        cdds=sorted(cdds)
        #data=np.loadtxt(path+cdds[0])
        #np.savetxt(path+cdds[0][:-3]+'csv',data,header='Frequency,Velocity,Velstd',delimiter=',')
        targets=[]
        #k=0
        for x in cdds:
            if not obtain_coords:
                fig5, ax5 = plt.subplots(2, 1)
                ax5[1].plot(coords_uni[:, 0], coords_uni[:, 1], 'b^')
                stas = coords[n]
                dist = np.round(((stas[0] - stas[2]) ** 2 + (stas[1] - stas[3]) ** 2) ** 0.5, 2)
            if obtain_coords:
                path_ori='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/'
                coordenadas_todas = np.loadtxt(path_ori + 'coords_all_sets.txt')
                coord_split = np.split(coordenadas_todas, 2, axis=1)
                coord_uni = np.unique(np.vstack((coord_split[0], coord_split[1])), axis=0)
                fig5, ax5 = plt.subplots(2, 1)
                f=open(path+x[:-3]+'csv')
                line=f.readlines()[1]
                coords_a_x=np.round(float(line.split(',')[1][1:]),2)
                coords_a_y=np.round(float(line.split(',')[2][:-1]),2)
                coords_b_x=np.round(float(line.split(',')[3][1:]),2)
                coords_b_y=np.round(float(line.split(',')[4][:-2]),2)
                stack=np.hstack((coords_a_x, coords_a_y, coords_b_x, coords_b_y))
                dist = np.round(((stack[0] - stack[2]) ** 2 + (stack[1] - stack[3]) ** 2) ** 0.5, 2)
                print(dist)
                try:
                    coords=np.vstack((coords,stack))
                except:
                    coords=stack


            target = swprepost.Target.from_csv(path+x[:-3]+'csv')
            #fig, axs = plot_target(target)
            domain = "wavelength"       # "frequency" or "wavelength", "wavelength" is recommended
            res_type = "log"            # "log" or 'linear', "log" is recommended.
            pmin = target.wavelength[-1]                  # Minimum value after resampling in units of domain
            pmax = target.wavelength[0]                  # Maximum value after resampling in units of domain
            pn = 20                     # Number of samples, 20-30 points are recommended.
            target.easy_resample(pmin=pmin, pmax=pmax, pn=pn, res_type=res_type, domain=domain, inplace=True)
            targets.append(target)
            freq=target.frequency
            vel=target.velocity
            velstd=target.velstd
            wv=target.wavelength
            fv=np.round(np.vstack((freq,vel,velstd,wv)).T,3)
            np.savetxt(path_tomo+'/cdd'+str(n).zfill(3)+'.dat',fv)
            np.savetxt(path_tomo+'/cdd'+str(n).zfill(3)+'.csv',fv[:,:3],header='Frequency,Velocity,Velstd',delimiter=',')

            ax5[0].plot(fv[:, 0], fv[:, 1], 'k', linewidth=2)
            ax5[0].plot(fv[:, 0], fv[:, 1] + 1.96 * fv[:, 2], '--r', linewidth=2)
            ax5[0].plot(fv[:, 0], fv[:, 1] - 1.96 * fv[:, 2], '--r', linewidth=2)
            ax5[0].set_xlabel('Frequency [Hz]')
            ax5[0].set_ylabel('Phase Velocity [m/s]')
            if not obtain_coords:
                ax5[1].plot([coords[n, 0], coords[n, 2]], [coords[n, 1], coords[n, 3]], 'r-^')
            if obtain_coords:
                ax5[1].plot([stack[0], stack[2]], [stack[1], stack[3]], 'r-^')

            fig5.suptitle('Interstation Distance ' + str(dist))
            plt.savefig(path_tomo + '/cdd_' + str(n).zfill(3) + '.png', format='png', dpi=300, bbox_inches='tight')
            plt.close()
            n += 1
if obtain_coords:
    np.savetxt(path_tomo+'/cdd_coordinates.dat',coords)

if check_reliable:
    files=os.listdir(path_tomo)
    files=sorted([x for x in files if '.dat' in x])[:-2]
    estados=np.loadtxt(path_estado+'estados_nuevos.txt').astype(int)
    for i in range(len(files)):
            data=np.loadtxt(path_tomo+'/'+files[i])
            try:
                 freq_stack=np.vstack((freq_stack,data[:,0]))
                 vel_stack=np.vstack((vel_stack,data[:,1]))
                 vel_std_stack=np.vstack((vel_std_stack,data[:,2]))
                 wv_stack=np.vstack((wv_stack,data[:,3]))


            except:
                freq_stack = data[:,0]
                vel_stack = data[:, 1]
                vel_std_stack = data[:, 2]
                wv_stack = data[:, 3]
    fig,ax=plt.subplots(2,1,figsize=(10,8))
    for i in range(len(files)):
    #for i in range(65,79):
        if estados[i]==1:
            ax[0].plot(freq_stack[i],vel_stack[i],'g-')
            ax[1].plot(wv_stack[i], vel_stack[i], 'g-')

        if estados[i] == 2:
            ax[0].plot(freq_stack[i],vel_stack[i],'y-')
            ax[1].plot(wv_stack[i], vel_stack[i], 'y-')
        if estados[i] == 3:
            ax[0].plot(freq_stack[i], vel_stack[i], 'r-')
            ax[1].plot(wv_stack[i], vel_stack[i], 'r-')

        # if estados[i] == 0:
        #     ax[0].plot(freq_stack[i],vel_stack[i],'r-')
        #     ax[1].plot(wv_stack[i], vel_stack[i], 'r-')
    #ax[0].set_xlim([1.5,30])
    ax[0].set_ylim([100,400])
    ax[1].set_ylim([100,400])
    ax[1].set_xlim([-5,200])
    ax[0].set_ylabel('Velocity [m/s]')
    ax[1].set_ylabel('Velocity [m/s]')
    ax[0].set_xlabel('Frequency [Hz]')

    ax[1].set_xlabel('Wavelength [m]')
    #plt.savefig( )
    plt.savefig(path_tomo+'/all_dispersion_curves.png', format='png', dpi=300, bbox_inches='tight', pad_inches=0.1)
if merge_clones:
    coords_sorted=np.loadtxt(path_tomo+'/cdd_coordinates_sorted.dat')
    f1,c1=freq_stack[0],vel_stack[0]
    f2,c2=freq_stack[1],vel_stack[1]
    f3,c3=freq_stack[83],vel_stack[83]

    x1a,x2a,y1a,y2a=coords_sorted[0,0],coords_sorted[0,2],coords_sorted[0,1],coords_sorted[0,3]
    x1b,x2b,y1b,y2b=coords_sorted[1,0],coords_sorted[1,2],coords_sorted[1,1],coords_sorted[1,3]
    x1c,x2c,y1c,y2c=coords_sorted[83,0],coords_sorted[83,2],coords_sorted[83,1],coords_sorted[83,3]

    midx = [np.mean([ f1[i], f2[i], f3[i]  ]) for i in range(len(f1))]
    midy = [np.mean([ c1[i], c2[i], c3[i]  ]) for i in range(len(f1))]
    cx = [np.mean([f1[i], f2[i], f3[i]]) for i in range(len(f1))]
    cy = [np.mean([c1[i], c2[i], c3[i]]) for i in range(len(f1))]
    #plt.plot(f1,c1,'b-')
    #plt.plot(f2,c2,'b-')
    #plt.plot(f3,c3,'b-')

    #plt.plot(midx,midy,'r-')
    lineas = []
    with open('clones.txt', 'r') as fileobject:
        lines = fileobject.readlines()
        for line in lines:
            linea = line.strip().split('\t')
            lista = [int(x) for x in linea]
            lineas.append(lista)
    clones_new=[]
    for clones in lineas:
        clones=[x for x in clones if estados[x] != 3 ]
        clones_new.append(clones)
    for clones in clones_new:
        midx=[np.mean(freq_stack[clones,i]) for i in range(freq_stack.shape[1])]
        midy=[np.mean(vel_stack[clones,i]) for i in range(freq_stack.shape[1])]
        midstd = [np.max(vel_std_stack[clones,i]) for i in range(freq_stack.shape[1])]

        cx1=np.mean(coords_sorted[clones,0])
        cx2=np.mean(coords_sorted[clones,2])
        cy1 = np.mean(coords_sorted[clones, 1])
        cy2 = np.mean(coords_sorted[clones, 3])
        merg=[cx1,cy1,cx2,cy2]
        try:
            coords_merged=np.vstack((coords_merged,merg))
        except:
            coords_merged=np.array(merg)

        try:
            freqs_merged=np.vstack((freqs_merged,midx))
            vels_merged=np.vstack((vels_merged,midy))
            vels_std_merged=np.vstack((vels_std_merged,midstd))

        except NameError:
            freqs_merged=midx
            vels_merged=midy
            vels_std_merged=midstd
    plt.figure(33)
    for i in range(len(coords_merged)):
        plt.plot([coords_merged[i, 0], coords_merged[i, 2]], [coords_merged[i, 1], coords_merged[i, 3]], '-ok')
    n=0
    plt.figure()
    for f,v in zip(freqs_merged,vels_merged):
        plt.plot(f,v,'b')
        plt.text(f[0],v[0],str(clones_new[n][0]))
        n+=1
        #print(clones,clones_new)


if check_by_single:
    lineas=[]
    with open('clones.txt', 'r') as fileobject:
        lines = fileobject.readlines()
        for line in lines:
            linea = line.strip().split('\t')
            lista = [int(x) for x in linea]
            for elem in lista:
                lineas.append(elem)
    clonadas=sorted(set(lineas))
    fig, ax = plt.subplots()
    q=0
    for i in range(len(files)):
        if i not in clonadas:
            if estados[i]==1:
                ax.plot(freq_stack[i], vel_stack[i],'g')
                ax.text(freq_stack[i, 0], vel_stack[i, 0], str(i))
                q+=1
                freqs_merged = np.vstack((freqs_merged, freq_stack[i]))
                vels_merged = np.vstack((vels_merged, vel_stack[i]))
                coords_merged = np.vstack((coords_merged, coords_sorted[i]))
                vels_std_merged=np.vstack((vels_std_merged,vel_std_stack[i]))
            if estados[i] == 2:
                ax.plot(freq_stack[i], vel_stack[i], 'r')
                ax.text(freq_stack[i, 0], vel_stack[i, 0], str(i))
                q+=1
                freqs_merged = np.vstack((freqs_merged, freq_stack[i]))
                vels_merged = np.vstack((vels_merged, vel_stack[i]))
                coords_merged = np.vstack((coords_merged, coords_sorted[i]))
                vels_std_merged=np.vstack((vels_std_merged,vel_std_stack[i]))

        print(q)

if todas:
    fig3,ax3=plt.subplots()
    for i in range(len(freqs_merged)):

        fig1, ax1 = plt.subplots(figsize=(10,5))
        #ax[0].plot(freqs_merged[i],vels_merged[i],'k-')
        ax1.plot(freqs_merged[i], vels_merged[i], 'k', linewidth=2)
        ax1.plot(freqs_merged[i], vels_merged[i] + 1.96 * vels_std_merged[i], '--r', linewidth=2)
        ax1.plot(freqs_merged[i], vels_merged[i] - 1.96 * vels_std_merged[i], '--r', linewidth=2)
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel('Phase Velocity [m/s]')
        ax1.set_ylim([100,325])



        dist=((coords_merged[i, 0] - coords_merged[i, 2])**2 + (coords_merged[i, 1] - coords_merged[i, 3])**2)**0.5
        fig1.suptitle('Interstation Distance ' + str(np.round(dist,2)))
        plt.savefig(path_tomo + '/INV/cdd_' + str(i).zfill(3) + '.png', format='png', dpi=300, bbox_inches='tight')
        plt.close()

        fig2, ax2 = plt.subplots()
        ax2.set_xlim([-5,30])
        ax2.set_ylim([-24,30])
        ax2.set_xlabel('X Distance [m]')
        ax2.set_ylabel('Y Distance [m]')
        ax2.plot([coords_merged[i, 0], coords_merged[i, 2]], [coords_merged[i, 1], coords_merged[i, 3]], 'r-^', zorder=2)
        for j in range(len(freqs_merged)):
            ax2.plot([coords_merged[j, 0], coords_merged[j, 2]], [coords_merged[j, 1], coords_merged[j, 3]], 'y^',
                     zorder=1)

        np.savetxt(path_tomo + '/INV/cdd_' + str(i).zfill(3) + '.dat',np.vstack((freqs_merged[i],vels_merged[i])).T)
        fig2.suptitle('Interstation Distance ' + str(np.round(dist,2)))
        plt.savefig(path_tomo + '/INV/cdd_' + str(i).zfill(3) + '_coords.png', format='png', dpi=300, bbox_inches='tight')
        plt.close()

        ax3.plot(freqs_merged[i], vels_merged[i], 'b', linewidth=1.5)
        ax3.set_xlabel('Frequency [Hz]')
        ax3.set_ylabel('Phase Velocity [m/s]')
        ax3.set_ylim([100, 325])
    plt.savefig(path_tomo + '/INV/cdd_' + str(i).zfill(3) + '_todas.png', format='png', dpi=300, bbox_inches='tight')

        #ax[0].plot(freqs_merged[i],vels_merged[i],'k-')
        #ax[1].plot([coords_merged[i, 0], coords_merged[i, 2]], [coords_merged[i, 1], coords_merged[i, 3]],'rv-')
    np.savetxt(path_tomo+'/INV/cdd_coordinates.dat',coords_merged)

if check_by_clone:
    clones=[]
    lineas=[]
    estados = estados[np.nonzero(estados)]
    np.savetxt(path_estado+'estados_nuevos.txt',estados,fmt='%i')
    with open('clones.txt','r') as fileobject:
        lines=fileobject.readlines()
        for line in lines:
            linea=line.strip().split('\t')
            lista=[int(x) for x in linea]
            lineas.append(lista)

    for clones in lineas[0:10]:
        fig,ax=plt.subplots(2,1)
        for curva in clones:
            ax[0].plot(freq_stack[curva],vel_stack[curva])
            ax[0].text(freq_stack[curva,0],vel_stack[curva,0],str(curva))
            ax[1].plot([coords[curva,0],coords[curva,2]],[coords[curva,1],coords[curva,3]],'rv-')

if check_by_location:
    coords = np.loadtxt(path_tomo + '/cdd_coordinates.dat')
    coords_new=np.zeros(coords.shape)
    fig,ax=plt.subplots()
    cdds=os.listdir(path_tomo)
    cdds=sorted([x for x  in cdds if '.dat' in x and '.txt' not in x])[:-1]
    cid=[x[3:6] for x in cdds]
    ax.plot(coords_uni[:,0],coords_uni[:,1],'ro-')
    estados=estados[np.nonzero(estados)]
    for n,x in enumerate(coords):

            if x[2] > x[0]:
                coords_new[n]=coords[n]

            if x[2] < x[0]:
                coords_new[n, 0] = coords[n, 2]
                coords_new[n, 1] = coords[n, 3]
                coords_new[n, 2] = coords[n, 0]
                coords_new[n, 3] = coords[n, 1]

            if x[2]==x[0]:
                if x[3]>x[1]:
                    coords_new[n] = coords[n]
                if x[1]<x[3]:
                    coords_new[n, 0] = coords[n, 2]
                    coords_new[n, 1] = coords[n, 3]
                    coords_new[n, 2] = coords[n, 0]
                    coords_new[n, 3] = coords[n, 1]
    np.savetxt(path_tomo+'/cdd_coordinates_sorted.dat',coords_new)
