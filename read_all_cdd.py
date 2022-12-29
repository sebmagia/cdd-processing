import numpy as np
import swprepost
import os
import matplotlib.pyplot as plt
folders=['AJUSTES_DIC21/','AJUSTES_HS21/','AJUSTES_MAR22/','AJUSTES_SEPT22/']
path_tomo='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO'
def plot_target(target):
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6, 3), dpi=150)
    target.plot(x="frequency", y="velocity", ax=axs[0])
    target.plot(x="wavelength", y="velocity", ax=axs[1])
    axs[1].set_ylabel("")
    axs[1].legend()
    return (fig, axs)

resample=True
#all_models=True
n=0
fig, ax = plt.subplots()

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
        for x in cdds:
            f=open(path+x[:-3]+'csv')
            line=f.readlines()[1]
            coords_a_x=np.round(float(line.split(',')[1][1:]),2)
            coords_a_y=np.round(float(line.split(',')[2][:-1]),2)
            coords_b_x=np.round(float(line.split(',')[3][1:]),2)
            coords_b_y=np.round(float(line.split(',')[4][:-2]),2)
            stack=np.hstack((coords_a_x, coords_a_y, coords_b_x, coords_b_y))
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
            np.savetxt(path_tomo+'/cdd'+str(n).zfill(2)+'.dat',fv)
            n+=1

            ax.plot(fv[:,0] , fv[:,1], 'k', linewidth=0.5)
            #ax.text(wnews[i][-1] / (2 * np.pi), cps[i][-1], 'c' + str(i).zfill(2), fontsize=10)
            #target.to_csv(path+x[:-4]+'_rs.csv')
            #fig, axs = plot_target(target)

#for c in coords:
 #   plt.plot([c[0],c[2]],[c[1],c[3]],'r^-')
fig,ax=plot_target(targets[0])