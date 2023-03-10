import matplotlib.pyplot as plt
import numpy as np
import os
path_gm='/home/doctor/Doctor/Magister/Tesis/databases/process_data/GROUND_MODELS_FEB/'
path_wm='/home/doctor/Doctor/Magister/Tesis/databases/process_data/WEIGHTED_MODELS_FEB/'
gms=sorted(os.listdir(path_gm))
gms=[x for x in  gms if 'SIGMA' not in x]
gms_sigma=sorted(os.listdir(path_gm))
gms_sigma=[x for x in gms_sigma if 'SIGMA' in x]

sigma1=np.loadtxt(path_gm+gms_sigma[0])
sigma2=np.loadtxt(path_gm+gms_sigma[1])
sigma12=np.vstack((sigma1[:,1],sigma2[:,1])).T
fig,ax=plt.subplots()
ax.plot(sigma1[:,1],sigma1[:,0],'r')
ax.plot(sigma2[:,1],sigma2[:,0],'b')

