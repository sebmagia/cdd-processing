
import numpy as np
import matplotlib.pyplot as plt
fs_SA='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Dic21/frecuencias_SA.txt'
cp_SA='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Dic21/curvas_SA.txt'

fs_SC='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Mar22/frecuencias_SC.txt'
cp_SC='/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/MSEEDS/Mar22/curvas_SC.txt'

fsA=np.loadtxt(fs_SA)
fsC=np.loadtxt(fs_SC)
cp_SA=np.loadtxt(cp_SA)
cp_SC=np.loadtxt(cp_SC)

for i in range(len(fsA)):
    plt.plot(fsA[:,i], cp_SA[:,i], 'r', linewidth=0.5)

for i in range(len(fsC)):
    plt.plot(fsC[:,i], cp_SA[:,i], 'k', linewidth=0.5)
plt.xlabel('Frequency [Hz]')
plt.ylabel('$C_R(f)$')
plt.savefig('ALL_AJUSTES.png', format='png', dpi=300, bbox_inches='tight')
