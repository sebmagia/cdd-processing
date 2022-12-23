import pandas as pd
import numpy as np
import os
path='AJUSTES/'
files=sorted(os.listdir(path))
files=[x for x in files if '.dat' in x]

for f in files:
    data=np.loadtxt(path+f)
    #df=pd.DataFrame(data,columns=['Frequency','Velocity','Velstd'])
    df=pd.DataFrame(data)

    df.to_csv((path+f[:-4]+'.csv'),header=['#Frequency','Velocity','Velstd'],index=False)