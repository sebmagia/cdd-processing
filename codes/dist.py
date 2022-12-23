import numpy as np
import itertools
import matplotlib.pyplot as plt
def dist(L):
    iters = list(itertools.combinations(L, 2))
    dists=np.zeros(len(iters))

    for n,x in enumerate(iters):
        print(x)
        dists[n]=np.linalg.norm(np.array(x[0])-np.array(x[1]))

    return dists,iters

## 0.0  2 metros corrido del punto 1 original
v1=(0,0)
v2=(4,0)
v3=(11,2)
v4=(1,6)
v5=(6,7)
v6=(13,8)
v7=(1,11)
v8=(8,12)
v9=(12,10)
v10=(2,16)
v11=(6,17)
v12=(12,15)
L=[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12]
plt.plot(np.array(L)[:,0],np.array(L)[:,1],'go')
d,xx=dist(L)
# v1=(4,0)
# v2=(9,1)
# v3=(16,2)
# v4=(5,4)
# v5=(10,7)
# v6=(17,7)
# v7=(3,10)
# v8=(11,11)
# v9=(15,10)
# v10=(3,14)
# v11=(9,16)
# v12=(18,15)
# v13=(5,20)
# v14=(12,21)
# v15=(17,18)
#v16=(19,19)
# v1=(3,3)
# v2=(8,4)
# v3=(13,3)
# v4=(19,5)
# v5=(5,7)
# v6=(9,7)
# v7=(13,8)
# v8=(10,7)
# v9=(4,14)
# v10=(10,13)
# v11=(16,15)
# v12=(20,16)
# v13=(4,19)
# v14=(9,18)
# v15=(15,20)
# v16=(19,19)
#L=[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15]
#plt.plot(np.array(L)[:,0],np.array(L)[:,1],'go')
#d,xx=dist(L)