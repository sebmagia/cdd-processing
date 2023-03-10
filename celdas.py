
import numpy as np
import matplotlib.pyplot as plt
path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO/'

# coords = np.loadtxt(path_cdd + 'cdd_coordinates.dat')
# fig, ax = plt.subplots()
# #dist = ((coords[num, 0] - coords[num, 2]) ** 2 + (coords[num, 1] - coords[num, 3]) ** 2) ** 0.5
# for x in coords:
#     ax.plot([x[0], x[2]], [x[1], x[3]], 'k-^')


def grid(x0,y0,LX,LY,dx,dy):
    x=np.arange(x0,x0+LX+dx,dx)
    y=np.arange(y0,y0+LY+dy,dy)

    return x,y

x0=0.0
y0=0.0
LX=40
LY=20
dx=5
dy=4

x,y=grid(x0,y0,LX,LY,dx,dy)
xx,yy=np.meshgrid(x,y)

P1=(8,16)
P2=(5,18)
plt.plot(xx,yy,'go')
plt.plot(P1,P2,'rv-')
#ax.plot([coords[num, 0], coords[num, 2]], [coords[num, 1], coords[num, 3]], 'r-^')