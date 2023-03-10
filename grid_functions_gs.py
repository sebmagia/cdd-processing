import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as mpl
import os
from shapely.geometry import Point, Polygon, LineString, box
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely.figures import SIZE, BLUE, GRAY, RED, YELLOW, BLACK, set_limits
import swprepost
#import matplotlib
#from grid_functions import create_grid

path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS/'
#coords = np.loadtxt(path_cdd + 'cdd_coordinates.dat')
coords=np.loadtxt(path_cdd+'cdd_coordinates.dat')
lines=[]
rectangles=[]
#fig, ax = plt.subplots()

for pair in coords:
    line = LineString([pair[0:2], pair[2:4]])
    lines.append(line)
    rectangle=line.buffer(0.5)
    rectangles.append(rectangle)
    #plot_line(line)
    plot_polygon(rectangle)


reculos=np.empty((len(rectangles),0)).tolist()
for n,r1 in enumerate(rectangles):
    for m,r2 in enumerate(rectangles):
        if r1.intersects(r2):
            amax = max(r1.area, r2.area)
            apor=r1.intersection(r2).area/amax
            if apor>=0.9:
                reculos[n].append(m)

lol=[x for x in reculos if len(x) > 1]
    #plot_polygon(rectangle)
#fig,ax=plt.subplots()
    #plot_line(line,ax=ax)
#plt.show()
# #estado=np.loadtxt(path_cdd+'estado.txt').astype(int)
# #coords=np.array([x for x,y in zip(coords,estado) if y != 0])
# #good=np.array([0,1,3,4,5,6,8,9,12,13,14,15,16,17,19,20,21,22,23,24,30,32,35,37,38,41,46])
# #coords=coords[good,:]
# xx = np.hstack((coords[:, 0], coords[:, 2]))
# xmin = np.min(xx)
# xmax = np.max(xx)
# yy = np.hstack((coords[:, 1], coords[:, 3]))
# ymin = np.min(yy)
# ymax = np.max(yy)
# print(xmin,xmax,ymin,ymax)
#
# xmin=-5
# xmax=30
# ymin=-22
# ymax=25
# wide=10
# length=10
# cols = np.arange(xmin-wide, xmax + wide,wide)
# rows = np.arange(ymin-length, ymax + length, length)
#
#     #cols= np.arange(xmin-wide,xmax+wide,20)
#     #rows= np.arange(ymin-length, ymax + length, 10)
# cols,rows,rect_isec,weights, indexes, polygons=create_grid(coords,xmin,xmax,ymin,ymax)