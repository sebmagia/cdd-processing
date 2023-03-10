import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, LineString, box
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely.figures import SIZE, BLUE, GRAY, RED, YELLOW, BLACK, set_limits
rec1=box(0.0, 0.0, 2.0, 2.0, ccw=True)
rec2=box(2.0, 0.0, 4.0, 2.0, ccw=True)
rec3=box(0.0, 2.0, 2.0, 4.0, ccw=True)
rec4=box(2.0, 2.0, 4.0, 4.0, ccw=True)

line1=LineString([(1, 1), (2.1, 2.7)])
line2=LineString([(2.3, 1.4), (1.8, 2.6)])
line3=LineString([(1.4, 3.6), (1.2, 2.5)])
line4=LineString([(1.16, 2.1), (3.4, 3.7)])


fig = plt.figure(1, figsize=SIZE, dpi=90)
#fig = plt.figure(1, dpi=90)

ax = fig.add_subplot(111)
rectangles=[rec1,rec2,rec3,rec4]
rect_isec=[[],[],[],[]]
lines=[line1,line2,line3,line4]
colors=[RED,RED,RED]
for rectangle in rectangles:
    plot_polygon(rectangle, ax=ax, color=YELLOW)
for n,line in enumerate(lines):
    plot_line(line, ax=ax, color=BLACK)
    for k,rectangle in enumerate(rectangles):
        bool_intersection=line.intersects(rectangle)
        intersection=line.intersection(rectangle)
        print(intersection)
        if bool_intersection and intersection.geom_type is 'LineString':
            plot_polygon(rectangle,ax=ax,color=RED)
            #inter=line.intersection(rectangle).coords.xy
            #isec_r=np.array([inter[0].tolist(),inter[1].tolist()])
            rect_isec[k].append(intersection.length)
            #p1=isec[0].tolist()
            #p2=isec[1].tolist()

            #rect_isec[k].append()
            print('rectanglex ' +str(k).zfill(2) + ' is intersected in ', intersection)
        if bool_intersection and intersection.geom_type is 'Point':
            print('rectangleb ' +str(k).zfill(2) + ' is intersected in ', intersection)
        else:
            print(intersection)
weights=[[],[],[],[]]
for n,rec in enumerate(rect_isec):
    if len(rec)==1:
        weights[n].append(1.0)
    else:
        weights[n]=[x/sum(rec) for x in rec]

from shapely.ops import polygonize

lines = [

    ((0, 0), (1, 1)),

    ((0, 0), (0, 1)),

    ((0, 1), (1, 1)),

    ((1, 1), (1, 0)),

    ((1, 0), (0, 0))

    ]
lines = [ ((0, 0), (1, 0)),((1, 0), (1, 1)),((1, 1), (0, 1)),((1, 1), (2, 1)),((2, 1), (2, 2)),((2, 2), (1, 2)) ]
XD=list(polygonize(lines))
fig = plt.figure(1, figsize=SIZE, dpi=90)
#fig = plt.figure(1, dpi=90)

ax = fig.add_subplot(111)
for poly in XD:
    plot_polygon(poly, ax=ax, color=BLACK)
#weights=[[],[],[],[]]
#for x in rect_isec:
 #   if len(x) ==1 and np.l:

#plot_polygon(rec1,ax=ax,color=YELLOW)
#plot_polygon(rec2,ax=ax,color=RED)
#plot_polygon(rec3,ax=ax,color=RED)
#plot_polygon(rec4,ax=ax,color=RED)
#plot_line(line,ax=ax,color=BLACK)

#print(line.intersection(rec1))
#print(line.intersection(rec2))
#print(line.intersection(rec3))
#print(line.intersection(rec4))
#plt.show()

# 1: valid polygon
# ax = fig.add_subplot(121)
#
# ext = [(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)]
# int = [(1, 0), (0.5, 0.5), (1, 1), (1.5, 0.5), (1, 0)][::-1]
# polygon = Polygon(ext, [int])
#
# plot_polygon(polygon, ax=ax, add_points=False, color=BLUE)
# plot_points(polygon, ax=ax, color=GRAY, alpha=0.7)
#
# ax.set_title('a) valid')
#
# set_limits(ax, -1, 3, -1, 3)
#
# #2: invalid self-touching ring
# ax = fig.add_subplot(122)
# ext = [(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)]
# int = [(1, 0), (0, 1), (0.5, 1.5), (1.5, 0.5), (1, 0)][::-1]
# polygon = Polygon(ext, [int])
#
# plot_polygon(polygon, ax=ax, add_points=False, color=RED)
# plot_points(polygon, ax=ax, color=GRAY, alpha=0.7)
#
# ax.set_title('b) invalid')
#
# set_limits(ax, -1, 3, -1, 3)
#
# plt.show()
