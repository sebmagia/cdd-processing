import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from shapely.geometry import Point, Polygon, LineString, box
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely.figures import SIZE, BLUE, GRAY, RED, YELLOW, BLACK, set_limits
from shapely.affinity import rotate, translate
import json
import math
import pyproj
import scipy.io
path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/INV/'
path_dmg='/home/doctor/Doctor/Magister/Tesis/databases/process_data/Damages/'
coords=np.loadtxt(path_cdd+'cdd_coordinates.dat')
#coords2=coords[[0,3,8,27,55,91]]
coords=np.delete(coords,[0,3,8,27,55,91],axis=0)
fig,ax=plt.subplots()
#ax.set_facecolor((144/255, 238/255, 144/255))
ax.set_facecolor((0.95,0.95,0.9))

#ax.axhspan(ymin = -30, ymax = 30, facecolor='palegreen', alpha=0.5)

for coord in coords:
     ax.plot([coord[0], coord[2]], [coord[1], coord[3]], 'r^', linewidth=1.5)

## coordenadas edificio
with open(path_dmg+'ejecta.geojson') as f:
    ejecta = json.load(f)
ejectas=[x['geometry']['coordinates'][0] for x in ejecta['features']]
#coords_eje=ejecta['features'][0]['geometry']['coordinates'][0]
with open(path_dmg+'cracking.geojson') as f:
    cracking = json.load(f)
## utm18
utm_zone=18
utm_proj = pyproj.Proj(proj='utm', zone=utm_zone, datum='WGS84')
## convertir coordenadas de ejecta y crackin de latlon a utm
shp_eje_xy=[]
eje_polys=[]
for shp_latlon in ejectas:
    shp_xy=[pyproj.transform(pyproj.Proj(proj='latlong', datum='WGS84'), utm_proj, ex, ey) for ex,ey in shp_latlon]
    shp_eje_xy.append([(x-shp_xy[0][0],y-shp_xy[0][1]) for x,y in shp_xy])
    eje_poly=Polygon(shp_eje_xy[-1])
    eje_polys.append(eje_poly)
    #plot_polygon(eje_poly, ax=ax, add_points=False, color='yellow')

#xxx=[(x-xx[0][0],y-xx[0][1]) for x,y in xx]
#eje_poly = Polygon(xxx)
dx=11.93
dy=25.65
SE = (-6,-5)
NE = (-6, -5 +dy)
NW = (-6 - dx, -5 +dy )
SW = (-6 -dx , -5)
riesco = Polygon([SE, NE, NW, SW])

plot_polygon(riesco, ax=ax, add_points = False, color='grey')
## poly 0
angle=30
centroid=eje_polys[0].centroid
rotated_poly=rotate(eje_polys[0],angle=angle,origin=centroid)
translated_poly=translate(rotated_poly,5.0,7.5)
#plot_polygon(eje_polys[0], ax=ax, add_points = False, color='yellow')
#plot_polygon(rotated_poly, ax=ax, add_points = False, color='red')
plot_polygon(translated_poly, ax=ax, add_points = False, color='#FFDB58')
## poly 1
angle=30
centroid=eje_polys[1].centroid
rotated_poly=rotate(eje_polys[1],angle=angle,origin=centroid)
translated_poly=translate(rotated_poly,20.0,-13.0)
#plot_polygon(eje_polys[1], ax=ax, add_points = False, color='yellow')
#plot_polygon(rotated_poly, ax=ax, add_points = False, color='red')
plot_polygon(translated_poly, ax=ax, add_points = False, color='#FFDB58')
## poly 2
angle=0
centroid=eje_polys[2].centroid
rotated_poly=rotate(eje_polys[2],angle=angle,origin=centroid)
translated_poly=translate(rotated_poly,-7.5,21.5)
#plot_polygon(eje_polys[2], ax=ax, add_points = False, color='yellow')
#plot_polygon(rotated_poly, ax=ax, add_points = False, color='red')
plot_polygon(translated_poly, ax=ax, add_points = False, color='#FFDB58')
## poly 3
angle=30
centroid=eje_polys[3].centroid
rotated_poly=rotate(eje_polys[3],angle=angle,origin=centroid)
translated_poly=translate(rotated_poly,-4.5,11.0)
#plot_polygon(eje_polys[3], ax=ax, add_points = False, color='yellow')
#plot_polygon(rotated_poly, ax=ax, add_points = False, color='red')
plot_polygon(translated_poly, ax=ax, add_points = False, color='#FFDB58')

## poly 4
angle=30
centroid=eje_polys[4].centroid
rotated_poly=rotate(eje_polys[4],angle=angle,origin=centroid)
translated_poly=translate(rotated_poly,20,16)
#plot_polygon(eje_polys[4], ax=ax, add_points = False, color='yellow')
#plot_polygon(rotated_poly, ax=ax, add_points = False, color='red')
plot_polygon(translated_poly, ax=ax, add_points = False, color='#FFDB58',label='pan')

#ax.plot(SE[0],SE[1],'ko')
#ax.plot(NE[0],NE[1],'ko')
#ax.plot(NW[0],NW[1],'ko')
#ax.plot(SW[0],SW[1],'ko')
SE = (27,-5)
NE = (27, -5 +dy)
NW = (27 + dx, -5 +dy )
SW = (27 + dx , -5)
errazuriz = Polygon([SE, NE, NW, SW])
plot_polygon(errazuriz, ax=ax, add_points = False, color='grey')
ax.text(-18.7,-6.0,'6.9cm',fontname='Arial',fontsize=12,fontweight='bold')
ax.text(-18.7+11.43,-6.0,'12.8cm',fontname='Arial',fontsize=12,fontweight='bold')
ax.text(-18.7+11.43,-6.0+25.65+1.1,'34.5cm',fontname='Arial',fontsize=12,fontweight='bold')
ax.text(-18.7,-6.0+25.65+1.1,'19.7cm',fontname='Arial',fontsize=12,fontweight='bold')


ax.text(-19+45.5,-6.0,'7.6cm',fontname='Arial',fontsize=12,fontweight='bold')
ax.text(-19+45.5+11.43,-6.0,'8.8cm',fontname='Arial',fontsize=12,fontweight='bold')
ax.text(-19+45.5+11.43,-6.0+25.65+1.1,'9.0cm',fontname='Arial',fontsize=12,fontweight='bold')
ax.text(-19+45.5,-6.0+25.65+1.1,'10.7cm',fontname='Arial',fontsize=12,fontweight='bold')


ax.text(-13.7,8.0, 'Riesco',fontname='Arial',fontsize=16,fontweight='bold')
ax.text(-13.7+44.5,8.0, 'Errázuriz',fontname='Arial',fontsize=16,fontweight='bold')
# Add the North direction arrow to the plot
north_x, north_y = 0.5, 0.9

arrow_start = (-18.5, 28.5)
arrow_end = (-17.5, 27.5)
dx = arrow_end[0] - arrow_start[0]
dy = arrow_end[1] - arrow_start[1]
arrow_angle = np.arctan2(dy, dx) * 180 / np.pi
arrow_props = dict(facecolor='black', arrowstyle="<-", lw=2.5)
ax.annotate('N', xy=arrow_end, xytext=arrow_start, arrowprops=arrow_props, fontsize=18, ha='center', va='center', rotation=45)

#spt = scipy.io.loadmat(path_dmg+'st_CPT_IIT.mat')

#ax.plot(SE[0],SE[1],'ko')
#ax.plot(NE[0],NE[1],'ko')
#ax.plot(NW[0],NW[1],'ko')
#ax.plot(SW[0],SW[1],'ko')
#for coord in coords2:
     #ax.plot([coord[0], coord[2]], [coord[1], coord[3]], 'k^', linewidth=1.5)
