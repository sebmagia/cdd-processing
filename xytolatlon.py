import pyproj
import numpy as np
import matplotlib.pyplot as plt
path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/INV/'

# Define the projection used by Google Maps
google_proj = pyproj.Proj(proj='merc', datum='WGS84')

# Define the lat,lon coordinates of the first point
lat, lon = -36.79101, -73.08175
x0,y0= 0.0, -12.0
# Define the (x,y) coordinates of the other points
#xy_points = [(100, 100), (200, 200), (300, 300)]
xy_points=np.loadtxt(path_cdd+'cdd_coordinates.dat')
fig,ax=plt.subplots()
for co in xy_points:
     ax.plot([co[0], co[2]], [co[1], co[3]], 'ro', linewidth=1.5)
xy_points=np.delete(xy_points,[8,55],axis=0)
xy_sp=np.split(xy_points,2,axis=1)
xy_points2=np.vstack(xy_sp)
xy_points3=np.unique(xy_points2,axis=0)
#xy_points3=np.delete(xy_points3,5,axis=0)

fig,ax=plt.subplots()
plt.plot(xy_points3[:,0],xy_points3[:,1],'go')

utm_zone = int((lon + 180) // 6) + 1

# Define the UTM projection for the zone
utm_proj = pyproj.Proj(proj='utm', zone=utm_zone, datum='WGS84')

# Convert the lat,lon coordinates to UTM
x, y = pyproj.transform(pyproj.Proj(proj='latlong', datum='WGS84'), utm_proj, lon, lat)
xypr=xy_points3-np.array([x0,y0])
XYPR=xypr + np.array([x,y])
# Print the UTM coordinates
print(x, y)
fig,ax=plt.subplots()
plt.plot(XYPR[:,0],XYPR[:,1],'go')

lonpr, latpr = pyproj.transform(utm_proj, pyproj.Proj(proj='latlong', datum='WGS84'), XYPR[:,0], XYPR[:,1])

fig,ax=plt.subplots()
plt.plot(lonpr,latpr,'go')
np.savetxt('/home/doctor/Doctor/Magister/Tesis/Papers/Ordenados/Liquefaction_Maule/Damages/coords_lonlat.txt',np.vstack((lonpr,latpr)).T)
# xy_sp=np.split(xy_points,2,axis=1)
# xy_points=np.vstack((xy_sp[0],xy_sp[1]))
# xy_points=np.unique(xy_points,axis=0)
# xy_points=np.delete(xy_points,5,axis=0)
# # Define the UTM zone for the point (based on its longitude)
# utm_zone = int((lon + 180) // 6) + 1
#
# # Define the UTM projection for the zone
# utm_proj = pyproj.Proj(proj='utm', zone=utm_zone, datum='WGS84')
#
# # Convert the lat,lon coordinates to UTM
# x, y = pyproj.transform(pyproj.Proj(proj='latlong', datum='WGS84'), utm_proj, lon, lat)
# xypr=xy_points-np.array([x0,y0])
# XYPR=xypr + np.array([x,y])
# # Print the UTM coordinates
# print(x, y)