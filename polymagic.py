#import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
from shapely.geometry import Point, Polygon, LineString, box
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely.figures import SIZE, BLUE, GRAY, RED, YELLOW, BLACK, set_limits
import matplotlib.pyplot as plt
#points = gpd.read_file('points.shp')
xmin=0
xmax=100
ymin=0
ymax=100
length = 10
wide = 5

cols = list(np.arange(xmin, xmax + wide, wide))
rows = list(np.arange(ymin, ymax + length, length))
fig,ax=plt.subplots()

polygons = []
for x in cols[:-1]:
    for y in rows[:-1]:
        polygons.append(Polygon([(x,y), (x+wide, y), (x+wide, y+length), (x, y+length)]))
        plot_polygon(polygons[-1],ax=ax)