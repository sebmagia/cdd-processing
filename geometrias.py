import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry
import descartes

circle = shapely.geometry.Point(5.0, 0.0).buffer(10.0)
clip_poly = shapely.geometry.Polygon([[-9.5, -2], [2, 2], [3, 4], [-1, 3]])
clipped_shape = circle.difference(clip_poly)

line = shapely.geometry.LineString([[-10, -5], [15, 5]])
line2 = shapely.geometry.LineString([[-10, -5], [-5, 0], [2, 3]])

print('Blue line intersects clipped shape:', line.intersects(clipped_shape))
print('Green line intersects clipped shape:', line2.intersects(clipped_shape))

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(*np.array(line).T, color='blue', linewidth=3, solid_capstyle='round')
ax.plot(*np.array(line2).T, color='green', linewidth=3, solid_capstyle='round')
ax.add_patch(descartes.PolygonPatch(clipped_shape, fc='blue', alpha=0.5))
ax.axis('equal')

plt.show()


rec1=shapely.geometry.box(0.0, 0.0, 2.0, 2.0, ccw=True)
rec2=shapely.geometry.box(2.0, 0.0, 4.0, 2.0, ccw=True)
rec3=shapely.geometry.box(0.0, 2.0, 2.0, 4.0, ccw=True)
rec4=shapely.geometry.box(2.0, 2.0, 4.0, 4.0, ccw=True)

line=shapely.geometry.LineString([(1, 1), (3, 3)])