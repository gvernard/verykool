import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.spatial import Voronoi, voronoi_plot_2d



fig, ax = plt.subplots()


t    = np.loadtxt('triangles.dat')
rows = t.shape[0]
tris = []
for i in range(0,rows):
    polygon = Polygon([[t[i,0],t[i,1]],[t[i,2],t[i,3]],[t[i,4],t[i,5]]],True)
    tris.append(polygon)

p = PatchCollection(tris,cmap=matplotlib.cm.Set1,alpha=0.1,linewidth=0)
colors = 100*np.random.rand(len(tris))
p.set_array(np.array(colors))
p = PatchCollection(tris,facecolor='none')

ax.add_collection(p)


#v    = np.loadtxt('voronoi.dat')
#rows = v.shape[0]
#points = []
#for i in range(0,rows):
#    points.append([v[i,0],v[i,1]])

#vor = Voronoi(points)
##voronoi_plot_2d(vor)

#plt.plot(vor.vertices[:, 0], vor.vertices[:, 1], '*')




ax.set_aspect('equal')
myrange = 8
plt.ylim(-myrange,myrange)
plt.xlim(-myrange,myrange)


#plt.savefig('triangles.png')
plt.show()
