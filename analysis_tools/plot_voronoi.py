import numpy as np
import sys

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import colorConverter


print matplotlib.__version__


path = sys.argv[1]


fig, ax = plt.subplots()


f        = open(path+'vkl_voronoi.dat')
#f        = open('voronoi.dat')
content  = [x.strip('\n') for x in f.readlines()]
cells    = []
colors   = []
polygons = []

#cmap = matplotlib.cm.get_cmap('Spectral')
cmap = matplotlib.cm.get_cmap('viridis')
#cmap = matplotlib.cm.get_cmap('Spectral')

for c in content:
    cell = np.fromstring(c,sep=' ')
    colors.append( cmap(cell[0]) )
    length = cell.shape[0]
    points = []
    for i in range(1,length,2):
        points.append([cell[i],cell[i+1]])

    polygon = Polygon(points,True)
    polygons.append(polygon)


#col = PatchCollection(polygons,cmap=plt.get_cmap('viridis'),alpha=1.0)
col = PatchCollection(polygons,alpha=1.0)
col.set_color(colors)
#col.set_edgecolor('none')

#colors = 100*np.random.rand(len(polygons))
#p = PatchCollection(polygons,facecolor='none')





ax.set_aspect('equal')
ax.add_collection(col)
ax = plt.Axes(fig, [0., 0., 1., 1.])
#plt.axis('off')



#myrange = 1.
myrange = 2.0
plt.ylim(-myrange,myrange)
plt.xlim(-myrange,myrange)


#plt.savefig('voronoi.png',frameon=False,bbox_inches='tight',pad_inches=-1)
#plt.savefig('voronoi.png',bbox_inches='tight')
plt.show()


