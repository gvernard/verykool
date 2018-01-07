import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
import matplotlib.path as mplPath

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)






## make up data points
#np.random.seed(1234)
#points2 = np.random.rand(15, 2)



v    = np.loadtxt('voronoi.dat')
rows = v.shape[0]
points = np.empty([rows,2])
colors = np.empty([rows])
for i in range(0,rows):
    points[i,0] = v[i,0]
    points[i,1] = v[i,1]
    colors[i]   = v[i,2]

#print max(colors),min(colors)
print colors


# compute Voronoi tesselation
vor = Voronoi(points)


# plot
regions, vertices = voronoi_finite_polygons_2d(vor)

# colorize
#cmap=plt.cm.PuRd
cmap=plt.cm.gray
for i,region in enumerate(regions):
    polygon = vertices[region]
    bbPath = mplPath.Path(polygon)

    index = 0
    for j,point in enumerate(points):
        if bbPath.contains_point(point):
            index = j
            break

    plt.fill(*zip(*polygon), color=cmap(colors[index]), alpha=1.)



plt.plot(points[:,0], points[:,1],'ko',ms=0.4)
#plt.xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
#plt.ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)
plt.xlim(-0.5,0.5)
plt.ylim(-0.5,0.5)

plt.show()
