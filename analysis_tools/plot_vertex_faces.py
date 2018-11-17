import numpy as np
import sys
import matplotlib.pyplot as plt





path = sys.argv[1]
case = sys.argv[2]
index = sys.argv[3]





xa,ya,xb,yb,xc,yc = np.loadtxt(path+case+"output/"+index+"_vertex_faces.dat",unpack=True)

x0 = xa[0]
y0 = ya[0]



for i in range(0,len(xa)):
    xa[i] = xa[i] - x0
    xb[i] = xb[i] - x0
    xc[i] = xc[i] - x0

    ya[i] = ya[i] - y0
    yb[i] = yb[i] - y0
    yc[i] = yc[i] - y0


plt.figure()
#for i in range(0,len(xa)):
#    plt.plot([xa[i],xb[i],xc[i]],[ya[i],yb[i],yc[i]])
#i=0
#plt.plot([xa[i],xb[i],xc[i],xa[i]],[ya[i],yb[i],yc[i],ya[i]])


plt.plot(np.append(xc,xc[0]),np.append(yc,yc[0]))
plt.scatter([xa,xb,xc],[ya,yb,yc])


plt.axhline(y=0,color='black')
plt.axvline(x=0,color='black')
plt.show()
