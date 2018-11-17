import sys
import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

file1 = sys.argv[1]
file2 = sys.argv[2]

data1 = fits.getdata(file1,ext=0)
data1 = data1[::-1,:]

data2 = fits.getdata(file2,ext=0)
data2 = data2[::-1,:]


dum = np.subtract(data1,data2)
# Also, I need to invert the array on the y-axis to output the correct fits file
diff = np.flipud(dum)
limit = np.amax([np.amax(diff),np.abs(np.amin(diff))])


hdu = fits.PrimaryHDU(diff)
hdu.writeto('diff.fits',overwrite=True)



fig    = plt.figure(figsize=(12,8))
myplot = fig.add_subplot(111)
myplot.set_title('MODEL')
myplot.set_xlabel('arcsec')
myplot.set_ylabel('arcsec')
im = myplot.imshow(diff,interpolation='none',cmap='Spectral',vmin=-limit,vmax=limit)

plt.gca().set_aspect('equal',adjustable='box')
# Plot colorbar
divider = make_axes_locatable(myplot)
cax = divider.append_axes("right",size="5%",pad=0.05)
fig.colorbar(im,cax=cax,format="%6.4f")

plt.tight_layout()
plt.savefig('diff.pdf',bbox_inches='tight')
