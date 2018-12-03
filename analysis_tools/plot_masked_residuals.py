# Import the astropy fits tools
from astropy.io import fits
import sys
import numpy as np

# arg 1: full path to the residuals .fits file
# arg 2: full path to the mask .fits file

header = fits.getheader(sys.argv[1],ext=0)
data   = fits.getdata(sys.argv[1],ext=0)
mask   = fits.getdata(sys.argv[2],ext=0)

newdata = np.multiply(mask,data)

newhdu  = fits.PrimaryHDU(newdata,header=header)
newhdu.writeto('masked_residual.fits',clobber=True)
