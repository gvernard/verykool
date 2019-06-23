from astropy.io import fits
import sys
import numpy as np
import imageio

# Always using the first channel of a colored image.
# Does not work well with transparent png images.

infile  = sys.argv[1]
outfile = sys.argv[2]

# Read the image file
outdata = imageio.imread(infile)
if len(outdata.shape) == 3:
    outdata = np.flipud(outdata[:,:,0])
    
# Create Header and add keys
header = fits.Header()
header.set('WIDTH',4.0)
header['WIDTH'] = (4.0,'width of the image in arcsec')
header.set('HEIGHT',4.0)
header['HEIGHT'] = (4.0,'height of the image in arcsec')

# Write the output
newhdu = fits.PrimaryHDU(outdata,header=header)
newhdu.writeto(outfile,overwrite=True)

# Make sure to close the file
#hdulist.close()
