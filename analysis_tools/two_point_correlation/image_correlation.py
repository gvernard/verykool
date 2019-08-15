import sys
import numpy as np
import math
import numpy.fft
from astropy.io import fits
import matplotlib.pyplot as plt


fname = sys.argv[1]
image = fits.getdata(fname,ext=0)
dum   = image[::-1,:]
image = np.flipud(dum)

image_width = 1.0 # arcsec
dpix = image_width/image.shape[1] # pixel size

# Fourier forward and inverse transforms should give the same image as the input
#fouriertf = np.fft.fft2(image)
#original  = np.fft.ifft2(fouriertf)



# Fourier transform image
fouriertf = np.fft.fft2(image,norm="ortho")

# Power spectrum (the square of the signal)
absval2 = fouriertf.real**2 + fouriertf.imag**2

# Correlation function (the inverse fourier transform of the power spectrum)
complex_corr = np.fft.fftshift(np.fft.ifft2(absval2,norm="ortho"))
corr = complex_corr.real

newhdu = fits.PrimaryHDU(corr)
newhdu.writeto('matrix_strue.fits',overwrite=True)


# Bin the 2D correlation function into radial bins
Nbins = 100
rmax = 2.0
rmin = 0.0
dr = (rmax-rmin)/Nbins
bins = np.arange(rmin,rmax,dr)
vals = np.zeros(len(bins))
counts = np.zeros(len(bins))

Nx = corr.shape[1]
Ny = corr.shape[0]
for i in range(0,Ny):
    for j in range(0,Nx):
        r = math.hypot((j-Nx/2.0)*dpix,(i-Ny/2.0)*dpix)
        if r < rmax and i != j:
            index = int(math.floor(r/dr))
            vals[index] += corr[i][j]
            counts[index] += 1

for i in range(0,Nbins):
    if counts[i] > 0:
        vals[i] /= counts[i]


np.savetxt("image_corr.dat",np.c_[bins,vals])
