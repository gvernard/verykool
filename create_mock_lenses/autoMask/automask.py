#!/usr/bin/env python
import sys
import json
import re
import numpy as np

from scipy.signal import convolve
from scipy.signal import fftconvolve
from astropy.io import fits

#
# addMachine needs to have run first to produce the noise.dat file.
#

json_input = sys.argv[1]
f          = open(json_input,'r')
input_str  = f.read()
input_str  = re.sub(r'\\\n', '', input_str)
input_str  = re.sub(r'//.*\n', '', input_str)
myjson     = json.loads(input_str)

output     = myjson['output']
pix_x      = myjson['iplane']['pix_x']
pix_y      = myjson['iplane']['pix_y']
width      = myjson['iplane']['width']
height     = myjson['iplane']['height']
smear      = myjson['smear']
threshold  = myjson['threshold']
if 'noise_flag' in myjson:
    noise_flag = myjson['noise_flag']
else:
    noise_flag = 'none'

(nx,ny) = (pix_x,pix_y)
(dx,dy) = (float(width)/float(pix_x),float(height)/float(pix_y))

smear_pix = smear # smear is measured in pixel sizes


if( noise_flag == 'none' ):
    threshold_noise = threshold
elif( noise_flag == 'uniform' ):
    # uniform noise is just the sigma of the gaussian used to generated it, here we set the threshold to a few times this value
    sigma = np.loadtxt(output+'noise.dat')
    threshold_noise = threshold*sigma
else:
    threshold_noise = 0.0
#print threshold_noise


# create smearing kernel for convolution
b = 1.0/(2*smear_pix**2)
c = 1.0/(2*np.pi*smear_pix**2)
d = np.zeros( (nx,ny) )
for i in range(nx):
    for j in range(ny):
        d[i][j] = c*np.exp(-b*(((i-nx/2)*dx)**2+((j-ny/2)*dy)**2))
p = (d/(d.sum()))


# read lensed image
hdulist = fits.open(output+'image.fits')
image   = hdulist[0].data
hdulist.close()


# convolve with smearing kernel
out_image = fftconvolve(image,p,mode='same')

'''
fft_image = np.fft.fft2(image)
fft_p     = np.fft.fft2(p)
fft_conv  = np.empty_like(p,dtype=np.complex128)
rows = fft_image.shape[0]
cols = fft_image.shape[1]
for j in range(0,cols):
    for i in range(0,rows):
        real = fft_image[i][j].real*fft_p[i][j].real - fft_image[i][j].imag*fft_p[i][j].imag
        imag = fft_image[i][j].real*fft_p[i][j].imag + fft_image[i][j].imag*fft_p[i][j].real
        fft_conv[i][j] = real + imag*1j
out = np.fft.ifft2(fft_conv)
tmp = np.fft.ifftshift(out)
out_image = tmp.real
'''


#hdu   = fits.PrimaryHDU(out_image)
#hdu.writeto(output+'pyconvolved.fits',clobber=True)


# create and write mask
mask = np.zeros( (nx,ny) ) #output
for i in range(nx): #iterating over input array
    for j in range(ny):
        #mask[i][j] = 1
        if( out_image[i][j] > threshold_noise ):
            mask[i][j] = 1

hdu = fits.PrimaryHDU(mask)
hdu.writeto(output+'mask.fits',overwrite=True)


##------ Masking central 7x7 pixels -------
#dxx = (nx-7) / 2 #centre  x
#dyy = (ny-7) / 2 #centre  y
#for ii in xrange ( 7 ): #iterating over the central 7x7 pixels
#    for jj in xrange ( 7 ):
#        new [ ii + dxx ] [ jj + dyy ] = 1
        
        
