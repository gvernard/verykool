#################################################################################################
#Importing necessary packages
#################################################################################################

import sys
import numpy as np
import math as m
from numpy.fft import fft, fft2, fftn, fftfreq, ifft, fftshift
from astropy.io import fits
import matplotlib.pyplot as plt

#################################################################################################
#1D power spectra of the image
#################################################################################################


fname = sys.argv[1]
myfits = fits.getdata(fname,ext=0)
dum  = myfits[::-1,:]
#data = np.flipud(dum)
data=dum
header = fits.getheader(fname)
#Lx = header['NAXIS1']
#Ly = header['NAXIS1']
Lx = 4.05
Ly = 4.05



steps = 30 #number of bins in which to determine the power spectra
pspecs = []  #The list to be filled with power spectra 
l_list = []  #The list to be filled with l values	

#FFT & absolute value squared
fouriertf = np.fft.fftshift(np.fft.fft2(data))
absval2 = fouriertf.real**2 + fouriertf.imag**2
#Dimensions of the Fourier transformed image
nx = np.shape (absval2) [0]
ny = np.shape (absval2) [1]

steplist = range (steps + 1)
pspeclist = [] #Array that will be filled with the powerspectrum value at every l in the grid

#Axes in Fourier space
lxlist = np.arange ((-nx/2.)/Lx,(nx/2.)/Lx, 1./Lx)
lylist = np.arange ((-ny/2.)/Ly,(ny/2.)/Ly, 1./Ly)
lmax = np.sqrt( np.max(np.abs(lxlist))**2. + np.max(np.abs(lylist))**2. )

#Generating the azimuthal averaged power spectra
for step in range (steps):
        bin = []
	
        for x in range (nx):
                for y in range (ny):
                        lx = lxlist [x]
                        ly = lylist [y]
                        l = m.sqrt(lx**2. + ly**2.)
                        if steplist[step]*lmax/steps < l <= steplist [step + 1]*lmax/steps:
                                bin = np.append(bin, absval2[y][x])
        pspeclist = np.append(pspeclist, np.mean(bin)) #Adding each bin value to pspeclist
pspecs.append(pspeclist) #Adding the entire powerspectrum to the list
	
#The frequency list with the k-values at the middle of each bin with unit (k/2*pi) arcsec inverse
l_list = np.linspace(lmax/(2.*steps),lmax - lmax/(2.*steps),steps)


#--------write the power spectrum-------------------------
table=np.array([l_list,pspecs[0]]).T
np.savetxt("ps.dat",table,fmt='%.8e',delimiter=' ',newline='\n')

#--------plot the power spectrum--------------------------
plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.01,10000)

plt.plot(l_list, pspecs[0])
plt.savefig('ps.png',bbox_inches='tight')
#plt.show()
