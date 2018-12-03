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
image = fits.getdata(fname,ext=0)
dum  = image[::-1,:]
#data = np.flipud(dum)
image=dum
header = fits.getheader(fname)
Lx = header['WIDTH']
Ly = header['HEIGHT']
#Lx = 4.05
#Ly = 4.05


if len(sys.argv) == 3:
        mask = fits.getdata(sys.argv[2],ext=0)
        data = np.multiply(np.flipud(mask),image)
else:
        data = image




steps = 50 #number of bins in which to determine the power spectra
pspecs = []  #The list to be filled with power spectra 
l_list = []  #The list to be filled with l values	

#FFT & absolute value squared
fouriertf = np.fft.fftshift(np.fft.fft2(data))
absval2 = fouriertf.real**2 + fouriertf.imag**2
#Dimensions of the Fourier transformed image
nx = np.shape (absval2) [0]
ny = np.shape (absval2) [1]

steplist = range(steps + 1)
#print len(steplist),steplist
#dum = np.logspace(1,np.log10(steps),steps)
#steplist = np.insert(dum,0,0)
#print len(steplist),steplist
pspeclist = [] #Array that will be filled with the powerspectrum value at every l in the grid
pspecerror = [] #Array that will be filled with the error of the powerspectrum value at every l in the grid


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
        if len(bin) > 0:
                pspeclist = np.append(pspeclist, np.mean(bin)) #Adding each bin value to pspeclist
                pspecerror = np.append(pspecerror, np.mean(bin)/np.sqrt(len(bin))) #adding the error: P(k)/Npixels per bin (sample variance)
                #pspecerror = np.append(pspecerror, np.std(bin)) #adding the error: P(k)/Npixels per bin (sample variance)
        else:
                pspeclist = np.append(pspeclist,0.0)
                pspecerror = np.append(pspecerror,0.0)
        #print len(bin)


#The frequency list with the k-values at the middle of each bin with unit (k/2*pi) arcsec inverse
l_list = np.linspace(lmax/(2.*steps),lmax - lmax/(2.*steps),steps)

#Filtering out zero elements
ind = np.nonzero(pspeclist)
bins = l_list[ind]
vals = pspeclist[ind]
errs = pspecerror[ind]



#--------write the power spectrum-------------------------
table=np.array([bins,vals,errs]).T
np.savetxt("ps.dat",table,fmt='%.8e',delimiter=' ',newline='\n')

#--------plot the power spectrum--------------------------
plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.01,10000)

#plt.plot(l_list, pspecs[0])
plt.errorbar(bins,vals,yerr=errs)
plt.savefig('ps.png',bbox_inches='tight')
#plt.show()
