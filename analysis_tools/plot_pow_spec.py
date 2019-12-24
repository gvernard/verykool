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
header = fits.getheader(fname)
Lx = header['WIDTH']
Ly = header['HEIGHT']
nx = header['NAXIS1']
ny = header['NAXIS1']

#print Lx,Ly,nx,ny

if len(sys.argv) == 3:
        mask = fits.getdata(sys.argv[2],ext=0)
        data = np.multiply(mask,image)
else:
        data = image

#newhdu = fits.PrimaryHDU(data)
#newhdu.writeto('res.fits',overwrite=True)


Nsteps = 10  # number of bins in which to determine the power spectra
pspecs = []  # The list to be filled with power spectra 

#FFT & absolute value squared
fouriertf = np.fft.fftshift(np.fft.fft2(data))
##fouriertf *= 1.0/((Lx*Ly)*(nx*ny))
absval2 = fouriertf.real**2 + fouriertf.imag**2


#Axes in Fourier space
kxlist = np.arange( (-nx/2.0)/Lx, (nx/2.0)/Lx, 1.0/Lx )
kylist = np.arange( (-ny/2.0)/Ly, (ny/2.0)/Ly, 1.0/Ly )
#kmax   = np.sqrt( np.max(np.abs(lxlist))**2.0 + np.max(np.abs(lylist))**2.0 )
kmax   = np.sqrt( pow((nx/2.0)/Lx,2.0) + pow((ny/2.0)/Ly,2.0) )
kmin   = np.sqrt( pow(1.0/Lx,2.0) + pow(1.0/Ly,2.0) )


steplist = np.zeros(Nsteps+1)
k_list = np.zeros(Nsteps)



# Equi-distant linear radial bins
dlin = (lmax-lmin)/Nsteps
steplist[0] = lmin
for i in range(1,Nsteps+1):
        steplist[i] = kmin + i*dlin
        k_list[i-1] = kmin + i*dlin - dlin/2.0
#l_list = np.linspace(lmax/(2.0*Nsteps),lmax - lmax/(2.0*Nsteps),Nsteps)


## Equi-distant logarithmic radial bins
#dlog   = (np.log10(lmax) - np.log10(lmin))/Nsteps
#logmin = np.log10(lmin)
#steplist[0] = lmin
#for i in range(1,Nsteps+1):
#        steplist[i] = pow(10,logmin + i*dlog)
#        l_list[i-1] = pow(10,logmin + i*dlog - dlog/2.0)
        
        
        
        
        
        
#print len(steplist),steplist
#dum = np.logspace(1,np.log10(Nsteps),steps)
#steplist = np.insert(dum,0,0)
#print len(steplist),steplist
pspeclist = [] #Array that will be filled with the powerspectrum value at every l in the grid
pspecerror = [] #Array that will be filled with the error of the powerspectrum value at every l in the grid

#Generating the azimuthal averaged power spectra
for step in range(Nsteps):
        bin = []
	
        for x in range(nx):
                for y in range(ny):
                        kx = kxlist[x]
                        ky = kylist[y]
                        k  = m.sqrt( kx**2 + ky**2 )
                        if steplist[step] < k <= steplist[step+1]:
                                bin = np.append(bin, absval2[y][x])
        
        N = len(bin)

        if N > 0:
                mean = np.mean(bin)
                std  = np.std(bin)
                pspeclist  = np.append(pspeclist,mean)                        # Adding each bin mean value to pspeclist
                #pspecerror    = np.append(pspecerror,mean/np.sqrt(N))           # Adding the error
                pspecerror = np.append(pspecerror,std)

        else:
                pspeclist  = np.append(pspeclist,0.0)
                pspecerror = np.append(pspecerror,0.0)
        #print len(bin)

bins = k_list
vals = pspeclist
errs = pspecerror



#--------write the power spectrum-------------------------
table=np.array([bins,vals,errs]).T
np.savetxt("ps.dat",table,fmt='%.8e',delimiter=' ',newline='\n')

#--------plot the power spectrum--------------------------
plt.figure()
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('k')
plt.ylabel('P(k)')

#plt.ylim(0.1,20)
#plt.xlim(0.6,20)
plt.ylim(0.1,2.5)
plt.xlim(0.0,16)


# Filtering out zero elements
ind = np.nonzero(pspeclist)
bins = k_list[ind]
vals = pspeclist[ind]
errs = pspecerror[ind]

#plt.plot(bins,vals)
plt.errorbar(bins,vals,yerr=[errs,errs])

plt.savefig('ps.pdf',bbox_inches='tight')
#plt.show()
