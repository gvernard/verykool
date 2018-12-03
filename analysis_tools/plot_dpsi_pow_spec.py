import os
from os.path import isfile, join
import sys
from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt



path = sys.argv[1]
run  = sys.argv[2]
if len(sys.argv) > 3:
    step = sys.argv[3]
    out_path  = path+run+'output/'+str(step)+'_'
else:
    step = ''
    out_path  = path+run+'output/'




fig = plt.figure(figsize=(10,6.2))
plt.xscale('log')
plt.yscale('log')

# reconstructed perturbations (need to create gridded source first)
os.system("python plot_pow_spec.py "+out_path+"perturbations_vkl_source.fits "+path+run+"output/mask_source.fits")
#os.system("python plot_pow_spec.py "+out_path+"perturbations_vkl_source.fits")
x,y,e = np.loadtxt("ps.dat",unpack=True)
plt.errorbar(x,y,yerr=e,label='reconstruction',color='blue')

# true perturbations
if os.path.isfile(path+"data/perturbations.fits"):
    os.system("python plot_pow_spec.py "+path+"data/perturbations.fits "+path+"data/mask.fits")
    #os.system("python plot_pow_spec.py "+path+"data/perturbations.fits")
    x,y,e = np.loadtxt("ps.dat",unpack=True)
    plt.errorbar(x,y,yerr=e,label='truth',color='black')

plt.legend()
plt.xlabel(r'$k [arcsec^{-1}]$',fontsize=17)
plt.ylabel(r'$P(k)$',fontsize=17)

plt.tight_layout()
plt.savefig('pert_ps.pdf',bbox_inches='tight')
#plt.show()

os.remove("ps.dat")
os.remove("ps.png")

