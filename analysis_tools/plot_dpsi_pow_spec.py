import os
from os.path import isfile, join
import sys
from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


path = sys.argv[1]
run  = sys.argv[2]
if len(sys.argv) == 4:
    lmodel = str(sys.argv[3])
    step = ''
    out_path = path + run + "output/" + lmodel
elif len(sys.argv) == 5:
    lmodel = str(sys.argv[3])
    step = str(sys.argv[4])
    out_path = path + run + "output/" + step + "_" + lmodel
else:
    print "Either 3 or 4 command line arguments required: path, run, lmodel, <step>"
    print len(sys.argv)," provided, exiting!!!"
    sys.exit()


fig = plt.figure(figsize=(10,6.2))
plt.xscale('log')
plt.yscale('log')

# reconstructed perturbations (need to create gridded source first)
os.system("python plot_pow_spec.py "+out_path+"_dpsi.fits "+path+run+"output/pert_dpsi_mask.fits")
#os.system("python plot_pow_spec.py "+out_path+"perturbations_vkl_source.fits")
x,y,e = np.loadtxt("ps.dat",unpack=True)
plt.errorbar(x,y,yerr=e,label='reconstruction',color='red')

# true perturbations
if os.path.isfile(path+"data/dpsi.fits"):
    os.system("python plot_pow_spec.py "+path+"data/dpsi.fits "+path+"data/mask.fits")
    #os.system("python plot_pow_spec.py "+path+"data/perturbations.fits")
    x,y,e = np.loadtxt("ps.dat",unpack=True)
    plt.errorbar(x,y,yerr=e,label='truth',color='black')

plt.legend(fontsize=17)
plt.xlabel(r'$k [arcsec^{-1}]$',fontsize=17)
plt.ylabel(r'$P(k)$',fontsize=17)
plt.xlim(0.4,5)

plt.tight_layout()
plt.savefig('dpsi_ps.pdf',bbox_inches='tight')
#plt.show()

os.remove("ps.dat")
os.remove("ps.pdf")

