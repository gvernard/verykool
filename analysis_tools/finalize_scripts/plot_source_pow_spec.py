import os
from os.path import isfile, join
import sys
from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def filterZero(x,y,e):
    ind = np.nonzero(y)
    bins = x[ind]
    vals = y[ind]
    errs = e[ind]
    return bins,vals,errs



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

# reconstructed source (need to create gridded source first)
os.system("./create_gridded_source "+path+" "+run+" 200 "+lmodel+" "+step)
os.system("python plot_pow_spec.py source_model.fits")
x,y,e = np.loadtxt("ps.dat",unpack=True)
x,y,e = filterZero(x,y,e)
plt.errorbar(x,y,yerr=e,label='reconstruction',color='blue')

# true source
if os.path.isfile(path+"data/source.fits"):
    os.system("python plot_pow_spec.py "+path+"data/source.fits")
    x,y,e = np.loadtxt("ps.dat",unpack=True)
    x,y,e = filterZero(x,y,e)
    plt.errorbar(x,y,yerr=e,label='truth',color='black')

plt.legend()
plt.xlabel(r'$k [arcsec^{-1}]$',fontsize=17)
plt.ylabel(r'$P(k)$',fontsize=17)

plt.tight_layout()
plt.savefig('source_ps.pdf',bbox_inches='tight')
#plt.show()

os.remove("ps.dat")
os.remove("ps.pdf")
#os.remove("source_model.fits")
