#
# Code to produce GRF realizations
#
import numpy as np
import sys
import json
from astropy.io import fits

import grf_functions as grf


json_input = sys.argv[1]
f          = open(json_input,'r')
input_str  = f.read()
myjson     = json.loads(input_str)


#----- Input parameters -----------------------------------
var_dpsi = np.power(10,float(myjson['logvar']))    # variance
pow_dpsi = float(myjson['slope'])                  # power-law exponent
Lx_dpsi  = float(myjson['width'])                  # field of view (in arc-sec)
Ly_dpsi  = float(myjson['height'])                 # field of view (in arc-sec)
nx_dpsi  = int(myjson['Nwidth'])                   # number of pixels
nx_dpsi  = int(myjson['Nheight'])                  # number of pixels
output   = myjson['output']
if 'seed' in myjson:
    np.random.seed(int(myjson['seed']))

L_dpsi = Lx_dpsi
n_dpsi = nx_dpsi




#------ Simulate GRF image --------------------------------
nx = n_dpsi
ny = n_dpsi
Lx = L_dpsi
Ly = L_dpsi

resolution = float(Lx)/float(nx)
X = np.linspace(-Lx/2.0,Lx/2.0,nx)
Y = np.linspace(-Ly/2.0,Ly/2.0,ny)
x,y = np.meshgrid(X,Y,sparse='True')

grpar      = np.asarray([nx,ny,Lx,Ly,var_dpsi,pow_dpsi])
gauss_rand = grf.gauss_rand_2d(grpar)


hdr = fits.Header()
hdr.set('WIDTH',Lx)
hdr['WIDTH'] = (Lx,'width of the image in arcsec')
hdr.set('HEIGHT',Ly)
hdr['HEIGHT'] = (Ly,'height of the image in arcsec')
hdu = fits.PrimaryHDU(gauss_rand,header=hdr)
hdu.writeto(output+'dpsi.fits',overwrite=True)
