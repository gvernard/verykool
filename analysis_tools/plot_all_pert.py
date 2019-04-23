import os.path
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import json
import decimal
import re
import math

from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import colorConverter
from astropy.io import fits


import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))




path = sys.argv[1]
run  = sys.argv[2]
if len(sys.argv) > 3:
    step = sys.argv[3]
    out_path  = path+run+'output/'+str(step)+'_'
else:
    step = ''
    out_path  = path+run+'output/'

mycmap_div    = 'Spectral'
mycmap_seq    = 'YlGnBu'


fig       = plt.figure(figsize=(12,16))
plt.rcParams.update({'axes.titlesize': 8})
plt.rcParams.update({'axes.labelsize': 8})
data_path = path+'data/'

if os.path.isfile(data_path+'source.fits'):
    hdulist = fits.open(data_path+'source.fits')
    srcwidth = float(hdulist[0].header['size']) # this is the total width of he fits image of the source
    srcrange = srcwidth/2.0
    hdulist.close()
else:
    srcrange = 0.5



#dum = path.split('/')
#f   = open(path+dum[-2]+'.json','r')
f   = open(path+run+'vkl_input.json','r')
input_str = f.read()
#input_str = re.sub(r'\\\n', '', input_str)
#input_str = re.sub(r'//.*\n', '', input_str)
#pars      = json.loads(input_str)
input_str = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str = re.sub(re.compile("//.*?\n" ),"",input_str)
pars      = json.loads(input_str)


#with open(path+dum[-2]+'.json') as json_input:
#    pars = json.load(json_input)

xmin = -pars['iplane']['siz_x']/2.
xmax =  pars['iplane']['siz_x']/2.
ymin = -pars['iplane']['siz_y']/2.
ymax =  pars['iplane']['siz_y']/2.

xticks = np.arange(math.ceil( -pars['iplane']['siz_x']/2.),math.floor( pars['iplane']['siz_x']/2.),1)
xticks = np.append(xticks,math.floor( pars['iplane']['siz_x']/2.))
yticks = np.arange(math.ceil( -pars['iplane']['siz_y']/2.),math.floor( pars['iplane']['siz_y']/2.),1)
yticks = np.append(yticks,math.floor( pars['iplane']['siz_y']/2.))
srcxticks = np.arange(-srcrange,srcrange,2*srcrange/4.0)
srcxticks = np.append(srcxticks,srcrange)
srcyticks = np.arange(-srcrange,srcrange,2*srcrange/4.0)
srcyticks = np.append(srcyticks,srcrange)






# Real/Mock image data and model. They need to share the same color scale
#########################################################################################################################

# Read data
im_data  = fits.getdata(data_path+'image.fits',ext=0)
im_data  = im_data[::-1,:]
# Read model
im_model = fits.getdata(out_path+'pert_model.fits',ext=0)
im_model = im_model[::-1,:]
# Read mask
if pars["maskpath"] == "0":
    mask = np.ones(im_model.shape)
else:
    mask = fits.getdata(data_path+'mask.fits',ext=0)

# Set color scale
limit = max([abs(np.amax(im_data)),abs(np.amin(im_data)),abs(np.amax(im_model)),abs(np.amin(im_model))])


# Plot data
data = fig.add_subplot(4,3,1)
data.set_title('DATA')
data.set_xlabel('arcsec')
data.set_ylabel('arcsec')
#im = data.imshow(im_data,interpolation='none',cmap=mycmap_seq,extent=[xmin,xmax,ymin,ymax],vmin=d_min,vmax=d_max,norm=MidpointNormalize(midpoint=0.0))
im = data.imshow(im_data,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
data.set_xticks(xticks)
data.set_yticks(yticks)
# Plot colorbar
divider = make_axes_locatable(data)
cax = divider.append_axes("right",size="5%",pad=0.05)
#fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
cb = fig.colorbar(im,cax=cax,format="%5.2f")
cb.set_ticks(matplotlib.ticker.MaxNLocator(nbins=10),update_ticks=True)
# Plot mask contour
if pars["maskpath"] != "0":
    data.contour(mask,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
    data.imshow(np.flipud(np.ma.masked_where(mask>0,mask)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])


# Plot model
model = fig.add_subplot(4,3,2)
model.set_title('MODEL')
model.set_xlabel('arcsec')
model.set_ylabel('arcsec')
#im = data.imshow(im_data,interpolation='none',cmap=mycmap_seq,extent=[xmin,xmax,ymin,ymax],vmin=d_min,vmax=d_max,norm=MidpointNormalize(midpoint=0.0))
im = model.imshow(im_model,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
model.set_xticks(xticks)
model.set_yticks(yticks)
# Plot colorbar
divider = make_axes_locatable(model)
cax = divider.append_axes("right",size="5%",pad=0.05)
#fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
fig.colorbar(im,cax=cax,format="%5.2f")
# Plot mask contour
if pars["maskpath"] != "0":
    model.contour(mask,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
    model.imshow(np.flipud(np.ma.masked_where(mask>0,mask)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])





# Residual image between data and model
#########################################################################################################################
# Read residuals
im_res = fits.getdata(out_path+'pert_residual.fits',ext=0)
im_res = im_res[::-1,:]
if pars["maskpath"] != "0":
    im_res = im_res*np.flipud(np.array(mask))

#tmp = np.flipud(im_res)
#hdu = fits.PrimaryHDU(tmp)
#hdu.writeto('new.fits',clobber=True)

# Set color scale
limit = max([abs(np.amax(im_res)),abs(np.amin(im_res))])

# Plot residuals
res = fig.add_subplot(4,3,3)
res.set_title('RESIDUAL')
res.set_xlabel('arcsec')
res.set_ylabel('arcsec')
im = res.imshow(im_res,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
res.set_xticks(xticks)
res.set_yticks(yticks)
# Plot colorbar
divider = make_axes_locatable(res)
cax = divider.append_axes("right",size="5%",pad=0.05)
#fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
#cax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
fig.colorbar(im,cax=cax,format="%5.2f")
if pars["maskpath"] != "0":
    res.contour(mask,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
    res.imshow(np.flipud(np.ma.masked_where(mask>0,mask)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])


# True source (in case of mock data) and reconstructed adaptive source. They need to share the same color scale
#########################################################################################################################
# Read true source if it exists
if os.path.isfile(data_path+'source.fits'):
    im_src = fits.getdata(data_path+'source.fits',ext=0)
    im_src = im_src[::-1,:]
    max_true = max([abs(np.amax(im_src)),abs(np.amin(im_src))])
else:
    max_true = 0

# Read reconstructed adaptive source
f        = open(out_path+'pert_source_voronoi.dat')
content  = [x.strip('\n') for x in f.readlines()]
zvalues  = []
polygons = []

dumcmap = matplotlib.cm.get_cmap(mycmap_seq)
for c in content:
    cell = np.fromstring(c,sep=' ')

    length = cell.shape[0]
    points = []
    for i in range(1,length,2):
        points.append([cell[i],cell[i+1]])

    polygon = Polygon(points,True)
    polygons.append(polygon)

    zvalues.append( cell[0] )



# Set color scale
limit = max([max_true,abs(max(zvalues)),abs(min(zvalues))])


# Plot real source if it exists
if os.path.isfile(data_path+'source.fits'):
    src = fig.add_subplot(4,3,4)
    src.set_title('TRUE SOURCE')
    src.set_xlabel('arcsec')
    src.set_ylabel('arcsec')
    im = src.imshow(im_src,interpolation='none',cmap=mycmap_div,extent=[-srcrange,srcrange,srcrange,-srcrange],vmin=-limit,vmax=limit)
    src.set_xticks(srcxticks)
    src.set_yticks(srcyticks)
    # Plot colorbar
    divider = make_axes_locatable(src)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    #fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
    fig.colorbar(im,cax=cax,format="%5.2f")
else:
    src = fig.add_subplot(4,3,4)
    src.axis('off')



# Plot reconstructed adaptive source
col = PatchCollection(polygons,alpha=1.0,cmap=mycmap_div)
col.set_clim(-limit,limit)
col.set_edgecolor('face')
col.set_array(np.array(zvalues))

vkl_src = fig.add_subplot(4,3,5)
vkl_src.add_collection(col)
vkl_src.set_ylim(-srcrange,srcrange)
vkl_src.set_xlim(-srcrange,srcrange)
vkl_src.title.set_text('RECONSTRUCTION')
vkl_src.set_xlabel('arcsec')
vkl_src.set_ylabel('arcsec')
vkl_src.set_xticks(srcxticks)
vkl_src.set_yticks(srcyticks)
plt.gca().set_aspect('equal',adjustable='box')
vkl_src.plot()
# Plot colorbar
divider = make_axes_locatable(vkl_src)
cax = divider.append_axes("right",size="5%",pad=0.05)
#plt.colorbar(col,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
plt.colorbar(col,cax=cax,format="%5.2f")







# Errors on reconstructed adaptive source
#########################################################################################################################
f        = open(out_path+'pert_source_voronoi_errors.dat')
content  = [x.strip('\n') for x in f.readlines()]
zvalues  = []
polygons = []

dumcmap = matplotlib.cm.get_cmap(mycmap_div)
for c in content:
    cell = np.fromstring(c,sep=' ')

    length = cell.shape[0]
    points = []
    for i in range(1,length,2):
        points.append([cell[i],cell[i+1]])

    polygon = Polygon(points,True)
    polygons.append(polygon)

    zvalues.append( cell[0] )

col = PatchCollection(polygons,alpha=1.0,cmap=mycmap_seq)
col.set_edgecolor('face')
col.set_array(np.array(zvalues))

vkl_err = fig.add_subplot(4,3,6)

vkl_err.add_collection(col)
vkl_err.set_ylim(-srcrange,srcrange)
vkl_err.set_xlim(-srcrange,srcrange)
vkl_err.title.set_text('ERRORS')
vkl_err.set_xlabel('arcsec')
vkl_err.set_ylabel('arcsec')
vkl_err.set_xticks(srcxticks)
vkl_err.set_yticks(srcyticks)
plt.gca().set_aspect('equal',adjustable='box')
vkl_err.plot()
# Plot colorbar
divider = make_axes_locatable(vkl_err)
cax = divider.append_axes("right",size="5%",pad=0.05)
#plt.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
plt.colorbar(col,cax=cax,format="%6.4f")







# True perturbations (in case of mock data) and reconstructed perturbations. They need to share the same color scale
#########################################################################################################################
# Read true perturbations file if it exists
if os.path.isfile(data_path+'perturbations.fits'):
    im_pert_true = fits.getdata(data_path+'perturbations.fits',ext=0)
    im_pert_true = im_pert_true[::-1,:]
#    im_pert_true /= np.power(6.05/121,2)
    max_pert_true = np.amax([abs(np.amax(im_pert_true)),abs(np.amin(im_pert_true))])
else:
    max_pert_true = 0

# Read reconstructed perturbations
if os.path.isfile(out_path+'pert_added_dpsi.fits'):
    im_pert = fits.getdata(out_path+'pert_added_dpsi.fits',ext=0)
    pert_title = "ADDED RECONSTRUCTION"
else:
    im_pert = fits.getdata(out_path+'pert_dpsi.fits',ext=0)
    pert_title = "RECONSTRUCTION"
im_pert = im_pert[::-1,:]
#im_pert /= np.power(6.05/40,2)

# Read mask
if pars["maskpath"] == "0":
    mask_dpsi = np.ones(im_pert.shape)
else:
    mask_dpsi = fits.getdata(out_path+'pert_dpsi_mask.fits',ext=0)


# Set color scale
#limit = np.amax([max_pert_true,abs(np.amax(im_pert*np.flipud(mask_dpsi))),abs(np.amin(im_pert*np.flipud(mask_dpsi)))])
limit = max_pert_true
#np.ma.masked_where(mask>0,mask)


# Plot original perturbations if they exist
if os.path.isfile(data_path+'perturbations.fits'):
    pert_true = fig.add_subplot(4,3,7)
    pert_true.set_title('TRUE PERTURBATIONS')
    pert_true.set_xlabel('arcsec')
    pert_true.set_ylabel('arcsec')
    im = pert_true.imshow(im_pert_true,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
    pert_true.set_xticks(xticks)
    pert_true.set_yticks(yticks)
    # Plot colorbar
    divider = make_axes_locatable(pert_true)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    #fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
    fig.colorbar(im,cax=cax,format="%5.2f")
    if pars["maskpath"] != "0":
        pert_true.contour(mask_dpsi,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
        pert_true.imshow(np.flipud(np.ma.masked_where(mask_dpsi>0,mask_dpsi)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])
else:
    pert_true = fig.add_subplot(4,3,7)
    pert_true.axis('off')



# Plot reconstructed perturbations
pert = fig.add_subplot(4,3,8)
pert.set_title(pert_title)
pert.set_xlabel('arcsec')
pert.set_ylabel('arcsec')
im = pert.imshow(im_pert,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
pert.set_xticks(xticks)
pert.set_yticks(yticks)
# Plot colorbar
divider = make_axes_locatable(pert)
cax = divider.append_axes("right",size="5%",pad=0.05)
#fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
#cax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
fig.colorbar(im,cax=cax,format="%5.2f")
if pars["maskpath"] != "0":
    pert.contour(mask_dpsi,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
    pert.imshow(np.flipud(np.ma.masked_where(mask_dpsi>0,mask_dpsi)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])


dum = fig.add_subplot(4,3,9)
dum.axis('off')


# True convergence (in case of mock data) and reconstructed convergence. They need to share the same color scale
#########################################################################################################################
# Read true convergence file if it exists
if os.path.isfile(data_path+'convergence.fits'):
    im_conv_true = fits.getdata(data_path+'convergence.fits',ext=0)
    im_conv_true = im_conv_true[::-1,:]
    max_conv_true = np.amax([abs(np.amax(im_conv_true)),abs(np.amin(im_conv_true))])
else:
    max_conv_true = 0

# Read reconstructed convergence
if os.path.isfile(out_path+'pert_added_convergence.fits'):
    im_conv = fits.getdata(out_path+'pert_added_convergence.fits',ext=0)
    conv_title = "ADDED CONVERGENCE"
else:
    im_conv = fits.getdata(out_path+'pert_convergence.fits',ext=0)
    conv_title = "CONVERGENCE"
im_conv = im_conv[::-1,:]

# Set color scale
limit = np.amax([max_conv_true,abs(np.amax(im_conv)),abs(np.amin(im_conv))])




# Plot original perturbations if they exists
if os.path.isfile(data_path+'convergence.fits'):
    conv_true = fig.add_subplot(4,3,10)
    conv_true.set_title('TRUE CONVERGENCE')
    conv_true.set_xlabel('arcsec')
    conv_true.set_ylabel('arcsec')
    im = conv_true.imshow(im_conv_true,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
    conv_true.set_xticks(xticks)
    conv_true.set_yticks(yticks)
    # Plot colorbar
    divider = make_axes_locatable(conv_true)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    #fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
    fig.colorbar(im,cax=cax,format="%5.2f")
    if pars["maskpath"] != "0":
        conv_true.contour(mask_dpsi,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
        conv_true.imshow(np.flipud(np.ma.masked_where(mask_dpsi>0,mask_dpsi)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])
else:
    conv_true = fig.add_subplot(4,3,10)
    conv_true.axis('off')



# Plot reconstructed convergence
conv = fig.add_subplot(4,3,11)
conv.set_title(conv_title)
conv.set_xlabel('arcsec')
conv.set_ylabel('arcsec')
im = conv.imshow(im_conv,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
conv.set_xticks(xticks)
conv.set_yticks(yticks)
# Plot colorbar
divider = make_axes_locatable(conv)
cax = divider.append_axes("right",size="5%",pad=0.05)
#fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
#cax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
fig.colorbar(im,cax=cax,format="%5.2f")
if pars["maskpath"] != "0":
    conv.contour(mask_dpsi,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
    conv.imshow(np.flipud(np.ma.masked_where(mask_dpsi>0,mask_dpsi)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])



dum = fig.add_subplot(4,3,12)
dum.axis('off')






fig.suptitle('step: '+str(step),fontsize=10)
plt.tight_layout(rect=[0, 0.0, 1, 0.95])
#plt.tight_layout()
plt.savefig('all.png',bbox_inches='tight')
plt.savefig('all.pdf',bbox_inches='tight')
#plt.show()
