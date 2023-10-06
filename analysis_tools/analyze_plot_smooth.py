import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json
import re
import math
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from astropy.io import fits



# Input args
path = sys.argv[1]
run  = sys.argv[2]

# Find model type
f   = open(os.path.join(path,run,'vkl_input.json'),'r')
input_str = f.read()
input_str = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str = re.sub(re.compile("//.*?\n" ),"",input_str)
inpars    = json.loads(input_str)
if inpars["parameter_model"] == "standard":
    lmodel = "smooth"
elif pars["parameter_model"] == "perturbations_standard":
    lmodel = "perts"
elif pars["parameter_model"] == "both":
    lmodel = "both"
else:
    print("Unknown parameter model: ",pars["parameter_model"])
    sys.exit()
    
# Make 'analysis' directory if not there
outpath = os.path.join(path,run,'analysis/')
if not os.path.exists(outpath):
    os.mkdir(outpath)

# Script specific paths and variables
if len(sys.argv) > 3:
    step = sys.argv[3]
    inpath  = os.path.join(path,run,'output',str(step)+'_')
    out_img_path = os.path.join(path,run,'analysis',str(step)+'_')
else:
    step = ''
    inpath = os.path.join(path,run,'output/')
    out_img_path = os.path.join(path,run,'analysis/')
if lmodel != "smooth":
    print("This routine can only print 'smooth' model results, this model is: ",lmodel)
    sys.exit()
data_path = os.path.join(path,'data/')
noise_file = inpars["noisepath"]

    
mycmap_div    = 'Spectral'
mycmap_seq    = 'YlGnBu'
fig       = plt.figure(figsize=(10,7))

if os.path.isfile(data_path+'source.fits'):
    hdulist   = fits.open(data_path+'source.fits')
    srcwidth  = float(hdulist[0].header['width']) # this is the total width of he fits image of the source
    srcheight = float(hdulist[0].header['height']) # this is the total height of he fits image of the source
    srcrange  = max(srcwidth,srcheight)/2.0
    hdulist.close()
else:
    srcrange = 1.2

#srcrange = 2.5




#with open(path+dum[-2]+'.json') as json_input:
#    pars = json.load(json_input)

xmin = -inpars['iplane']['siz_x']/2.
xmax =  inpars['iplane']['siz_x']/2.
ymin = -inpars['iplane']['siz_y']/2.
ymax =  inpars['iplane']['siz_y']/2.

xticks = np.arange(math.ceil( -inpars['iplane']['siz_x']/2.),math.floor( inpars['iplane']['siz_x']/2.),1)
xticks = np.append(xticks,math.floor( inpars['iplane']['siz_x']/2.))
yticks = np.arange(math.ceil( -inpars['iplane']['siz_y']/2.),math.floor( inpars['iplane']['siz_y']/2.),1)
yticks = np.append(yticks,math.floor( inpars['iplane']['siz_y']/2.))
srcxticks = np.arange(-srcrange,srcrange,2*srcrange/4.0)
srcxticks = np.append(srcxticks,srcrange)
srcyticks = np.arange(-srcrange,srcrange,2*srcrange/4.0)
srcyticks = np.append(srcyticks,srcrange)






# Real/Mock image data and model. They need to share the same color scale
#########################################################################################################################

# Read data
im_data  = fits.getdata(data_path+'data.fits',ext=0)
im_data  = im_data[::-1,:]
# Read model
im_model = fits.getdata(inpath+'smooth_model.fits',ext=0)
im_model = im_model[::-1,:]
# Read mask
if inpars["maskpath"] == "0":
    mask = np.ones(im_model.shape)
else:
    mask = fits.getdata(data_path+'mask.fits',ext=0)
    im_data = im_data*np.flipud(np.array(mask))
    im_model = im_model*np.flipud(np.array(mask))

# Set color scale
limit = max([abs(np.amax(im_data)),abs(np.amin(im_data)),abs(np.amax(im_model)),abs(np.amin(im_model))])


# Plot data
data = fig.add_subplot(231)
data.set_title('DATA')
data.set_xlabel('arcsec')
data.set_ylabel('arcsec')
im = data.imshow(im_data,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
data.set_xticks(xticks)
data.set_yticks(yticks)
# Plot colorbar
divider = make_axes_locatable(data)
cax = divider.append_axes("right",size="5%",pad=0.05)
#fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
fig.colorbar(im,cax=cax,format="%5.2f")
# Plot mask contour
if inpars["maskpath"] != "0":
    data.contour(mask,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
    data.imshow(np.flipud(np.ma.masked_where(mask>0,mask)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])


# Plot model
model = fig.add_subplot(232)
model.set_title('MODEL')
model.set_xlabel('arcsec')
model.set_ylabel('arcsec')
im = model.imshow(im_model,interpolation='none',cmap=mycmap_div,extent=[xmin,xmax,ymin,ymax],vmin=-limit,vmax=limit)
model.set_xticks(xticks)
model.set_yticks(yticks)
# Plot colorbar
divider = make_axes_locatable(model)
cax = divider.append_axes("right",size="5%",pad=0.05)
#fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
fig.colorbar(im,cax=cax,format="%5.2f")
# Plot mask contour
if inpars["maskpath"] != "0":
    model.contour(mask,levels=[0.5],extent=[xmin,xmax,ymin,ymax],colors='#001A47')
    model.imshow(np.flipud(np.ma.masked_where(mask>0,mask)),interpolation='none',cmap='Greys',extent=[xmin,xmax,ymin,ymax])





# Residual image between data and model
#########################################################################################################################
# Read residuals
im_res = fits.getdata(inpath+'smooth_residual.fits',ext=0)
im_res = im_res[::-1,:]
if inpars["maskpath"] != "0":
    im_res = im_res*np.flipud(np.array(mask))

noise_map = np.flipud(fits.getdata(os.path.join(path,noise_file),ext=0)) # This is the sigma map
#noise_map = np.sqrt(noise_map)
im_res = np.divide(im_res,noise_map)
    
#tmp = np.flipud(im_res)
#hdu = fits.PrimaryHDU(tmp)
#hdu.writeto('new.fits',clobber=True)

# Set color scale
limit = max([abs(np.amax(im_res)),abs(np.amin(im_res))])

# Plot residuals
res = fig.add_subplot(233)
#res.set_title('RESIDUAL')
res.set_title('norm RESIDUAL')
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
if inpars["maskpath"] != "0":
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
f        = open(inpath+'smooth_source_voronoi.dat')
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

    polygon = Polygon(points,closed=True)
    polygons.append(polygon)

    zvalues.append( cell[0] )



# Set color scale
#limit = max([max_true,abs(max(zvalues)),abs(min(zvalues))])
limit = 0.15



# Plot real source if it exists
if os.path.isfile(data_path+'source.fits'):
    src = fig.add_subplot(234)
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
    src = fig.add_subplot(234)
    src.axis('off')



# Plot reconstructed adaptive source
col = PatchCollection(polygons,alpha=1.0,cmap=mycmap_div)
col.set_clim(-limit,limit)
col.set_edgecolor('face')
#col.set_edgecolor('black')
col.set_array(np.array(zvalues))

vkl_src = fig.add_subplot(235)
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
f        = open(inpath+'smooth_source_voronoi_errors.dat')
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

    polygon = Polygon(points,closed=True)
    polygons.append(polygon)

    zvalues.append( cell[0] )

col = PatchCollection(polygons,alpha=1.0,cmap=mycmap_seq)
col.set_edgecolor('face')
col.set_array(np.array(zvalues))

vkl_err = fig.add_subplot(236)

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





plt.tight_layout()
#plt.savefig('all.pdf',bbox_inches='tight')
plt.savefig(out_img_path+"all.png",bbox_inches='tight')
#plt.show()
