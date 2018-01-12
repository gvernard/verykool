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
from astropy.visualization import astropy_mpl_style

path = sys.argv[1]
run  = sys.argv[2]
if len(sys.argv) > 3:
    step = sys.argv[3]
    out_path  = path+run+'output/'+str(step)+'_'
else:
    out_path  = path+run+'output/'+str(step)+'_'

mycmap    = 'Spectral'
fig       = plt.figure(figsize=(12,8))
data_path = path+'data/'

hdulist = fits.open(data_path+'source.fits')
srcrange = float(hdulist[0].header['size'])
hdulist.close()
#srcrange = 0.4


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
srcxticks = np.arange(-srcrange,srcrange,1)
srcxticks = np.append(srcxticks,srcrange)
srcyticks = np.arange(-srcrange,srcrange,1)
srcyticks = np.append(srcyticks,srcrange)



# Real/Mock image data and model
im_data  = fits.getdata(data_path+'image.fits',ext=0)
im_data  = im_data[::-1,:]
im_model = fits.getdata(out_path+'vkl_image.fits',ext=0)
im_model = im_model[::-1,:]
mask     = fits.getdata(data_path+'mask.fits',ext=0)

d_max  = np.amax(im_data)
d_min  = np.amin(im_data)
m_max  = np.amax(im_model)
m_min  = np.amin(im_model)
dm_max = np.maximum(d_max,m_max)
dm_min = np.minimum(d_min,m_min)

data = fig.add_subplot(231)
data.set_title('DATA')
data.set_xlabel('arcsec')
data.set_ylabel('arcsec')
im = data.imshow(im_data,interpolation='nearest',cmap=mycmap,extent=[xmin,xmax,ymin,ymax],vmin=dm_min,vmax=dm_max)
data.set_xticks(xticks)
data.set_yticks(yticks)
divider = make_axes_locatable(data)
cax = divider.append_axes("right",size="5%",pad=0.05)
fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
data.contour(mask,levels=[0],extent=[xmin,xmax,ymin,ymax],colors='#001A47')

model = fig.add_subplot(232)
model.set_title('MODEL')
model.set_xlabel('arcsec')
model.set_ylabel('arcsec')
im = model.imshow(im_model,interpolation='nearest',cmap=mycmap,extent=[xmin,xmax,ymin,ymax],vmin=dm_min,vmax=dm_max)
model.set_xticks(xticks)
model.set_yticks(yticks)
divider = make_axes_locatable(model)
cax = divider.append_axes("right",size="5%",pad=0.05)
fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
model.contour(mask,levels=[0],extent=[xmin,xmax,ymin,ymax],colors='#001A47')


# Residual image between data and model
im_res = fits.getdata(out_path+'vkl_residual.fits',ext=0)
im_res = im_res[::-1,:]
res = fig.add_subplot(233)
res.set_title('RESIDUAL')
res.set_xlabel('arcsec')
res.set_ylabel('arcsec')
im = res.imshow(im_res,interpolation='nearest',cmap=mycmap,extent=[xmin,xmax,ymin,ymax])
res.set_xticks(xticks)
res.set_yticks(yticks)
divider = make_axes_locatable(res)
cax = divider.append_axes("right",size="5%",pad=0.05)
fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")



# True source (in case of mock data)
if os.path.isfile(data_path+'source.fits'):
    im_src = fits.getdata(data_path+'source.fits',ext=0)
    im_src = im_src[::-1,:]
    src = fig.add_subplot(234)
    src.set_title('TRUE SOURCE')
    src.set_xlabel('arcsec')
    src.set_ylabel('arcsec')
    im = src.imshow(im_src,interpolation='nearest',cmap=mycmap,extent=[-srcrange,srcrange,srcrange,-srcrange])
    src.set_xticks(srcxticks)
    src.set_yticks(srcyticks)
    divider = make_axes_locatable(src)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    fig.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")
else:
    src = fig.add_subplot(234)
    src.axis('off')



# Reconstructed adaptive source
f        = open(out_path+'vkl_voronoi.dat')
content  = [x.strip('\n') for x in f.readlines()]
cells    = []
colors   = []
polygons = []

dumcmap = matplotlib.cm.get_cmap(mycmap)
for c in content:
    cell = np.fromstring(c,sep=' ')
    colors.append( dumcmap(cell[0]) )
    length = cell.shape[0]
    points = []
    for i in range(1,length,2):
        points.append([cell[i],cell[i+1]])

    polygon = Polygon(points,True)
    polygons.append(polygon)


col = PatchCollection(polygons,alpha=1.0)
col.set_color(colors)


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
divider = make_axes_locatable(vkl_src)
cax = divider.append_axes("right",size="5%",pad=0.05)
plt.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")






# Errors on reconstructed adaptive source
f        = open(out_path+'vkl_voronoi_errors.dat')
content  = [x.strip('\n') for x in f.readlines()]
cells    = []
colors   = []
polygons = []

dumcmap = matplotlib.cm.get_cmap(mycmap)
for c in content:
    cell = np.fromstring(c,sep=' ')
    colors.append( dumcmap(cell[0]) )
    length = cell.shape[0]
    points = []
    for i in range(1,length,2):
        points.append([cell[i],cell[i+1]])

    polygon = Polygon(points,True)
    polygons.append(polygon)


col = PatchCollection(polygons,alpha=1.0)
col.set_color(colors)


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
divider = make_axes_locatable(vkl_err)
cax = divider.append_axes("right",size="5%",pad=0.05)
plt.colorbar(im,cax=cax,ticks=MultipleLocator(0.2),format="%4.1f")





plt.tight_layout()
plt.savefig('all.png',bbox_inches='tight')
#plt.show()
