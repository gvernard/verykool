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



data  = fits.getdata('/net/argo/data/users/gvernard/RESULTS/VeryKooL_DATA/mysim6/data/mask.fits',ext=0)
data  = data[::-1,:]

CS = plt.contour(data,levels=[0])



plt.show()
