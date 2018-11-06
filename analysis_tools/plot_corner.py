from getdist import plots,chains,mcsamples
import numpy as np
import sys
import json
import os.path
import matplotlib
import re
from shutil import copyfile

matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman'],'size':25})
matplotlib.rc('text', usetex=True)

path = sys.argv[1]
run  = sys.argv[2]
rootdir = path+run+"output/"
if len(sys.argv) > 3:
    step = str(sys.argv[3])
else:
    step = ''
    

if os.path.isfile(rootdir+str(step)+"_corner.txt"):
    copyfile(rootdir+str(step)+"_corner.txt",rootdir+"plt_corner.txt")
elif os.path.isfile(rootdir+"corner.txt"):
    copyfile(rootdir+"corner.txt",rootdir+"plt_corner.txt")
else:
    sys.exit()




### get the active parameter names from verykool
###################################################################################################################
names,labels = np.genfromtxt(rootdir+'plt_corner.paramnames',dtype='str',unpack=True)



### Create the corner plot
###################################################################################################################
#g = plots.getSubplotPlotter(chain_dir=rootdir,analysis_settings={'ignore_rows': 110},subplot_size=4)
g = plots.getSubplotPlotter(chain_dir=rootdir,analysis_settings={"ignore_rows": 0.0,"smooth_scale_2D":0.9,"mult_bias_correction_order":1},subplot_size=4)
##g.settings.axes_fontsize = 25
##g.settings.lab_fontsize = 25
g.settings.rcSizes()
#g.triangle_plot(['plt_corner'],names,filled=True,smooth_scale_2D=3,smooth_scale_1D=3,mult_bias_correction_order=-1)
g.triangle_plot(['plt_corner'],names,filled=True)




### Plot the true values if available (mock data) or do nothing (if true data)
###################################################################################################################
#if os.path.isfile(rootdir+'data/params.json'):

#tmp = rootdir.split('/')
if os.path.isfile(path+run+'/vkl_input.json'):

    # get all the parameter names and values from the true data
    true_params = {}

    input_file = path+run+'/vkl_input.json'
    f          = open(input_file,'r')
    input_str  = f.read()
    input_str  = re.sub('//.*?\n|/\*.*?\*/','',input_str,flags=re.S)
    true_pars  = json.loads(input_str)
    

    for param in true_pars["physical"]["nlpars"]:
        true_params[param['nam']] = param['val']

    for lens in true_pars["lenses"]:
        for param in true_pars["lenses"][lens]["nlpars"]:
            true_params[lens+'_'+param['nam']] = param['val']

    for i in range(0,len(names)):
        if names[i] not in ["lambda","ampl","sdev"]:
            for ax in g.subplots[i:,i]:
                ax.axvline(true_params[names[i]],color='black',ls='solid',alpha=0.5)
            for ax in g.subplots[i,:i]:
                ax.axhline(true_params[names[i]],color='black',ls='solid',alpha=0.5)




### Export image
g.export('corner.pdf')
