from getdist import plots,chains,mcsamples
import numpy as np
import sys
import json
import os.path
import matplotlib
from shutil import copyfile

matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman'],'size':25})
matplotlib.rc('text', usetex=True)

path = sys.argv[1]
run  = sys.argv[2]
rootdir = path+run+"output/"
if len(sys.argv) > 3:
    step = sys.argv[3]
    copyfile(rootdir+str(step)+"_corner.txt",rootdir+"plt_corner.txt")
else:
    copyfile(rootdir+"corner.txt",rootdir+"plt_corner.txt")



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
if os.path.isfile(rootdir+'/vkl_input.json'):

    # get all the parameter names and values from the true data
    true_params = {}

    with open(rootdir+'/vkl_input.json') as json_input:    
        true_pars = json.load(json_input)

    for param in true_pars["physical"]:
        true_params[param['nam']] = param['val']

    if (len(true_pars["lenses"]) == 1) :
        for param in true_pars["lenses"][0]["nlpars"]:
            true_params[param['nam']] = param['val']
    else:
        for i in range(0,len(true_pars["lenses"])):
            for param in true_pars["lenses"][i]["nlpars"]:
                true_params[param['nam']+'_'+str(i)] = param['val']


    # match the active parameters to their true values
    final_params = np.empty([0])
    final_values = np.empty([0])
    for name in used_pars:
        if name == 'lambda':
            final_params = np.append(final_params,name)
            final_values = np.append(final_values,0) 
        else:
            final_params = np.append(final_params,name)
            final_values = np.append(final_values,true_params[name])
    

    # add vertical lines at the location of the true values (exept "lambda")
    for i in range(0,len(final_params)):
        if final_params[i] != 'lambda':
            for ax in g.subplots[i:,i]:
                ax.axvline(final_values[i],color='black',ls='solid',alpha=0.5)
        
    # add horizontal lines at the location of the true values (exept "lambda")
    for i in range(0,len(final_params)):
        if final_params[i] != 'lambda':
            for ax in g.subplots[i,:i]:
                ax.axhline(final_values[i],color='black',ls='solid',alpha=0.5)





### Export image
g.export('corner.pdf')
