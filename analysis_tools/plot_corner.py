from getdist import plots,chains,MCSamples,loadMCSamples
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
if pars["parameter_model"] == "standard":
    lmodel = "smooth"
elif pars["parameter_model"] == "perturbations_standard":
    lmodel = "perts"
elif pars["parameter_model"] == "both":
    lmodel = "both"
else:
    print("Unknown parameter model: ",pars["parameter_model"])
    sys.exit()


if len(sys.argv) > 3:
    step = sys.argv[3]
    out_path = path + run + "output/" + str(step) + "_" + lmodel
else:
    step = ''
    out_path = path + run + "output/" + lmodel



if os.path.isfile(out_path+"_postdist.txt"):
    copyfile(out_path+"_postdist.txt","plt_corner.txt")
    copyfile(path+run+"output/"+lmodel+"_postdist.paramnames","plt_corner.paramnames")
    copyfile(path+run+"output/"+lmodel+"_postdist.ranges","plt_corner.ranges")
else:
    print("Postdist not found, exiting")
    sys.exit()



### get the active parameter names from verykool
###################################################################################################################
names,labels = np.genfromtxt('plt_corner.paramnames',dtype='str',unpack=True)
labels = [w.replace('_', '-') for w in labels]
#labels = [w.strip() for w in labels]
#names  = [w.strip() for w in names]

dum,mins,maxs = np.genfromtxt('plt_corner.ranges',unpack=True)
mapping = {}
for i in range(0,len(names)):
    mapping[names[i]] = (mins[i],maxs[i])

ola = np.loadtxt("./plt_corner.txt")
pars = ola[:,2:]
#mypost = ola[:,0]
mypost = ola[:,1]
#for i in range(0,len(mypost)):
#    mypost[i] = i


'''
tmp = np.copy(pars[:,2])
pars[:,2] = pars[:,1]
pars[:,1] = pars[:,0]
pars[:,0] = tmp
tmp = names[2]
names[2] = names[1]
names[1] = names[0]
names[0] = tmp
tmp = labels[2]
labels[2] = labels[1]
labels[1] = labels[0]
labels[0] = tmp

tmp = np.copy(pars[:,0])
pars[:,0] = np.log10(tmp)
'''


# ### Change to log scale
# ###################################################################################################################
# for i in range(0,len(names)):
#     if names[i] in ["lambda","lambda_s","lambda_dpsi","sdev","sdev_s","sdev_dpsi"]:
#         labels[i] = "log("+labels[i]+")"
#         tmp = np.copy(pars[:,i])
#         pars[:,i] = np.log10(tmp)
#         mins[i] = np.log10(abs(mins[i]))
#         maxs[i] = np.log10(abs(maxs[i]))
#         mapping[names[i]] = (mins[i],maxs[i])



#dyo = MCSamples(samples=pars,weights=wei,loglikes=loglike,names=names,labels=labels)
#ena = MCSamples(samples=pars,names=names,labels=labels,settings={"ignore_rows": 0.0,"smooth_scale_2D":0.3,"mult_bias_correction_order":1})
ena = MCSamples(samples=pars,names=names,labels=labels,settings={"ignore_rows": 0.0})
ena.reweightAddingLogLikes(mypost)







### Create the corner plot
###################################################################################################################
#ena = loadMCSamples("./plt_corner",settings={"ignore_rows": 0.0,"smooth_scale_2D":0.9,"mult_bias_correction_order":1})
#ena = loadMCSamples("./plt_corner",settings={"ignore_rows": 0.0,"smooth_scale_2D":0.3,"mult_bias_correction_order":1})

for i in range(0,len(ena.getParamNames().names)):
    ena.getParamNames().names[i].label = labels[i]

g = plots.getSubplotPlotter(subplot_size=4)
g.settings.rcSizes()
g.triangle_plot([ena],param_limits=mapping,filled=True)
#g.triangle_plot([ena],filled=True)


best = ena.getTable(columns=1,include_bestfit=False)
print(best.tableTex())
stats = ena.getMargeStats()
print(stats)


### Plot the true values if available (mock data) or do nothing (if true data)
###################################################################################################################

# if os.path.isfile(path+run+'/vkl_input.json'):
#     # get all the parameter names and values from the true data
#     true_params = {}

#     input_file = path+run+'/vkl_input.json'
#     f          = open(input_file,'r')
#     input_str  = f.read()
#     input_str  = re.sub('//.*?\n|/\*.*?\*/','',input_str,flags=re.S)
#     true_pars  = json.loads(input_str)
    

#     for param in true_pars["physical"]["nlpars"]:
#         true_params[param['nam']] = param['val']

#     for lens in true_pars["lenses"]:
#         for param in true_pars["lenses"][lens]["nlpars"]:
#             true_params[lens+'_'+param['nam']] = param['val']

#     for i in range(0,len(names)):
#         if names[i] not in ["lambda","lambda_s","lambda_dpsi","sdev","sdev_s","sdev_dpsi"]:
#             for ax in g.subplots[i:,i]:
#                 ax.axvline(true_params[names[i]],color='black',ls='solid',alpha=0.5)
#             for ax in g.subplots[i,:i]:
#                 ax.axhline(true_params[names[i]],color='black',ls='solid',alpha=0.5)

'''
for i in range(0,len(names)):
    if names[i] in ["lambda","lambda_s","lambda_dpsi","sdev"]:

        for ax in g.subplots[i+1:,i]:
            ax.set_xscale("log")
            ax.get_xaxis().set_ticks([])
        ax = g.subplots[-1,i]
        ax.set_xscale("log")

        for ax in g.subplots[i,1:i]:
            ax.set_yscale("log")
            ax.get_yaxis().set_ticks([])
        ax = g.subplots[i,0]
        ax.set_yscale("log")
'''

### Export image
g.export('corner.pdf')
#os.remove('plt_corner.txt')
#os.remove('plt_corner.paramnames')
#os.remove('plt_corner.ranges')
#os.remove('plt_corner.py_mcsamples')
