import re 
import json
import sys
import corner
import numpy as np
import matplotlib
import matplotlib.pylab as plt
matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman'],'size':20})
matplotlib.rc('text', usetex=True)

rootdir = sys.argv[1]

par_names  = np.genfromtxt(rootdir+'output/corner.paramnames',usecols=[0],dtype='str')
par_ranges = np.genfromtxt(rootdir+'output/corner.ranges',usecols=[1,2])
dum        = np.genfromtxt(rootdir+'output/corner.txt')
samples = dum[:,2:]
weights = dum[:,1]





# get all the parameter names and values from the true data
true_params = {}

#with open(rootdir+'data/params.json') as json_input:    
#tmp = rootdir.split('/')
#input_file = rootdir+'/'+tmp[-2]+'.json'
input_file = rootdir+'/vkl_input.json'
f          = open(input_file,'r')
input_str  = f.read()
input_str  = re.sub(r'\\\n', '', input_str)
input_str  = re.sub(r'//.*\n', '', input_str)
true_pars  = json.loads(input_str)


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
true_values = np.empty([0])
for name in par_names:
    if name == 'lambda':
        final_params = np.append(final_params,name)
        true_values = np.append(true_values,None)
    else:
        final_params = np.append(final_params,name)
        true_values = np.append(true_values,true_params[name])









figure = corner.corner(samples,weights=weights,range=par_ranges,labels=par_names,bins=20,max_n_ticks=5,plot_datapoints=False,plot_density=False,truths=true_values,smooth=2,fill_contours=True,squeeze=True)

# If i don't draw the figure the tick labels are not populated and will be emtpy
#figure.canvas.draw()


#fig,ax = plt.subplots()
ax = figure.get_axes()


for a in ax:
    dum = a.get_yticklabels() # just a way to check if the plot exists
    if len(dum)>0:
        ymin,ymax = a.get_ylim()
        step = (ymax-ymin)/4
        ticks = np.arange(ymin+step,ymax,step)
        a.set_yticks(ticks)
        
    
    dum = a.get_xticklabels() # just a way to check if the plot exists
    if len(dum)>0:
        xmin,xmax = a.get_xlim()
        step = (xmax-xmin)/4
        ticks = np.arange(xmin+step,xmax,step)
        a.set_xticks(ticks)





#    a.set_ticklabels(str(ticks))

'''
    # get all the labels of Y axis
    labels = [item.get_text() for item in a.get_yticklabels()]
    if len(labels) > 0:
        # remove the first and the last labels
        labels[-1] = ''
        labels[-2] = ''
        # set these new labels

    a.set_yticklabels(labels)


    # get all the labels of X axis
    labels = [item.get_text() for item in a.get_xticklabels()]
    if len(labels) > 0:
        # remove the first and the last labels
        labels[-1] = ''
        labels[-2] = ''
        # set these new labels
        a.set_xticklabels(labels)
'''


plt.savefig('corner.png')
