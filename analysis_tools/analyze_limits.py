import numpy as np
import sys
import json
import os
import re
from getdist import MCSamples


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
inpath = os.path.join(path,run,'output/' + lmodel)
if len(sys.argv) > 3:
    step = sys.argv[3]
    inpath_t  = os.path.join(path,run,'output',str(step) + '_' + lmodel)
    outpath = os.path.join(outpath,str(step) + '_')
else:
    step = ''
    inpath_t = inpath




# Posterior distributions
ola = np.loadtxt(inpath_t+'_postdist.txt')
part_a = ola[:,2:4]
part_b = ola[:,4:9]
part_c = ola[:,9:]
tmp = np.column_stack((ola[:,:2],part_b,part_a,part_c))
pars = tmp[:,2:]
wei  = tmp[:,0]
mypost = tmp[:,1]




old_lines = [line.rstrip('\n') for line in open(inpath+'_postdist.paramnames')]
tmp = ['a' for i in range(0,len(old_lines)-1)]
tmp[0:5] = old_lines[2:7]
tmp[5:7] = old_lines[0:2]
tmp[7:9] = old_lines[7:9]
tmp[9:11] = old_lines[9:11]

names = []
labels = []
for line in tmp:
    dum = re.sub(r'\s+', ',', line).split(',')
    names.append(dum[0])
    labels.append(dum[1])
labels = [w.replace('_', '-') for w in labels]



ena = MCSamples(samples=pars,names=names,labels=labels,settings={"ignore_rows": 0.0,"smooth_scale_1D":0.5,"mult_bias_correction_order":1})
#min_non_zero = np.min(mypost[np.nonzero(mypost)])
#mypost = np.where(mypost<min_non_zero,min_non_zero,mypost)
ena.reweightAddingLogLikes(mypost)



mystats = ena.getMargeStats()

json_out = {}
for name in names:
    mypar = mystats.parWithName(name)
    json_out[name] = {
        "mean": mypar.mean,
        "s1_lower": mypar.limits[0].lower,
        "s1_upper": mypar.limits[0].upper,
        "s2_lower": mypar.limits[1].lower,
        "s2_upper": mypar.limits[1].upper,
        "s3_lower": mypar.limits[2].lower,
        "s3_upper": mypar.limits[2].upper        
    }


#print(myjson["collapsed_active"])
with open(outpath+'limits.json','w') as outfile:
    json.dump(json_out,outfile,indent=4)
