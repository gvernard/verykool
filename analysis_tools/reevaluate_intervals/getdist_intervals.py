from getdist import plots,chains,MCSamples,loadMCSamples
import numpy as np
import sys
import json
import os.path
import matplotlib
import re
from shutil import copyfile
import re

case = sys.argv[1]
lmodel = sys.argv[2]
step = sys.argv[3]

# Get json to replace
input_file = case+'output/'+str(step)+'_'+str(lmodel)+'_output.json'
f          = open(input_file,'r')
input_str  = f.read()
input_str  = re.sub('//.*?\n|/\*.*?\*/','',input_str,flags=re.S)
myjson  = json.loads(input_str)





names,labels = np.genfromtxt('tmp_postdist.paramnames',dtype='str',unpack=True)
labels = [w.replace('_', '-') for w in labels]

ola = np.loadtxt("tmp_postdist.txt")
pars = ola[:,2:]
wei  = ola[:,0]
mypost = ola[:,1]


#ena = MCSamples(samples=pars,names=names,labels=labels,settings={"ignore_rows": 0.0,"fine_bins_1D":800,"smooth_scale_1D":0.5,"mult_bias_correction_order":1})
ena = MCSamples(samples=pars,names=names,labels=labels,settings={"ignore_rows": 0.0,"smooth_scale_1D":0.5,"mult_bias_correction_order":1})

min_non_zero = np.min(mypost[np.nonzero(mypost)])
mypost = np.where(mypost<min_non_zero,min_non_zero,mypost)
ena.reweightAddingLogLikes(mypost)


#print(ena.getTable(limit=1).tableTex())
mystats = ena.getMargeStats()
print(mystats)
for p in names:
    #orig = ena.getInlineLatex(p,limit=1)
    #print('-----',orig)
    mypar = mystats.parWithName(p)

    if p in ["lambda","lambda_s","lambda_dpsi","sdev","sdev_s","sdev_dpsi"]:
        mean     = mypar.mean
        s1_lower = mypar.mean - mypar.err
        s1_upper = mypar.mean + mypar.err
        # mean     = np.power(10.0,mypar.mean)
        # s1_lower = np.power(10.0,mypar.limits[0].lower)
        # s1_upper = np.power(10.0,mypar.limits[0].upper)
    else:
        mean     = mypar.mean
        s1_lower = mypar.limits[0].lower
        s1_upper = mypar.limits[0].upper

        
    #print(p,mypar.limits[0].limitType(),mean,s1_upper,s1_lower)
    myjson["collapsed_active"][p]["mean"] = mean
    myjson["collapsed_active"][p]["s1_low"] = s1_lower
    myjson["collapsed_active"][p]["s1_high"] = s1_upper

#print(myjson["collapsed_active"])
with open(case+'output/getdist_output.json','w') as outfile:
    json.dump(myjson,outfile,indent=4)
