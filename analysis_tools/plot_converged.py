import os
import sys
import re
import numpy as np
import json
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


path = sys.argv[1]
run  = sys.argv[2]

# Determine the number of iterations
files = [f for f in os.listdir(path+run+"output/") if re.match(r'[0-9]+_vkl_output.json',f)]
ind = np.zeros(len(files))
for i in range(0,len(files)):
    dum = files[i].split('_')
    ind[i] = int(dum[0])

iters = range(1,int(max(ind))+1)
#print iters


# Read active parameter ranges and names
names,dum_mins,dum_maxs = np.loadtxt(path+run+'output/plt_corner.ranges',dtype={'names': ('names','min','max'),'formats':('|S15',np.float,np.float)},unpack=True)
mins = {}
maxs = {}
for i in range(0,len(names)):
    mins[names[i]] = dum_mins[i]
    maxs[names[i]] = dum_maxs[i]


# split parameters to active (linear) and regularization (log)
f         = open(path+run+'output/1_vkl_output.json','r')
input_str = f.read()
myjson0   = json.loads(input_str)

term_names = ["like","chi2","reg"]
phys_names = [ f["nam"] for f in myjson0["json_active"]["physical"] ]
reg_names = [ f["nam"] for f in myjson0["json_active"]["reg"] ]
lenses_names = {}
for lens in myjson0["json_active"]["lenses"]:
    lenses_names[lens] = [ f["nam"] for f in myjson0["json_active"]["lenses"][lens] ]




terms  = np.zeros(shape=(len(term_names),len(iters)))
phys   = np.zeros(shape=(len(phys_names),len(iters)))
reg    = np.zeros(shape=(len(reg_names),len(iters)))
lenses = {}
for lens in lenses_names:
    lenses[lens] = np.zeros(shape=(len(lenses_names[lens]),len(iters)))

for i in iters:
    f         = open(path+run+'output/'+str(i)+'_vkl_output.json','r')
    input_str = f.read()
    myjson    = json.loads(input_str)

    for j in range(0,len(term_names)):
        key = term_names[j]
        terms[j][i-1] = float(myjson["terms"][key])

    for j in range(0,len(phys_names)):
        key = phys_names[j]
        phys[j][i-1] = (float(myjson["full_active"][key]) - mins[key])/(maxs[key] - mins[key])

    for j in range(0,len(reg_names)):
        key = reg_names[j]
        if key[:6] == "lambda":
            reg[j][i-1] = np.log10(float(myjson["full_active"][key]))
        else:
            reg[j][i-1] = (float(myjson["full_active"][key]) - mins[key])/(maxs[key] - mins[key])

    for lens in lenses_names:
        for j in range(0,len(lenses_names[lens])):
            key = lens+"_"+lenses_names[lens][j]
            lenses[lens][j][i-1] = (float(myjson["full_active"][key]) - mins[key])/(maxs[key] - mins[key])




fig = plt.figure(figsize=(8,20))
gs = gridspec.GridSpec(3+len(lenses_names),1)


# plot likelihood
pterms = plt.subplot(gs[0])
for j in range(0,len(term_names)):
    pterms.plot(iters,terms[j],label=term_names[j])
pterms.legend(loc='upper left')
pterms.set_ylabel("evidence")
pterms.set_xlabel("iteration")


# plot regularization parameters
preg = plt.subplot(gs[1])
for j in range(0,len(reg_names)):
    preg.plot(iters,reg[j],label=reg_names[j])
preg.legend(loc='upper left')
preg.set_ylabel("log parameters")
preg.set_xlabel("iteration")


# plot physical parameters (gamma,phi_gamma) normalized in (0,1)
pphys = plt.subplot(gs[2])
for j in range(0,len(phys_names)):
    pphys.plot(iters,phys[j],label=phys_names[j])
pphys.set_ylim(0,1)
pphys.legend(loc='upper left')
pphys.set_ylabel("physical")
pphys.set_xlabel("iteration")


# plot mass model parameters for each lens normalized in (0,1)
i = 0
for lens in lenses_names:
    dum = plt.subplot(gs[3+i])
    for j in range(0,len(lenses_names[lens])):
        dum.plot(iters,lenses[lens][j],label=lenses_names[lens][j])
    dum.set_ylim(0,1)
    dum.legend(loc='upper left')
    dum.set_ylabel(lens)
    dum.set_xlabel("iteration")
    i += 1



plt.tight_layout()
plt.savefig('converged.pdf',bbox_inches='tight')
