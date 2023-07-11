import os
import sys
import json
import common_funcs
import glob
import shutil


case = sys.argv[1] # full path to the case

if len(sys.argv) > 2:
    mode = sys.argv[2] # mode (map,maxlike,mean)
    if mode not in ["map","maxlike","mean"]:
        print("Second command line argument must be one of: map, maxlike, mean")
        sys.exit()
else:
    mode = "maxlike"
print(mode)    

j_input  = common_funcs.get_json(case+"vkl_input.json",True)


# Find name
pmodel = j_input['parameter_model']
if pmodel == 'both':
    pmodel_name = 'both'
elif pmodel == 'perturbations_standard':
    pmodel_name = 'pert'
elif pmodel == 'standard':
    pmodel_name = 'smooth'


# Find timestep
list_to_check = ['minimizer_output.json']
existing = [os.path.basename(x) for x in glob.glob(case+'output/'+pmodel_name+'_*')]
missing = []
for myfile in list_to_check:
    if pmodel_name+'_'+myfile not in existing:
        missing.append(pmodel_name+'_'+myfile)
    
if len(missing)>0:
    # find timestep
    dum = [os.path.basename(x) for x in glob.glob(case+'output/*_'+pmodel_name+'_model.fits')]
    timesteps = []
    for i in range(0,len(dum)):
        dumdum = dum[i].split('_')
        timesteps.append(int(dumdum[0]))
    timesteps.sort()
    prefix = str(timesteps[-1])+'_'+pmodel_name
else:
    prefix = pmodel_name
print(prefix)

j_output = common_funcs.get_json(case+"output/"+prefix+"_minimizer_output.json",True)
active = j_output["parameters"]


### Set minimizer to test and nproc to 1
j_input["nproc"] = 1
j_input["minimizer"]["type"] = "test"



### Check and update individual groups of parameters

# Lenses
for key,lens in j_input["lenses"].items():
    for obj in lens["nlpars"]:
        full_name = key+"_"+obj["nam"]
        if full_name in active:
            obj["val"] = active[full_name][mode]

# Shear magnitude and angle
for obj in j_input["physical"]["nlpars"]:
    full_name = obj["nam"]
    if full_name in active:
        obj["val"] = active[full_name][mode]

# Source regularization
if "reg_s" in j_input:
    if pmodel in ['standard','smooth']:
        regs_input = j_input["reg_s"]["nlpars"]
    else:
        regs_input = j_input["perturbations"]["reg_s"]["nlpars"]
    for obj in regs_input:
        full_name = obj["nam"]+"_s"
        if full_name in active:
            obj["val"] = active[full_name][mode]

# Perturbations regularization
if "reg_dpsi" in j_input:
    regdpsi_input = j_input["perturbations"]["reg_dpsi"]["nlpars"]
    for obj in regdpsi_input:
        full_name = obj["nam"]+"_dpsi"
        if full_name in active:
            obj["val"] = active[full_name][mode]




### Create the new case dir at the same path as the input case
dum = os.path.normpath(case).split('/')
case_name = os.path.basename(dum[-1])
case_path = '/'.join(dum[:-1])+'/'
new_name = case_name+'_'+mode

if not os.path.exists(case_path+new_name):
    os.makedirs(case_path+new_name)

f = open(case_path+new_name+"/vkl_input.json","w")
json.dump(j_input,f,indent=4,separators=(',',': '))
f.close()

if os.path.exists(case+'coolest_initialization.json'):
    shutil.copyfile(case+'coolest_initialization.json',case_path+new_name+'/coolest_initialization.json')
