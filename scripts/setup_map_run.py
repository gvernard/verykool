import os
import sys
import json
import common_funcs
import glob



case = sys.argv[1] # full path to the case



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
list_to_check = ['output.json']
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

j_output = common_funcs.get_json(case+"output/"+prefix+"_output.json",True)
active = j_output["json_active"]

#keys = j_output['collapsed_active'].keys()
#maps = {}
#for key in keys:
#    maps[key] = round(j_output['collapsed_active'][key]["map"],5)
#print(maps)


### Set minimizer to test and nproc to 1
j_input["nproc"] = 1
j_input["minimizer"]["type"] = "test"



### Check and update individual groups of parameters

# Lenses
for key,lens in active["lenses"].items():
    if len(lens) != 0:
        lenses_input = j_input["lenses"][key]["nlpars"]
        for obj in lens:
            for in_obj in lenses_input:
                if obj["nam"] == in_obj["nam"]:
                    in_obj["val"] = obj["map"]
                    break
    
# shear magnitude and angle
if "physical" in active:
    for obj in active["physical"]:
        for in_obj in j_input["physical"]["nlpars"]:
            if obj["nam"] == in_obj["nam"]:
                in_obj["val"] = obj["map"]
                break

# source regularization            
if "reg_s" in active:
    if pmodel == 'standard':
        regs_input = j_input["reg_s"]["nlpars"]
    else:
        regs_input = j_input["perturbations"]["reg_s"]["nlpars"]
    for obj in active["reg_s"]:
        for in_obj in regs_input:
            dum = obj["nam"].split('_')
            name = dum[0]
            if in_obj["nam"] == name:
                in_obj["val"] = obj["map"]
                break

# perturbations regularization
if "reg_dpsi" in active:
    regdpsi_input = j_input["perturbations"]["reg_dpsi"]["nlpars"]
    for obj in active["reg_dpsi"]:
        for in_obj in regdpsi_input:
            dum = obj["nam"].split('_')
            name = dum[0]
            if in_obj["nam"] == name:
                in_obj["val"] = obj["map"]
                break


# Create the new case dir at the same path as the input case
dum = os.path.normpath(case).split('/')
case_name = os.path.basename(dum[-1])
case_path = '/'.join(dum[:-1])+'/'
new_name = case_name+'_map'

if not os.path.exists(case_path+new_name):
    os.makedirs(case_path+new_name)

f = open(case_path+new_name+"/vkl_input.json","w")
json.dump(j_input,f,indent=4,separators=(',',': '))
f.close()

