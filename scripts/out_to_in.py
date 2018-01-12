import os
import sys
import json
import re
from collections import OrderedDict
from functools import reduce  # forward compatibility for Python 3
import operator

def main():
    path         = sys.argv[1]
    run_old      = sys.argv[2]
    run_new      = sys.argv[3]
    changes_file = sys.argv[4]
    



    # Check if directories exist and create new ones accordingly
    #########################################################################
    if not os.path.isdir(path+run_old):
        print "Old run: '" + run_old + "' does not exist!"
        sys.exit()

    if os.path.isdir(path+run_new):
        answer = ''
        while answer != 'y' and answer != 'n':
            answer = raw_input("New run: '" + run_new + "' already exists, overwrite it? (y/n): ")
        if answer == 'n':
            sys.exit();
    else:
        os.mkdir(path+run_new)




    # Read the various files
    #########################################################################
    # Read vkl_input.json
    options = get_json(path+run_old+"vkl_input.json",True)

    # Read 'best' parameters from run_old
    vkl_output = get_json(path+run_old+"output/vkl_output.json",False)
    best = vkl_output["json_active"]
  
    # Read changes to be made to vkl_input.json
    changes = get_json(changes_file,True)




    # Update all the parameter values (the old active parameters)
    #########################################################################
    for list_name in best:
        if list_name == "lenses":
            for i,lens in enumerate(best["lenses"]):
                if "priors" in changes:
                    update_nlpar_with_prior(options["lenses"][str(lens)]["nlpars"],best["lenses"][lens],changes["priors"])
                else:
                    update_nlpar(options["lenses"][str(lens)]["nlpars"],best["lenses"][lens])
        else:
            if "priors" in changes:
                update_nlpar_with_prior(options[list_name]["nlpars"],best[list_name],changes["priors"])
            else:
                update_nlpar(options[list_name]["nlpars"],best[list_name])

    
    # Loop over the first level of changes and apply them
    #########################################################################
    # For the list of nlpars
    if "nlpars" in changes:
        for list_name in changes["nlpars"]:
            if list_name == "lenses":
                for lens in changes["nlpars"]["lenses"]:
                    if lens in options["lenses"]:
                        update_nlpar_list(changes["nlpars"]["lenses"][lens],options["lenses"][lens]["nlpars"])
                    else:
                        print "Invalid lens to change: ",lens," !"
            else:
                update_nlpar_list(changes["nlpars"][list_name],options[list_name]["nlpars"])

    # For all the other first level options in vkl_input.json
    for change in changes:
        if change not in ["priors","nlpars"]:
            if change in options:
                options[change] = changes[change]
            else:
                print "Invalid option to change: ",change," !"
                sys.exit()




    # Write output file (which will be the new input)
    #########################################################################
    f = open(path+run_new+"vkl_input.json","w")
    json.dump(options,f,indent=4,separators=(',',': '))
    f.close()











'''
def flatten_json(y):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '-.-')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '-.-')
                i += 1
        else:
            out[name[:-3]] = x

    flatten(y)
    return out

def update_recursively(keylist,target):
    for key,val in keylist.iteritems():
        dum = key.split("-.-")
        for i,z in enumerate(dum):
            try:
                int(z)
            except ValueError:
                continue
            dum[i] = int(z)
            
        setInDict(target,dum,val)

def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value

def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)
'''



def get_json(filename,ordered=True):
    f = open(filename,"r")
    dum = f.read()
    f.close()
    dum = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",dum)
    dum = re.sub(re.compile("//.*?\n" ),"",dum)

    if ordered:
        json_obj = json.loads(dum,object_pairs_hook=OrderedDict)
    else:
        json_obj = json.loads(dum)

    return json_obj




def update_nlpar_list(new_list,old_list):
    for new_par in new_list:
        for old_par in old_list:
            if new_par["nam"] == old_par["nam"]:
                for x in new_par:
                    old_par[x] = new_par[x]
                break



def update_nlpar(old,new):
    for par_new in new:
        for par_old in old:
            if par_old["nam"] == par_new["nam"]:
                par_old["val"] = round(par_new["val"],5)
                #par_old["err"] = par_new["err"]
                break


def update_nlpar_with_prior(old,new,prior):
    for par_new in new:
        for par_old in old:
            if par_old["nam"] == par_new["nam"]:
                par_old["val"] = round(par_new["val"],5)
                if prior["type"] == "gauss":
                    par_old["pri"]["type"] = "gauss"
                    par_old["pri"]["mean"] = round(par_new["val"],5)
                    par_old["pri"]["sdev"] = round((par_new["max"] - par_new["min"])*prior["factor"],5)
                elif prior["type"] == "uni":
                    par_old["pri"]["type"] = "uni"
                    new_ran = (par_old["max"] - par_old["min"])*prior["factor"]
                    new_max = round(par_old["val"] + new_ran/2.0,5)
                    if new_max < par_old["max"]:
                        par_old["max"] = new_max
                    new_min = round(par_old["val"] - new_ran/2.0,5)
                    if new_min > par_old["min"]:
                        par_old["min"] = new_min
                break





if __name__ == "__main__":
    main()

