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







    # Read vkl_input.json
    options = get_json(path+run_old+"vkl_input.json",True)
    
    

    # Read 'best' parameters from run_old
    vkl_output = get_json(path+run_old+"output/vkl_output.json",False)
    best = vkl_output["json_active"]

    # Update each nlpar list with parameter values from 'best'
    if "physical" in best:
        update_nlpar_list(options["physical"],best["physical"])
    if "other" in best:
        update_nlpar_list(options["reg"]["nlpars"],best["other"])

    if "lenses" in best:
        for i,lens in enumerate(best["lenses"]):
            if len(lens) != 0:
                update_nlpar_list(options["lenses"][i]["nlpars"],lens)



    # Read changes to be made to the updated file
    changes = get_json(changes_file,False)

    # Update the json data recursively
    key_list = flatten_json(changes)
    for key,val in key_list.iteritems():

        dum = key.split("-.-")
        for i,z in enumerate(dum):
            try:
                int(z)
            except ValueError:
                continue
            dum[i] = int(z)

        setInDict(options,dum,val)

    # Force new minimizer, reread in correct order
    if "minimizer" in changes:
        changes_dum = get_json(changes_file,True) # reread in correct order
        options["minimizer"] = changes_dum["minimizer"]




    # Write output file (which will be the new input)
    f = open(path+run_new+"vkl_input.json","w")
    json.dump(options,f,indent=4,separators=(',',': '))
    f.close()





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

def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)

def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value




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



def update_nlpar_list(old,new):
    for par_new in new:
        for par_old in old:
            if par_old["nam"] == par_new["nam"]:
                par_old["val"] = par_new["val"]
                #par_old["err"] = par_new["err"]
                break




if __name__ == "__main__":
    main()

