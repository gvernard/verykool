import json
import re
from collections import OrderedDict
from functools import reduce  # forward compatibility for Python 3
import operator




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
