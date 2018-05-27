import os
import sys
import json
import common_funcs

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
options = common_funcs.get_json(path+run_old+"vkl_input.json",True)

# Read changes to be made to vkl_input.json
changes = common_funcs.get_json(changes_file,True)



    
# Loop over the first level of changes and apply them
#########################################################################
# For the list of nlpars
if "nlpars" in changes:
    for list_name in changes["nlpars"]:
        if list_name == "lenses":
            for lens in changes["nlpars"]["lenses"]:
                if lens in options["lenses"]:
                    common_funcs.update_nlpar_list(changes["nlpars"]["lenses"][lens],options["lenses"][lens]["nlpars"])
                else:
                    print "Invalid lens to change: ",lens," !"
        else:
            common_funcs.update_nlpar_list(changes["nlpars"][list_name],options[list_name]["nlpars"])

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

