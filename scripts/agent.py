# Select which VeryKool version to run based on the minimizer type
import os
import sys
import cmd
import json
import re
import shutil
import subprocess




def main():
    vkl_dir       = "/Users/users/gvernard/myGit/VeryKooL/"
    cosmo_lib_dir = "/net/argo/data/users/gvernard/myLibraries/cosmosis/"
    conda_env     = "/net/argo/data/users/gvernard/myLibraries/anaconda_envs/cosmosis_env"
    
    original = ['test','multinest','simplex']
    cosmosis = ['cosmosis_test','cosmosis_metropolis','cosmosis_maxlike','cosmosis_minuit','cosmosis_emcee','cosmosis_multinest']
    
    path = sys.argv[1]
    run  = sys.argv[2]
    
    
    
    
    
    
    # Check/create relative directories and files
    ##############################################################################################
    if not os.path.exists(path+run+"vkl_input.json"):
        print "'vkl_input.json' file not found!"
        sys.exit()
        
    if os.path.isdir(path+run+"output"):
        answer = ''
        while answer != 'y' and answer != 'n':
            answer = raw_input("output directory exists, cleanup? (y/n): ")
            
        if answer == 'y':
            for item in os.listdir(path+run+"output"):
                item_path = os.path.join(path+run+"output",item)
                try:
                    if os.path.isfile(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path): shutil.rmtree(item_path)
                except Exception as e:
                    print(e)
    else:
        os.mkdir(path+run+"output")
                            



    # Read input options
    ##############################################################################################
    f = open(path+run+"vkl_input.json","r")
    input_str = f.read()
    input_str = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
    input_str = re.sub(re.compile("//.*?\n" ),"",input_str)
    options   = json.loads(input_str)
    
    if options["minimizer"]["type"] in original:
        mode = "original"
    elif options["minimizer"]["type"] in cosmosis:
        mode = "cosmosis"
        
            
        
        
    # Directory and file check based on the input options
    ##############################################################################################
    if mode == "cosmosis":
        create_pipeline(path,run,options,vkl_dir)
        subprocess.call([vkl_dir+"bin/createCosmosisValuesPriorsIni",path,run])




    # Execute main code
    ##############################################################################################
    if mode == "original":
        msg = "Executing original code, good luck:"
        cmd = "mpirun -np 1 " + vkl_dir + "bin/verykool " + path + " " + run
    else:
        msg = "Executing cosmosis code, good luck:"
        create_bash_script(path,run,options,cosmo_lib_dir,conda_env,vkl_dir)
        cmd = "bash " + path + run + "auto_run_cosmosis.sh"

        
        
        
    print msg
    os.system(cmd)
    print "Execution succesful (results may still be wrong though...)"




def create_bash_script(path,run,options,cosmo_lib_dir,conda_env,vkl_dir):
    f = open(path+run+"auto_run_cosmosis.sh","w")
    f.write("#!/bin/bash\n")
    f.write("# This is an automatically generated file!\n")
    f.write("source " + cosmo_lib_dir +"setup-my-cosmosis\n")
    f.write("source activate " + conda_env + "\n")
    #f.write("mpirun -np 1 cosmosis --mpi " + path + run + "cosmosis_pipeline.ini\n")
    f.write("cosmosis " + path + run + "cosmosis_pipeline.ini\n")
    f.close()


def create_pipeline(path,run,options,vkl_dir):
    f = open(path+run+"cosmosis_pipeline.ini","w")
    f.write("[runtime]\n")
    f.write("sampler = " + options["minimizer"]["type"][9:] + "\n")
    f.write("\n")

    f.write("["+options["minimizer"]["type"][9:]+"]\n")
    if options["minimizer"]["type"] in options:
        for key,value in options[options["minimizer"]["type"]].iteritems():
            f.write(key + " = " + str(value) + "\n")
    f.write("fatal_errors = T\n")
    f.write("\n")

    f.write("[output]\n")
    f.write("filename = " + path + run + "output/cosmosis_output.dat\n")
    f.write("format = text\n")
    f.write("\n")

    f.write("[pipeline]\n")
    f.write("modules = verykool\n")
    f.write("values = " + path + run + "cosmosis_values.ini\n")
    f.write("likelihoods = verykool\n")
    f.write("quiet = T\n")
    f.write("timing = F\n")
    f.write("debug = F\n")
    f.write("\n")

    f.write("[verykool]\n")
    f.write("file = " + vkl_dir + "lib/libverykool.so\n")
    f.write("path = " + path + "\n")
    f.write("run = " + run + "\n")
    f.close()


def create_values(path,run,options):
    f = open(path+run+"cosmosis_values.ini","w")
    f.write("[nlpars]\n")
    for nlpar in options["physical"]:
        f.write(nlpar["nam"] + " = " + str(nlpar["min"]) + " " + str(nlpar["val"]) + " " + str(nlpar["max"]) + "\n")
    for nlpar in options["reg"]["nlpars"]:
        f.write(nlpar["nam"] + " = " + str(nlpar["min"]) +  " " +str(nlpar["val"]) + " " + str(nlpar["max"]) + "\n")
    for i,lens in enumerate(options["lenses"]):
        for nlpar in lens["nlpars"]:
            f.write(str(i)+"_" + nlpar["nam"] + " = " + str(nlpar["min"]) + " " + str(nlpar["val"]) + " " + str(nlpar["max"]) + "\n")
    f.close()


def create_priors(path,run,options):
    f = open(path+run+"cosmosis_priors.ini","w")
    f.write("[nlpars]\n")
    for nlpar in options["physical"]:
        f.write(nlpar["nam"] + " = " + nlpar_prior(nlpar) + "\n")
    for nlpar in options["reg"]["nlpars"]:
        f.write(nlpar["nam"] + " = " + nlpar_prior(nlpar) + "\n")
    for i,lens in enumerate(options["lenses"]):
        for nlpar in lens["nlpars"]:
            f.write(str(i) + "_" + nlpar["nam"] + " = " + nlpar_prior(nlpar) + "\n")
    f.close()
def nlpar_prior(nlpar):
    out  = ""
    if nlpar["pri"]["type"] == "uni":
        out = "uniform " + str(nlpar["min"]) + " " + str(nlpar["max"])
    elif nlpar["pri"]["type"] == "gauss":
        out = "gauss " + str(nlpar["pri"]["mean"]) + " " + str(nlpar["pri"]["sdev"])
    elif nlpar["pri"]["type"] == "exp":
        out = "exponential " + str(nlpar["pri"]["beta"])
    return out





if __name__ == "__main__":
    main()

