#!/bin/bash

echo | cat > scripts/report.txt



path="VKL_MODELS/VKL_paper/koop2005/"

# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
#names=(`cd $path; ls`)
#delete=(data)
#runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
declare -a runs=("mn_curv_run" "mn_modg_run" "mn_gaus_run")

#for run in "${runs[@]}"; do
for (( i=0; i<${#runs[@]}; i++ )); do
    run=${runs[$i]}
    echo $run ...
    echo $run running | cat >> scripts/report.txt
    #python scripts/agent.py $path $run/
    echo $run done | cat >> scripts/report.txt
    echo $run done
done
