#!/bin/bash
path="RUNS/VKL_runs/s22_unfiltered_new/"


# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
#names=(`cd $path; ls`)
#delete=(data)
#runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
declare -a runs=("exp_evi" "gauss_evi" "mattern_evi")


echo | cat > scripts/report.txt

for run in "${runs[@]}"; do
#for (( i=0; i<${#runs[@]}; i++ )); do
    #run=${runs[$i]}
    echo $run ...
    echo $run running | cat >> scripts/report.txt
    bin/vkl_driver $path $run/
    echo $run done | cat >> scripts/report.txt
    echo $run done
done
