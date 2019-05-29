#!/bin/bash
path="RUNS/DR_exclusion_plot/slope_-3/"


# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
names=(`cd $path; ls`)
delete=(data)
runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
#declare -a runs=("mock_03" "mock_06" "mock_09" "mock_12" "mock_15" )
#declare -a runs=("model_5" "model_6" "model_7" )


echo | cat > scripts/report.txt

for run in "${runs[@]}"; do
#for (( i=0; i<${#runs[@]}; i++ )); do
    #run=${runs[$i]}
    echo $run ...
    echo $run running | cat >> scripts/report.txt
    #python scripts/agent.py $path $run/
    python scripts/agent.py  ${path}$run/ base_run/
    echo $run done | cat >> scripts/report.txt
    echo $run done
done
