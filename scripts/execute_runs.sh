#!/bin/bash
path="RUNS/DR_mock_models/slope_-5/"

# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
names=(`cd $path; ls`)
delete=(data)
runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
#declare -a runs=("mn_curv_run" "mn_modg_run" "mn_gaus_run")


echo | cat > scripts/report.txt

for run in "${runs[@]}"; do
#for (( i=0; i<${#runs[@]}; i++ )); do
    #run=${runs[$i]}
    echo $run ...
    echo $run running | cat >> scripts/report.txt
    #python scripts/agent.py $path $run/
    python scripts/agent.py  ${path}$run/ mn_run/
    echo $run done | cat >> scripts/report.txt
    echo $run done
done
