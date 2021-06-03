#!/bin/bash
path="RUNS/section_3.3/final_upd/"


# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
#names=(`cd $path; ls`)
#delete=(data)
#runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
declare -a runs=("fff_gauss_gauss_n3_map" "xff_gauss_curv_n3_map" "xff_gauss_gauss_n3_map" "xxf_gauss_gauss_n3_map")


echo | cat > scripts/report.txt

for run in "${runs[@]}"; do
#for (( i=0; i<${#runs[@]}; i++ )); do
    #run=${runs[$i]}
    echo $run ...
    echo $run running | cat >> scripts/report.txt
    python3 scripts/agent.py $path $run/
    #python scripts/agent.py  ${path}$run/ base_run/
    echo $run done | cat >> scripts/report.txt
    echo $run done
done
