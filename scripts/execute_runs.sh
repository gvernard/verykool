#!/bin/bash
path="RUNS/section_3.4/perturbed_complex_final/"


# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
#names=(`cd $path; ls`)
#delete=(data)
#runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
declare -a runs=("fff_exp_gauss_n3_map" "fxf_exp_gauss_n3_map" "smooth_exp_n1_map" "xff_exp_gauss_n3_map" "xfx_exp_gauss_n3_map" "xxf_exp_curv_n3_map" "xxf_exp_gauss_n3_map" "xxx_exp_gauss_n3_map" "xxx_exp_gauss_n2_map")


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
