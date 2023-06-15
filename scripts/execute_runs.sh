#!/bin/bash
path="RUNS/VKL_runs/sourcemag22_unfiltered_min/"


# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
#names=(`cd $path; ls`)
#delete=(data)
#runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
declare -a runs=("N_exp_n2_log_A" "N_exp_n2_log_B" "N_exp_n2_log_C" "N_exp_n2_log_D" "N_gau_n2_log_A" "N_gau_n2_log_B" "N_gau_n2_log_C" "N_gau_n2_log_D" "N_cur_n2_log_A" "N_cur_n2_log_B" "N_cur_n2_log_C" "N_cur_n2_log_D")


echo | cat > scripts/report.txt

mkdir $path/collected_plots

for run in "${runs[@]}"; do
#for (( i=0; i<${#runs[@]}; i++ )); do
    #run=${runs[$i]}
    echo $run ...
    echo $run running | cat >> scripts/report.txt
    python scripts/agent.py $path $run/
    #python scripts/agent.py  ${path}$run/ base_run/

    cp ${path}${run}/output/all.png $path/collected_plots/${run}.png
    
    echo $run done | cat >> scripts/report.txt
    echo $run done
done
