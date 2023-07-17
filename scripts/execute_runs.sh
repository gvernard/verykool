#!/bin/bash
path="/home/gvernardos/s22_unfiltered_new/"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # path to this script

# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
#names=(`cd $path; ls`)
#delete=(data)
#runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
declare -a runs=("curv_evi" "exp_evi" "gauss_evi" "mattern_evi")


echo | cat > ${SCRIPT_DIR}/report.txt

for run in "${runs[@]}"; do
#for (( i=0; i<${#runs[@]}; i++ )); do
    #run=${runs[$i]}
    echo $run ...
    echo $run running | cat >> ${SCRIPT_DIR}/report.txt
    #bin/vkl_driver $path $run/
    sbatch ${SCRIPT_DIR}/vkl_cluster_submit.sh ${path} ${run}/
    echo $run done | cat >> ${SCRIPT_DIR}/report.txt
    echo $run done
done
