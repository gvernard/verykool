#!/bin/bash

echo | cat > scripts/report.txt



path="VKL_MODELS/double_ring/dorotas_data/"

# Scan the given path and create alist with all the names (except 'data')
##names=(`find $path -type d -exec basename {} \;`)
#names=(`cd $path; ls`)
#delete=(data)
#runs=( "${names[@]/$delete}" )


# Custom array of run names that are found in 'path'
declare -a runs=("flip_horizontal_psf" "flip_vertical_psf")



for run in "${runs[@]}"; do
    echo $run ...
    echo $run running | cat >> scripts/report.txt
    python scripts/agent.py $path $run/
    echo $run done | cat >> scripts/report.txt
    echo $run done
done
