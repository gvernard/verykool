#!/bin/bash

suffix="_all.png"
for i in `seq 1 30`;
do
    python plot_all.py /net/argo/data/users/gvernard/RESULTS/VKL_MODELS/double_ring/new_perturbations/pert_3_4/ cov_n3/ $i
    mv all.png $i$suffix
done

