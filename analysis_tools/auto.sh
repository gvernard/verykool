#!/bin/bash

suffix="_all.png"
for i in `seq 1 30`;
do
    python plot_all.py /net/argo/data/users/gvernard/RESULTS/VeryKooL_DATA/double_ring/dorotas_data/ degen_test/ $i
    mv all.png $i$suffix
done

