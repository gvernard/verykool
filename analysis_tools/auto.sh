#!/bin/bash

suffix="_all.png"
for i in `seq 1 9`;
do
    python plot_all.py /net/argo/data/users/gvernard/RESULTS/VeryKooL_DATA/J0946+1006/ ada3/ $i
    mv all.png $i$suffix
done

