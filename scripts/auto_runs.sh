#!/bin/bash

path="/net/argo/data/users/gvernard/RESULTS/VeryKooL_DATA/double_ring/dorotas_data/"

touch ~/myGit/VeryKooL/scripts/report.txt
cd $path
for run in new_degen_test_[0-9]
do
#    echo $run
    echo $run running | cat >> ~/myGit/VeryKooL/scripts/report.txt 
    python ~/myGit/VeryKooL/scripts/agent.py $path $run/
    echo $run done | cat >> ~/myGit/VeryKooL/scripts/report.txt 
done
