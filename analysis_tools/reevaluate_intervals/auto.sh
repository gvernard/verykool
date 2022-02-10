#!/bin/bash

source /home/giorgos/anaconda3/etc/profile.d/conda.sh
conda activate vkl


#path="/home/giorgos/myData/VKL_paper/section_3.4/perturbed_complex_final2/"
#run="fff_exp_gauss_n3/"

#path="/home/giorgos/myData/VKL_paper/section_3.2/new_spiral2/"
#run="modgauss/"

#path="/home/giorgos/myData/VKL_paper/section_3.3/final_upd2/"
#run="xxf_gauss_gauss_n3/"

path=$1
run=$2

read lmodel step <<<$(php ../get_lmodel_step.php $path $run)
echo $lmodel $step
python swap_cols.py ${path}${run} ${lmodel} ${step}
python getdist_intervals.py ${path}${run} ${lmodel} ${step}

rm tmp_postdist.*
