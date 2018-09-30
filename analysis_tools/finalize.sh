#!/bin/bash

path=$1
run=$2
target=$run

if [ "$#" -eq 3 ]; then
    step=$3
    out_path=${path}${run}output/${step}_
else
    step=''
    out_path=${path}${run}output/
fi



mkdir -p $target
source activate /data/users/gvernard/myLibraries/anaconda_envs/new_basic

#########################################
echo ' >>> plotting images...'
python plot_all.py $path $run $step

#########################################
echo ' >>> making corner plot...'
python plot_corner.py $path $run $step

#########################################
echo ' >>> producing latex tables...'
php latex.php $path $run $step

#########################################
echo ' >>> producing latex document...'
pdflatex report.tex > /dev/null 
cp report.pdf ${path}${run}report.pdf

#########################################
echo ' >>> collecting and packaging code output...'
cp ${path}data/mask.fits $target/
cp ${path}data/psf.fits $target/
cp ${path}data/noise.* $target/
cp ${path}data/image.fits $target/
cp ${out_path}vkl_image.fits $target/model.fits
cp ${out_path}vkl_residual.fits $target/residual.fits

./create_gridded_source ${out_path} 200 # second argument is the resolution of the gridded source
mv vkl_source.fits $target/source.fits

php get_par_table.php $path $run $step



mv all.png $target
mv corner.pdf $target
mv report.pdf $target
mv table_pars.txt $target/




conda deactivate
echo "Done!"
