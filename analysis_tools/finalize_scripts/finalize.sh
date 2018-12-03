#!/bin/bash

path=$1
run=$2


#########################################
# Read or determine lmodel and step
if [ "$#" -eq 4 ]; then
    lmodel=$3
    step=$4
    out_path=${path}${run}output/${step}_${lmodel}_
else
    read lmodel step <<<$(php get_lmodel_step.php $path $run $3)
    out_path=${path}${run}output/${lmodel}_
fi
if [ "$lmodel" = "problem" ]; then
    echo "Could not determine lmodel and step, exiting."
    exit
fi



#########################################
# Create output dir
target=analyzed_cases/$run
mkdir -p $target
source activate /data/users/gvernard/myLibraries/anaconda_envs/new_basic


#########################################
echo ' >>> plotting images...'
if [ "$lmodel" = "smooth" ]; then
    python plot_all_smooth.py $path $run $step
else
    python plot_all_pert.py $path $run $step
fi


#########################################
echo ' >>> making corner plot...'
python plot_corner.py $path $run $lmodel $step


#########################################
echo ' >>> producing latex tables...'
php latex.php $path $run $lmodel $step


#########################################
echo ' >>> producing latex document...'
pdflatex report.tex > /dev/null 
cp report.pdf ${path}${run}report.pdf
mv report.pdf $target


#########################################
echo ' >>> creating regularly gridded source...'
./create_gridded_source $path $run 200 $lmodel $step # 3rd argument is the resolution of the gridded source
mv source_model.fits $target/source_model.fits


#########################################
echo ' >>> creating parameter table.txt...'
php get_par_table.php $path $run $lmodel $step
mv table_pars.txt $target/


#########################################
echo ' >>> creating masked residuals...'
if [ "$lmodel" = "pert" ]; then
    python plot_masked_residuals.py ${path}${run}output/smooth_residual.fits ${path}data/mask.fits
    mv masked_residual.fits $target/smooth_masked_residual.fits
fi
python plot_masked_residuals.py ${out_path}residual.fits ${path}data/mask.fits
mv masked_residual.fits $target/masked_residual.fits


#########################################
echo ' >>> collecting and packaging remaining files...'
# move the plots after the latex document is produced
mv all.png $target/
mv all.pdf $target/
if [ -f corner.pdf ]; then
    mv corner.pdf $target/
fi

cp ${path}data/mask.fits $target/
cp ${path}data/psf.fits $target/
cp ${path}data/noise.* $target/
cp ${path}data/image.fits $target/
cp ${out_path}model.fits $target/model.fits
cp ${out_path}residual.fits $target/residual.fits
if [ "$lmodel" = "pert" ]; then
    cp ${path}${run}output/smooth_model.fits $target/smooth_model.fits
    cp ${path}${run}output/smooth_residual.fits $target/smooth_residual.fits
fi



conda deactivate
echo "Done!"
