#!/bin/bash

# arg 1: path
# arg 2: run
# arg 3: json options for theoretical covariance matrices
# arg 4: true source

./theo_matrix_correlation $1$2/output/smooth_source_irregular.dat $3

python image_correlation.py $4
mv image_corr.dat corr_true.dat

../create_gridded_source $1 $2 400 smooth
python image_correlation.py source_model.fits
mv image_corr.dat corr_model.dat

python plot_corr.py
