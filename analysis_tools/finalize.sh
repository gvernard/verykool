#!/bin/bash
# command line argument 1: path to root directory of a run of VeryKooL code, e.g. /smth/smth/../mysimX/
# command line argument 2: step of the partial code output that is reported

#path2data=$1data/
#path2out=$1output/

#ds9 $path2data/image.fits -export png image.png -quit
#ds9 $path2out/vkl_image.fits -export png vkl_image.png -quit
#ds9 $path2out/vkl_residual.fits -export png vkl_residual.png -quit
#ds9 $path2data/source.fits -export png source.png -quit
#python plot_voronoi.py $path2out
#mv voronoi.png vkl_source.png
#convert vkl_source.png -resize 110x100 vkl_source.png

#montage image.png vkl_image.png vkl_residual.png source.png vkl_source.png -mode Unframe -tile 3x -background WHITE collage.png
#eog collage.png




#eog all.png &
#eog corner.png &


echo ' >>> plotting images...'
python plot_all.py $1 $2
echo ' >>> making corner plot...'
cp $1output/$2_corner.txt $1output/corner.txt
python plot_corner.py $1
echo ' >>> producing latex tables...'
php latex.php $1 $2
echo ' >>> producing latex document...'
pdflatex report.tex > /dev/null 

mv report.pdf $1

evince $1report.pdf

