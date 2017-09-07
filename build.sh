#!/bin/bash -f
mpic++ -std=c++11 -fopenmp -o verykool fastell.o imagePlane.cpp sourcePlane.cpp massModels.cpp tableAlgebra.cpp inputOutput.cpp minimizers.cpp verykool.cpp -lgfortran -lmultinest_mpi -lCCfits -lcfitsio -ljsoncpp -lgmp -lCGAL
#This command assumes that the environment variable CPATH has been set according to the right path.
#E.g.: export LD_LIBRARY_PATH=/Users/users/gvernard/myLibraries/eigen/include/eigen3/


#The following compiling command is the same as above but it enables profiling
#g++ -std=c++11 -fopenmp -pg -o verykool fastell.o imagePlane.cpp sourcePlane.cpp massModels.cpp tableAlgebra.cpp parseInput.cpp minimizers.cpp verykool.cpp -lgfortran -lmultinest_mpi -lCCfits -ljsoncpp



#nvcc -D_DEBUG -lcufft -lgomp -arch=compute_11 --compiler-options "-fno-strict-aliasing -Wall -fopenmp" -O3 -o convolution fourierFuncs.cu initFuncs.cpp profileFuncs.cpp semiRandom.cpp outputFuncs.cu imageFuncs.cpp lcurveFuncs.cpp convolution.cpp

