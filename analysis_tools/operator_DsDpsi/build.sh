#!/bin/bash
OBJ_DIR=../../obj

g++ -std=c++11 -g -frounding-math -I ../../inc -o create_operator $OBJ_DIR/inputOutput.o $OBJ_DIR/min_multinest.o $OBJ_DIR/min_test.o $OBJ_DIR/min_iterator.o $OBJ_DIR/fastell.o $OBJ_DIR/sourcePlane.o $OBJ_DIR/massModels.o $OBJ_DIR/imagePlane.o $OBJ_DIR/nonLinearPars.o $OBJ_DIR/sourceProfile.o $OBJ_DIR/eigenAlgebra.o $OBJ_DIR/likelihoodModels.o $OBJ_DIR/covKernels.o create_operator.cpp -ljsoncpp -lcfitsio -lCCfits -lgfortran -lgmp -lCGAL -lgsl -lgslcblas -lmultinest_mpi
