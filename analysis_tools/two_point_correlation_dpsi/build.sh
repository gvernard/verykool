#!/bin/bash
OBJ_DIR=../../obj

g++ -w -std=c++11 -g -frounding-math -I ../../inc -o theo_matrix_correlation_dpsi $OBJ_DIR/sourcePlane.o $OBJ_DIR/massModels.o $OBJ_DIR/fastell.o $OBJ_DIR/imagePlane.o $OBJ_DIR/nonLinearPars.o $OBJ_DIR/eigenAlgebra.o $OBJ_DIR/covKernels.o theo_matrix_correlation_dpsi.cpp -ljsoncpp -lgfortran -lcfitsio -lCCfits -lgmp -lCGAL
