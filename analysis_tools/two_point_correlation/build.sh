#!/bin/bash
OBJ_DIR=../../obj

g++ -std=c++11 -g -frounding-math -I ../../inc -o theo_matrix_correlation -I ../../create_mock_lenses/FProject $OBJ_DIR/fastell.o $OBJ_DIR/sourcePlane.o $OBJ_DIR/massModels.o $OBJ_DIR/imagePlane.o $OBJ_DIR/nonLinearPars.o $OBJ_DIR/sourceProfile.o $OBJ_DIR/eigenAlgebra.o $OBJ_DIR/likelihoodModels.o $OBJ_DIR/covKernels.o theo_matrix_correlation.cpp -ljsoncpp -lcfitsio -lCCfits -lgfortran -lgmp -lCGAL

# imagePlane.o, massModels.o, and fastell.o are needed by sourcePlane.o (no need to include the headers in main)
g++ -std=c++11 -g -frounding-math -I ../../inc -o reconstruction_matrix_correlation -I ../../create_mock_lenses/FProject $OBJ_DIR/sourcePlane.o $OBJ_DIR/sourceProfile.o $OBJ_DIR/imagePlane.o $OBJ_DIR/massModels.o $OBJ_DIR/fastell.o reconstruction_matrix_correlation.cpp -lcfitsio -lCCfits -lgfortran -lgmp -lCGAL
