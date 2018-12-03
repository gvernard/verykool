#!/bin/bash -f
SRC_PATH_VERYKOOL=/Users/users/gvernard/myGit/VeryKooL/src/common
SRC_PATH_FPROJECT=/Users/users/gvernard/myCodes/createMockStrongLenses/FProject
INC_PATH_VERYKOOL=/Users/users/gvernard/myGit/VeryKooL/inc
INC_PATH_FPROJECT=/Users/users/gvernard/myCodes/createMockStrongLenses/FProject

# fastell.o / imagePlane.cpp / sourcePlane.cpp / massModels.cpp / nonLinearPars.cpp are only for the code to compile

g++ -std=c++11 -o create_gridded_source /Users/users/gvernard/myGit/VeryKooL/obj/fastell.o -I$INC_PATH_VERYKOOL/ -I$INC_PATH_FPROJECT/ $SRC_PATH_VERYKOOL/sourcePlane.cpp $SRC_PATH_VERYKOOL/massModels.cpp $SRC_PATH_VERYKOOL/imagePlane.cpp $SRC_PATH_VERYKOOL/nonLinearPars.cpp $SRC_PATH_FPROJECT/sourceProfile.cpp create_gridded_source.cpp -ljsoncpp -lCCfits -lcfitsio -lgfortran -lgmp -lCGAL
