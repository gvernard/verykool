OBJ_DIR = obj
INC_DIR = inc
LIB_DIR = lib
BIN_DIR = bin

GPP = g++
MPICC = mpic++
FORTRAN = gfortran
FORTRAN_FLAGS = -fPIC


# Setting up verykool variables
VERYKOOL_SRC   = src/verykool
VERYKOOL_LIBS  = -lmultinest_mpi
VERYKOOL_FLAGS = -fopenmp -g -frounding-math

# Same for test version (no MPI)
VERYKOOL_SRC   = src/verykool
VERYKOOL_TEST_LIBS  = -lmultinest_mpi
VERYKOOL_TEST_FLAGS = -g -frounding-math

# Setting up cosmosis variables
COSMOSIS_SRC   = src/cosmosis
COSMOSIS_FLAGS = -fPIC -shared -g -frounding-math
COSMOSIS_LIBS  = -lcosmosis


# Setting up common variables
COMMON_SRC   = src/common
COMMON_FLAGS = -std=c++11 -fPIC -g -frounding-math
COMMON_LIBS  = -lgfortran -lCCfits -lcfitsio -ljsoncpp -lgmp -lCGAL
DEPS         = imagePlane.hpp inputOutput.hpp massModels.hpp sourcePlane.hpp eigenAlgebra.hpp nonLinearPars.hpp likelihoodModels.hpp covKernels.hpp
OBJ          = imagePlane.o   inputOutput.o   massModels.o   sourcePlane.o   eigenAlgebra.o   nonLinearPars.o   likelihoodModels.o   covKernels.o   fastell.o
COMMON_DEPS  = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
COMMON_OBJ   = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${COMMON_OBJ}])

MOCKS_DIR    = create_mock_lenses
FP_SRC  = initFuncs.cpp polyclip.cpp polygons.cpp sourceProfile.cpp
FP_DEPS = initFuncs.hpp polyclip.hpp polygons.hpp sourceProfile.hpp fitsHeader.hpp 
FP_OBJ  = initFuncs.o   polyclip.o   polygons.o   sourceProfile.o
FPROJECT_SRC  = $(patsubst %,$(MOCKS_DIR)/FProject/%,$(FP_SRC)) #Pad names with dir
FPROJECT_DEPS = $(patsubst %,$(MOCKS_DIR)/FProject/%,$(FP_DEPS)) #Pad names with dir
FPROJECT_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(FP_OBJ)) #Pad names with dir







all: verykool cosmosis other



# Compiling common code
$(OBJ_DIR)/%.o: $(COMMON_SRC)/%.cpp $(COMMON_DEPS)
	$(GPP) $(COMMON_FLAGS) -I inc -c -o $@ $<
$(OBJ_DIR)/%.o: $(COMMON_SRC)/%.f
	$(FORTRAN) $(FORTRAN_FLAGS) -c -o $@ $<


# Compiling verykool code
$(OBJ_DIR)/min_multinest.o: $(VERYKOOL_SRC)/min_multinest.cpp
	$(MPICC) $(COMMON_FLAGS) $(VERYKOOL_FLAGS) -I inc -c -o $@ $<
$(OBJ_DIR)/min_simplex.o: $(VERYKOOL_SRC)/min_simplex.cpp
	$(GPP) $(COMMON_FLAGS) $(VERYKOOL_FLAGS) -I inc -c -o $@ $<
$(OBJ_DIR)/min_test.o: $(VERYKOOL_SRC)/min_test.cpp
	$(GPP) $(COMMON_FLAGS) $(VERYKOOL_TEST_FLAGS) -I inc -c -o $@ $<
$(OBJ_DIR)/min_iterator.o: $(VERYKOOL_SRC)/min_iterator.cpp
	$(GPP) $(COMMON_FLAGS) $(VERYKOOL_TEST_FLAGS) -I inc -c -o $@ $<
$(OBJ_DIR)/verykool.o: $(VERYKOOL_SRC)/verykool.cpp
	$(MPICC) $(COMMON_FLAGS) $(VERYKOOL_FLAGS) -I inc -c -o $@ $<
$(OBJ_DIR)/verykool_test.o: $(VERYKOOL_SRC)/verykool_test.cpp
	$(GPP) $(COMMON_FLAGS) $(VERYKOOL_TEST_FLAGS) -I inc -c -o $@ $<


# Compiling comsosis code
$(OBJ_DIR)/min_cosmosis.o: $(COSMOSIS_SRC)/min_cosmosis.cpp
	$(GPP) $(COMMON_FLAGS) -I inc -c -o $@ $< 


# Compiling other code
$(OBJ_DIR)/createCosmosisValuesPriorsIni.o: src/other/createCosmosisValuesPriorsIni.cpp
	$(GPP) -std=c++11 -I inc -c -o $@ $< 
$(OBJ_DIR)/createEmceeStart.o: src/other/createEmceeStart.cpp
	$(GPP) -std=c++11 -I inc -c -o $@ $< 


# Compiling mocks code
$(OBJ_DIR)/addmachine.o: $(MOCKS_DIR)/addMachine/addmachine.cpp
	$(GPP) -std=c++11 -I inc -c -o $@ $<
$(OBJ_DIR)/automask.o: $(MOCKS_DIR)/autoMask/automask.cpp
	$(GPP) -std=c++11 -I inc -c -o $@ $<
$(OBJ_DIR)/kappa.o: $(MOCKS_DIR)/createPerturbations/kappa.cpp
	$(GPP) -std=c++11 -I inc -c -o $@ $<
$(OBJ_DIR)/initFuncs.o: $(MOCKS_DIR)/FProject/initFuncs.cpp
	$(GPP) -std=c++11 -c -o $@ $<
$(OBJ_DIR)/polyclip.o: $(MOCKS_DIR)/FProject/polyclip.cpp
	$(GPP) -std=c++11 -c -o $@ $<
$(OBJ_DIR)/polygons.o: $(MOCKS_DIR)/FProject/polygons.cpp
	$(GPP) -std=c++11 -I inc -c -o $@ $<
$(OBJ_DIR)/sourceProfile.o: $(MOCKS_DIR)/FProject/sourceProfile.cpp
	$(GPP) -std=c++11 -I inc -c -o $@ $<
$(OBJ_DIR)/fproject.o: $(MOCKS_DIR)/FProject/fproject.cpp $(FPROJECT_DEPS)
	$(GPP) -std=c++11 -I $(MOCKS_DIR)/FProject -I inc -c -o $@ $<






common: $(COMMON_OBJ)
	@echo ""


cosmosis: common $(OBJ_DIR)/min_cosmosis.o
	@echo ""
	$(CC) $(COSMOSIS_FLAGS) -I inc -o $(LIB_DIR)/libverykool.so $(COMMON_OBJ) $(OBJ_DIR)/min_cosmosis.o $(COMMON_LIBS) $(COSMOSIS_LIBS)
	@echo ""


verykool: common $(OBJ_DIR)/verykool.o $(OBJ_DIR)/min_multinest.o $(OBJ_DIR)/min_test.o $(OBJ_DIR)/min_iterator.o $(VERYKOOL_SRC)/minimizers.hpp verykool_test
	@echo ""
#	$(MPICC) $(VERYKOOL_FLAGS) -I inc -o $(BIN_DIR)/$@ $(COMMON_OBJ) $(OBJ_DIR)/verykool.o $(OBJ_DIR)/min_multinest.o $(OBJ_DIR)/min_simplex.o $(COMMON_LIBS) $(VERYKOOL_LIBS)
	$(MPICC) $(VERYKOOL_FLAGS) -I inc -o $(BIN_DIR)/$@ $(COMMON_OBJ) $(OBJ_DIR)/verykool.o $(OBJ_DIR)/min_multinest.o $(OBJ_DIR)/min_test.o $(OBJ_DIR)/min_iterator.o $(COMMON_LIBS) $(VERYKOOL_LIBS)
	@echo ""

verykool_test: common $(OBJ_DIR)/verykool_test.o $(OBJ_DIR)/min_multinest.o $(OBJ_DIR)/min_test.o $(OBJ_DIR)/min_iterator.o $(VERYKOOL_SRC)/minimizers.hpp
	@echo ""
	$(GPP) $(VERYKOOL_TEST_FLAGS) -I inc -o $(BIN_DIR)/$@ $(COMMON_OBJ) $(OBJ_DIR)/verykool_test.o $(OBJ_DIR)/min_multinest.o $(OBJ_DIR)/min_test.o $(OBJ_DIR)/min_iterator.o $(COMMON_LIBS) $(VERYKOOL_TEST_LIBS)
	@echo ""

other: common $(OBJ_DIR)/createCosmosisValuesPriorsIni.o $(OBJ_DIR)/createEmceeStart.o
	@echo ""
	$(GPP) -std=c++11 -I inc -o $(BIN_DIR)/createCosmosisValuesPriorsIni $(COMMON_OBJ) $(OBJ_DIR)/createCosmosisValuesPriorsIni.o $(COMMON_LIBS)
	$(GPP) -std=c++11 -I inc -o $(BIN_DIR)/createEmceeStart $(COMMON_OBJ) $(OBJ_DIR)/createEmceeStart.o $(COMMON_LIBS)
	@echo ""

mocks: common $(OBJ_DIR)/addmachine.o $(OBJ_DIR)/automask.o $(OBJ_DIR)/kappa.o $(OBJ_DIR)/fproject.o $(FPROJECT_OBJ)
	$(GPP) -std=c++11 -I inc -o $(MOCKS_DIR)/addMachine/addmachine $(OBJ_DIR)/addmachine.o $(OBJ_DIR)/imagePlane.o -ljsoncpp -lcfitsio -lCCfits -lfftw3
	$(GPP) -std=c++11 -I inc -o $(MOCKS_DIR)/autoMask/automask $(OBJ_DIR)/automask.o $(OBJ_DIR)/imagePlane.o -ljsoncpp -lcfitsio -lCCfits -lfftw3
	$(GPP) -std=c++11 -I inc -o $(MOCKS_DIR)/createPerturbations/kappa $(OBJ_DIR)/kappa.o $(OBJ_DIR)/imagePlane.o -ljsoncpp -lcfitsio -lCCfits
	$(GPP) -std=c++11 -I inc -I $(MOCKS_DIR)/FProject -o $(MOCKS_DIR)/FProject/fproject $(OBJ_DIR)/fproject.o $(OBJ_DIR)/fastell.o $(OBJ_DIR)/imagePlane.o $(OBJ_DIR)/nonLinearPars.o $(OBJ_DIR)/massModels.o $(OBJ_DIR)/sourcePlane.o $(FPROJECT_OBJ) -ljsoncpp -lcfitsio -lCCfits -lgfortran -lgmp -lCGAL

# src/common/source_plane.o
support: analysis_tools/create_gridded_source.cpp 
	$(GPP) -std=c++11 -I inc -o analysis_tools/create_gridded_source -I $(MOCKS_DIR)/FProject $(OBJ_DIR)/fastell.o $(OBJ_DIR)/sourcePlane.o $(OBJ_DIR)/massModels.o $(OBJ_DIR)/imagePlane.o $(OBJ_DIR)/nonLinearPars.o $(OBJ_DIR)/sourceProfile.o analysis_tools/create_gridded_source.cpp -ljsoncpp -lcfitsio -lCCfits -lgfortran -lgmp -lCGAL


clean:
	$(RM) -r $(OBJ_DIR)/* $(LIB_DIR)/*

clean_mocks:
	$(RM) $(MOCKS_DIR)/addMachine/addmachine $(MOCKS_DIR)/automask/automask $(MOCKS_DIR)/createPerturbations/kappa $(MOCKS_DIR)/FProject/fproject

