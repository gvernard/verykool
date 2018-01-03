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
VERYKOOL_LIBS  = -lmultinest_mpi -lsimplex -lMinuit2 
VERYKOOL_FLAGS = -fopenmp

# Setting up cosmosis variables
COSMOSIS_SRC   = src/cosmosis
COSMOSIS_FLAGS = -fPIC -shared
COSMOSIS_LIBS  = -lcosmosis


# Setting up common variables
COMMON_SRC   = src/common
COMMON_FLAGS = -std=c++11 -fPIC
COMMON_LIBS  = -lgfortran -lCCfits -lcfitsio -ljsoncpp -lgmp -lCGAL
DEPS         = imagePlane.hpp inputOutput.hpp massModels.hpp sourcePlane.hpp tableAlgebra.hpp nonLinearPars.hpp cov_kernels.hpp mainLogLike.hpp
OBJ          = imagePlane.o   inputOutput.o   massModels.o   sourcePlane.o   tableAlgebra.o   nonLinearPars.o   cov_kernels.o   fastell.o
COMMON_DEPS  = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
COMMON_OBJ   = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir



all: verykool cosmosis



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
$(OBJ_DIR)/verykool.o: $(VERYKOOL_SRC)/verykool.cpp
	$(MPICC) $(COMMON_FLAGS) $(VERYKOOL_FLAGS) -I inc -c -o $@ $<


# Compiling comsosis code
$(OBJ_DIR)/min_cosmosis.o: $(COSMOSIS_SRC)/min_cosmosis.cpp
	$(GPP) $(COMMON_FLAGS) -I inc -c -o $@ $< 



common: $(COMMON_OBJ)
	@echo ""
	@echo " >>> making common code: done!"
	@echo ""
	@echo ""


cosmosis: common $(OBJ_DIR)/min_cosmosis.o
	$(CC) $(COSMOSIS_FLAGS) -I inc -o $(LIB_DIR)/libverykool.so $(COMMON_OBJ) $(OBJ_DIR)/min_cosmosis.o $(COMMON_LIBS) $(COSMOSIS_LIBS)
	@echo ""
	@echo " >>> making cosmosis: done!"
	@echo ""
	@echo ""


verykool: common $(OBJ_DIR)/verykool.o $(OBJ_DIR)/min_multinest.o $(OBJ_DIR)/min_simplex.o $(VERYKOOL_SRC)/minimizers.hpp
	$(MPICC) $(VERYKOOL_FLAGS) -I inc -o $(BIN_DIR)/$@ $(COMMON_OBJ) $(OBJ_DIR)/verykool.o $(OBJ_DIR)/min_multinest.o $(OBJ_DIR)/min_simplex.o $(COMMON_LIBS) $(VERYKOOL_LIBS)
	@echo ""
	@echo " >>> making verykool: done!"
	@echo ""
	@echo ""


clean:
	$(RM) -r $(OBJ_DIR)/* $(LIB_DIR)/*


