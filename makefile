CC = mpic++
CFLAGS = -std=c++11 -fopenmp
FORTRAN = gfortran
FORTRAN_FLAGS = 

ODIR = obj
LIBS = -lgfortran -lmultinest_mpi -lCCfits -lcfitsio -ljsoncpp -lgmp -lCGAL -lMinuit2
DEPS = imagePlane.hpp inputOutput.hpp massModels.hpp sourcePlane.hpp tableAlgebra.hpp nonLinearPars.hpp minimizers.hpp
OBJ  = imagePlane.o   inputOutput.o   massModels.o   sourcePlane.o   tableAlgebra.o   nonLinearPars.o   mn_multinest.o   mn_simplex.o  verykool.o fastell.o

#Pad object names with object dir
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.f
	$(FORTRAN) -c -o $@ $< $(FORTRAN_FLAGS)

verykool: $(_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)


