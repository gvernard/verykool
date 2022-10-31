ACLOCAL_AMFLAGS=-I m4

AM_CXXFLAGS = -std=c++11 -fopenmp -g -frounding-math





BUILT_SOURCES = libcommon.la

lib_LTLIBRARIES = libcommon.la

libcommon_la_SOURCES= \
	src/verykool/min_iterator.cpp \
	src/verykool/min_multinest.cpp \
	src/verykool/min_test.cpp \
	src/common/covKernels.cpp \
	src/common/eigenAlgebra.cpp \
	src/common/imagePlane.cpp \
	src/common/inputOutput.cpp \
	src/common/likelihoodModels.cpp \
	src/common/massModels.cpp \
	src/common/nonLinearPars.cpp \
	src/common/sourcePlane.cpp \
	src/common/sourceProfile.cpp \
	src/common/fastell.f
libcommon_la_CXXFLAGS = -std=c++17 -fPIC -g -frounding-math
libcommon_la_CPPFLAGS = -Iinclude
libcommon_la_LIBADD = -lcfitsio -lCCfits -ljsoncpp -lgsl -lgslcblas -lboost_system -lmpfr -lgmp -lCGAL -l:libgfortran.so.5
#libcommon_la_LIBADD = -lcfitsio -lCCfits -ljsoncpp -lgsl -lgslcblas -lboost_system -lmpfr -lgmp -lCGAL  -lgfortran


bin_PROGRAMS = \
	verykool_test \
	verykool


verykool_SOURCES = \
	src/verykool/verykool.cpp \
verykool_CXXFLAGS = -std=c++11 -fopenmp -g -frounding-math
verykool_CPPFLAGS = -Iinclude
verykool_LDADD = -lcommon -lmultinest_mpi 


verykool_test_SOURCES = \
	src/verykool/verykool_test.cpp
verykool_test_CXXFLAGS = -std=c++11 -g -frounding-math
verykool_test_CPPFLAGS = -Iinclude
verykool_test_LDADD = -lcommon


install-exec-hook:
	rm $(bin_PROGRAMS) $(lib_LTLIBRARIES)
