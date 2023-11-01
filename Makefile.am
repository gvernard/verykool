ACLOCAL_AMFLAGS=-I m4

AM_CXXFLAGS = -std=c++11 -fopenmp -g -frounding-math





#BUILT_SOURCES = libgpu.la libcommon_gpu.la libcpu.la libcommon.la 
#lib_LTLIBRARIES = libgpu.la libcommon_gpu.la libcpu.la libcommon.la 

BUILT_SOURCES = libcpu.la libcommon.la 
lib_LTLIBRARIES = libcpu.la libcommon.la 




# Functions for the CPU
libcpu_la_SOURCES= \
	src/cpu/baseAlgebra.cpp
libcpu_la_CXXFLAGS = -std=c++17 -fPIC -g -frounding-math
libcpu_la_CPPFLAGS = -Iinclude




# # Functions for the GPU
# .cu.o:
# 	$(NVCC) -Isrc/gpu -Iinclude -o $@ -c $<

# .cu.lo:
# 	python $(top_srcdir)/cudalt.py $@ $(NVCC) --compiler-options=\" -Isrc/gpu -Iinclude $(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" -c $<

# LINK = $(LIBTOOL) --mode=link $(CXX) -o $@ $(CUDA_LDFLAGS) $(LDFLAGS) $(CUDA_LIBS) $(LIBS)
# libgpu_la_SOURCES= \
# 	src/gpu/baseAlgebra.cpp \
# 	src/gpu/cusparse_cholesky_solver.cu \
# 	src/gpu/cusparse_wrapper.cu
# libgpu_la_CXXFLAGS = -std=c++17 -fPIC -g -frounding-math
# libgpu_la_CPPFLAGS = -Iinclude




libcommon_la_SOURCES= \
	src/common/covKernels.cpp \
	src/common/imagePlane.cpp \
	src/common/likelihoodModels.cpp \
	src/common/eigenAlgebra.cpp \
	src/common/massModels.cpp \
	src/common/nonLinearPars.cpp \
	src/common/sourcePlane.cpp \
	src/common/sourceProfile.cpp \
	src/common/fastell.f
libcommon_la_CXXFLAGS = -std=c++17 -fPIC -g -frounding-math
libcommon_la_CPPFLAGS = -Iinclude
libcommon_la_LIBADD = -lcpu -lcfitsio -lCCfits -ljsoncpp -lgsl -lgslcblas -lboost_system -lmpfr -lgmp -lCGAL -l:libgfortran.so.5




# libcommon_gpu_la_SOURCES= \
# 	src/verykool/min_iterator.cpp \
# 	src/verykool/min_multinest.cpp \
# 	src/verykool/min_test.cpp \
# 	src/common/covKernels.cpp \
# 	src/common/imagePlane.cpp \
# 	src/common/inputOutput.cpp \
# 	src/common/likelihoodModels.cpp \
# 	src/common/eigenAlgebra.cpp \
# 	src/common/massModels.cpp \
# 	src/common/nonLinearPars.cpp \
# 	src/common/sourcePlane.cpp \
# 	src/common/sourceProfile.cpp \
# 	src/common/fastell.f
# libcommon_gpu_la_CXXFLAGS = -std=c++17 -fPIC -g -frounding-math
# libcommon_gpu_la_CPPFLAGS = -Iinclude
# libcommon_gpu_la_LIBADD = $(CUDA_LIBS) -lgpu -lcfitsio -lCCfits -ljsoncpp -lgsl -lgslcblas -lboost_system -lmpfr -lgmp -lCGAL -l:libgfortran.so.5




bin_PROGRAMS = \
	verykool \
	coolest_conversion


verykool_SOURCES = \
	src/verykool/inputOutput.cpp \
	src/verykool/min_test.cpp \
	src/verykool/min_multinest.cpp \
	src/verykool/min_iterator.cpp \
	src/verykool/verykool.cpp
verykool_CXXFLAGS = -std=c++11 -fopenmp -g -frounding-math
verykool_CPPFLAGS = -Iinclude
verykool_LDADD = -lmultinest_mpi -lcommon -lgsl


coolest_conversion_SOURCES = \
	analysis_tools/coolest_conversion.cpp
coolest_conversion_CXXFLAGS = -std=c++17 -g -frounding-math
coolest_conversion_CPPFLAGS = -Iinclude
coolest_conversion_LDADD = -lcommon




bin_SCRIPTS = \
	scripts/vkl_driver \
	scripts/vkl_analyzer
CLEANFILES = $(bin_SCRIPTS)




install-exec-hook:
	rm $(bin_PROGRAMS) $(bin_SCRIPTS)
