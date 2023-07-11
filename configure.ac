AC_INIT([verykool], [0.1], [gvernard@astro.rug.nl])
AC_LANG(C++)
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIRS([m4])
AC_PREFIX_DEFAULT($PWD)

MPIRUN=mpirun

AC_SUBST(VKL_HOME,$PWD)
AC_SUBST(MPIRUN,$MPIRUN)

AC_PROG_CXX([mpic++])
AC_PROG_F77
AC_PROG_MKDIR_P # MKDIR_P
AC_PROG_SED # SED

LT_INIT([disable-static])



AC_DEFUN([AC_INSTALLATION_INSTRUCTIONS],
	[
		name=$1
		echo ""
		echo "Program/library $(tput bold)$(tput setaf 1)${name} is missing$(tput sgr0), you need to install it using your package manager."
		echo "Alternatively, you can install ${name} from source using the following instructions:"
		echo $(tput setaf 6)
		cat ${VKL_HOME}/third_party/install_${name}.txt
		echo $(tput sgr0)
		echo "If you use the 'prefix' option while running ./configure above, then you need to pass:"
		echo "   $(tput setaf 2)--with-${name}=prefix$(tput sgr0)"
		echo "when running the VKL ./configure script."
		AC_MSG_ERROR("Program/library ${name} not found")
	]
)





ac_new_CPPFLAGS=""
ac_new_LDFLAGS=""
ac_new_LIBS=""

AC_ARG_WITH([cfitsio],
    [AS_HELP_STRING([--with-cfitsio=prefix],[give the path to the cfitsio library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lcfitsio"
    ],
    [
	AC_CHECK_LIB(cfitsio,main,[],[AC_INSTALLATION_INSTRUCTIONS(cfitsio)])
	AC_CHECK_HEADERS(fitsio.h,[],[AC_MSG_ERROR("you need to have the cfitsio header files installed.")])
    ])

AC_ARG_WITH([CCfits],
    [AS_HELP_STRING([--with-CCfits=prefix],[give the path to the CCfits library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lCCfits"
    ],
    [
	AC_CHECK_LIB(CCfits,main,[],[AC_INSTALLATION_INSTRUCTIONS(CCfits)])
	AC_CHECK_HEADERS(CCfits/CCfits.h,[],[AC_MSG_ERROR("you need to have the CCfits header files installed.")])
    ])

AC_ARG_WITH([jsoncpp],
    [AS_HELP_STRING([--with-jsoncpp=prefix],[give the path to the jsoncpp library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib64 -Wl,-rpath -Wl,${withval}/lib64"
	ac_new_LIBS+=" -ljsoncpp"
    ],
    [
	AC_CHECK_LIB(jsoncpp,main,[],[AC_INSTALLATION_INSTRUCTIONS(jsoncpp)])
	AC_CHECK_HEADERS(json/json.h,[],[AC_MSG_ERROR("you need to have the jsoncpp header files installed.")])
    ])

AC_ARG_WITH([jq],
    [AS_HELP_STRING([--with-jq=prefix],[give the path to the jq program binary file.])],
    [
	AC_SUBST(AM_JQ,${withval}/bin/jq)
    ],
    [
	AC_PATH_PROG(AM_JQ,jq,[no],[path=$PATH])
	if test "$AM_JQ" = no ; then
	   AC_INSTALLATION_INSTRUCTIONS(jq)
	fi
    ])

AC_ARG_WITH([boost],
    [AS_HELP_STRING([--with-boost=prefix],[give the path to the boost library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lboost_system"
    ],
    [
	AC_CHECK_LIB(boost_system,main,[],[AC_MSG_ERROR("you need to have the boost library installed.")])
	AC_CHECK_HEADERS(boost/limits.hpp,[],[AC_MSG_ERROR("you need to have the boost header files installed.")])
    ])

AC_ARG_WITH([mpfr],
    [AS_HELP_STRING([--with-mpfr=prefix],[give the path to the mpfr library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lmpfr"
    ],
    [
	AC_CHECK_LIB(mpfr,main,[],[AC_MSG_ERROR("you need to have the mpfr library installed.")])
	AC_CHECK_HEADERS(mpfr.h,[],[AC_MSG_ERROR("you need to have the mpfr header files installed.")])
    ])

AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=prefix],[give the path to the gmp library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lgmp"
    ],
    [
	AC_CHECK_LIB(gmp,main,[],[AC_MSG_ERROR("you need to have the gmp library installed.")])
	AC_CHECK_HEADERS(gmp.h,[],[AC_MSG_ERROR("you need to have the gmp header files installed.")])
    ])

AC_ARG_WITH([CGAL],
    [AS_HELP_STRING([--with-CGAL=prefix],[give the path to the CGAL library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib64 -Wl,-rpath -Wl,${withval}/lib64"
	ac_new_LIBS+=" -lCGAL"
    ],
    [
	AC_CHECK_LIB(CGAL,main,[],[AC_MSG_ERROR("you need to have the CGAL library installed.")])
	#AC_CHECK_HEADERS(CGAL/CGAL.h,[],[AC_MSG_ERROR("you need to have the CGAL header files installed.")])
    ])
    
AC_ARG_WITH([multinest],
    [AS_HELP_STRING([--with-multinest=prefix],[give the path to the multinest library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lmultinest_mpi"
    ],
    [
	AC_CHECK_LIB(multinest,main,[],[AC_MSG_ERROR("you need to have the multinest library installed.")])
	AC_CHECK_HEADERS(multinest.h,[],[AC_MSG_ERROR("you need to have the multinest header files installed.")])
    ])

AC_ARG_WITH([eigen],
    [AS_HELP_STRING([--with-eigen=prefix],[give the path to the eigen library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include/eigen3"
    ],
    [
	AC_CHECK_HEADERS([Eigen/Dense],[],[AC_MSG_ERROR("you need to have the Eigen header files installed.")])
    ])


AC_ARG_WITH([cuda],
   [  --with-cuda=PATH    prefix where cuda is installed [default=/usr/local/cuda]])

if test -n "$with_cuda"
then
   CUDA_CFLAGS="-I$with_cuda/include"
   CUDA_LIBS="-lcudart -lcusparse -lcusolver"
   CUDA_LDFLAGS="-L$with_cuda/lib64"
   NVCC="$with_cuda/bin/nvcc"
else
   CUDA_CFLAGS="-I/usr/local/cuda/include"
   CUDA_LIBS="-lcudart -lcusparse -lcusolver"
   CUDA_LDFLAGS="-L/usr/local/cuda/lib64"
   NVCC="nvcc"
fi

dnl echo $ac_new_CPPFLAGS
dnl echo $ac_new_LDFLAGS
dnl echo $ac_new_LIBS

CPPFLAGS+="$ac_new_CPPFLAGS"
LDFLAGS+="$ac_new_LDFLAGS"
LIBS+="$ac_new_LIBS"

dnl echo $CPPFLAGS
dnl echo $LDFLAGS
dnl echo $LIBS

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])


#CUDA_CFLAGS="-std=c++17 --tag=CXX --compiler-options '-fPIC' -Wno-deprecated-gpu-targets"
#CUDA_LIBS="-lcudart -lcusparse -lcusolver"
#NVCC="nvcc"

AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(CUDA_LDFLAGS)
AC_SUBST(NVCC)




AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([scripts/vkl_driver],[chmod +x scripts/vkl_driver])
AC_CONFIG_FILES([scripts/vkl_analyzer],[chmod +x scripts/vkl_analyzer])

AC_OUTPUT
