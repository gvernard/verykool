AC_INIT([verykool], [0.1], [gvernard@astro.rug.nl])
AC_LANG(C++)
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIRS([m4])
AC_PREFIX_DEFAULT($PWD)

AC_SUBST(VKL_HOME,$PWD)

CXX=mpic++

AC_PROG_CXX
AC_PROG_F77

LT_INIT([disable-static])



AC_DEFUN([AC_INSTALLATION_INSTRUCTIONS],
	[
		name=$1
		echo ""
		echo "Program/library $(tput bold)$(tput setaf 1)${name} is missing$(tput sgr0), you need to install it using your package manager."
		echo "Alternatively, you can install ${name} from source using the following instructions:"
		echo $(tput setaf 6)
		cat ${MOLET_HOME}/third_party/install_${name}.txt
		echo $(tput sgr0)
		echo "If you use the 'prefix' option while running ./configure above, then you need to pass:"
		echo "   $(tput setaf 2)--with-${name}=prefix$(tput sgr0)"
		echo "when running the MOLET ./configure script."
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

AC_ARG_WITH([boost],
    [AS_HELP_STRING([--with-boost=prefix],[give the path to the boost library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lboost_system"
    ],
    [
	AC_CHECK_LIB(boost,main,[],[AC_MSG_ERROR("you need to have the boost library installed.")])
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
	ac_new_LIBS+=" -lmultinest"
    ],
    [
	AC_CHECK_LIB(multinest,main,[],[AC_MSG_ERROR("you need to have the multinest library installed.")])
	AC_CHECK_HEADERS(multinest.h,[],[AC_MSG_ERROR("you need to have the multinest header files installed.")])
    ])

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


AC_CONFIG_FILES([Makefile])

AC_OUTPUT
