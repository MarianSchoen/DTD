AC_INIT([DTD], 1.0.0) dnl package name, version
AC_PREREQ(2.60)

# TODO could also check for preinstalled eigen...
AC_CHECK_PROG(found_git, git, yes, no)
if test [ x"${found_git}" == "xyes" ]; then
   AC_PATH_PROG(git, git)
   AC_MSG_NOTICE(checking out / updating eigen...)
   ${git} submodule update --init
   DTD_CPPFLAGS="$DTD_CPPFLAGS -I./eigen/"
fi

AC_LANG(C++)

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP

# reuse blas from the R installation (and assume that is ok)
BLAS_LIB=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
if test -n "${BLAS_LIB}"; then
   DTD_CPPFLAGS="$DTD_CPPFLAGS -DEIGEN_USE_BLAS"
   LIBS="$LIBS $BLAS_LIB"
fi

# check for openmp
AC_OPENMP
DTD_CPPFLAGS="$DTD_CPPFLAGS $OPENMP_CXXFLAGS -Wall"
LIBS="$LIBS $OPENMP_CXXFLAGS"

dnl substitute RODBC_CPPFLAGS and LIBS
AC_SUBST(DTD_CPPFLAGS)
AC_SUBST(LIBS)
dnl and do substitution in the src/Makevars.in and src/config.h
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
