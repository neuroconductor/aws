AC_INIT([aws], 1.9-0)


# find R home and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`;
CFLAGS=`"${RBIN}" CMD config CFLAGS`
FFLAGS=`"${RBIN}" CMD config FFLAGS`

# Checks for programs.
AC_PROG_CC
# Check for OpenMP
AC_OPENMP

OPENMP_FFLAGS=""
HAS_OPENMP=FALSE
# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp=no
if test -n "${OPENMP_CFLAGS}"; then
  AC_MSG_CHECKING([whether OpenMP will work in a package])
  AC_LANG_CONFTEST(
  [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
  PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
  AC_MSG_RESULT([${ac_pkg_openmp}])
fi

# if ${ac_pkg_openmp} = "yes" then we have OMP, otherwise it will be "no"
if test "${ac_pkg_openmp}" = no; then
  OPENMP_FFLAGS=""
  OPENMP_CFLAGS=""
  LGOMP_CFLAGS=""
else
  OPENMP_FFLAGS="-fopenmp"
  LGOMP_CFLAGS="-lgomp -lpthread"
  HAS_OPENMP=TRUE
fi

AC_SUBST(OPENMP_FFLAGS)
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(LGOMP_CFLAGS)
AC_SUBST(HAS_OPENMP)

AC_CONFIG_FILES([R/onLoad.r src/Makevars])
AC_OUTPUT

