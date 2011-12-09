AC_INIT([aws], 1.7-1)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

FFLAGS=`"${R_HOME}/bin/R" CMD config FFLAGS`
FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`

AC_ARG_ENABLE([openmp],
	AC_HELP_STRING([--enable-openmp],
		[specify whether OpenMP is enabled or not]),
	[USE_OPENMP="$enableval"])

HAS_OPENMP=FALSE
if test "${USE_OPENMP}" = yes ; then
  AC_CHECK_LIB(gomp, omp_set_num_threads,
                [FFLAGS="${FFLAGS} -fopenmp";
                 FLIBS="${FLIBS} -lgomp";
                 HAS_OPENMP=TRUE],
                [AC_MSG_ERROR("Cannot find OpenMP library.")])
fi

AC_SUBST(FLIBS)
AC_SUBST(FFLAGS)
AC_SUBST(HAS_OPENMP)

AC_CONFIG_FILES([R/onLoad.r src/Makevars])
AC_OUTPUT

