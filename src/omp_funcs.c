#include <omp.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/**   This function is a slightly modified version of
   function setCores in package spMC version 0.2.2
   written by Luca Sartore <drwolf85@gmail.com>
*/

void getNumCores(int *n) {
  /* Get the max number of CPU cores
          *n - num of cores */
  #ifdef _OPENMP
    *n = omp_get_num_procs();
  #else
    *n = 1;
  #endif
}

void getNumThreads(int *n) {
  /* Get the number of threads to use
          *n - num of threads */
  #ifdef _OPENMP
    #pragma omp parallel default(shared)
    {
      #pragma omp master
        *n = omp_get_num_threads();
    }
  #else
    *n = 1;
  #endif
}

void setNumThreads(int *n) {
  /* Set the number of threads to use
          *n - num of threads */
  #ifdef _OPENMP
    if (omp_get_num_procs() < *n) {
      *n = omp_get_num_procs();
      omp_set_num_threads(*n);
    }
    else if(*n > 0) {
      omp_set_num_threads(*n);
    }
    else {
      omp_set_num_threads(1);
      *n = 1;
    }
  #else
    *n = 0;
  #endif
}
