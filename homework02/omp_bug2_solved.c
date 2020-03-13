/******************************************************************************
* FILE: omp_bug2.c
* DESCRIPTION:
*   Another OpenMP program with a bug. 
* AUTHOR: Blaise Barney 
* LAST REVISED: 04/06/05 
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
int nthreads, tid;
double total;  //total should be double rather than float to avoid overflow.

/*** Spawn parallel region ***/
#pragma omp parallel private(tid) //tid should be a private variable between each thread.
  {
  /* Obtain thread number */
  tid = omp_get_thread_num();
  /* Only master thread does this */
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d is starting...\n",tid);
  
  total = 0;
  #pragma omp barrier

  /* do some work */
  #pragma omp for schedule(static,10) reduction(+: total)  //use reduction clause.
  for (long i=0; i<1000000; i++){
      total += i*1.0;
  }

  printf ("Thread %d is done! Total= %f\n",tid,total);

  } /*** End of parallel region ***/
}
