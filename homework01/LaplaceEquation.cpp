#include <stdio.h>
#include "utils.h"

#define N 100

double residual(double*u, double*f){
    double res = 0;
    for(int i = 0; i <N; i++){
        if(i == 0)
            res += (2*u[0] - 2*u[1] - f[0])*(2*u[0] - 2*u[1] - f[0]);
        else if(i == N-1)
            res += (2*u[N-1] - 2*u[N-2] - f[N-1])*(2*u[N-1] - 2*u[N-2] - f[N-1]);
        else
        {
            res += (-u[i-1] + 2*u[i] - u[i-1] - f[i])*(-u[i-1] + 2*u[i] - u[i-1] - f[i]);
        }
    }
    return res;
}

int main(int argc, char** argv) {
  double* u = (double*) malloc(N*sizeof(double));
  double* f = (double*) malloc(N*sizeof(double));

  for(int i = 0; i<N; i++){
      u[i] = 2;
      f[i] = 1;
  }

  std::cout<<residual(u,f);

  return 0;
}
