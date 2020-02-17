#include <stdio.h>
#include "utils.h"
#include <cmath>

double residual(double*u, double*f, double h, int N){
    double res = 0;
    double coef = 1.0/(h*h);
    for(int i = 0; i < N; i++){
        if(i == 0)
            res += (coef*(2 * u[0] - u[1]) - f[0])*(coef*(2 * u[0] - u[1]) - f[0]);
        else if(i == N-1)
            res += (coef*(2 * u[N-1] - u[N-2]) - f[N-1]) * (coef*(2 * u[N-1] - u[N-2]) - f[N-1]);
        else
        {
            res += (coef*(-u[i-1] + 2*u[i] - u[i+1]) - f[i])*(coef*(-u[i-1] + 2*u[i] - u[i+1]) - f[i]);
        }
    }
    return sqrt(res);
}

void jacobi(double* u, double* new_u, double* f, double h, int N){
    double coef = 1.0/(h*h);
    for(int i = 0; i<N; i++){
        if(i == 0)
            new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-u[1]));
        else if(i == N-1)
            new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-u[N-2]));
        else
            new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-u[i-1] - u[i+1]));
    }
}

void gaussSeidel(double* u, double* new_u, double* f, double h, int N){
    double coef = 1.0/(h*h);
    for(int i = 0; i<N; i++){
    if(i == 0)
        new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-u[1]));
    else if(i == N-1)
        new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-new_u[N-2]));
    else
        new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-new_u[i-1] - u[i+1]));
    }
}

void iteration(int N, std::string method, double* u){
    double* new_u = (double*) malloc(N*sizeof(double));
    double* temp;
    double* f = (double*) malloc(N*sizeof(double));
    double h = 1.0/(N+1);

    for(int i = 0; i<N; i++){
        u[i] = 0;
        f[i] = -1;
        new_u[i] = 0;
    }

    double initial_residual = residual(u,f,h,N);
    double current_residual = residual(u,f,h,N);
    int iteration = 0;
    while(initial_residual/current_residual >=1e-6 && iteration < 5000){
        if(method == "jacobi")
            jacobi(u, new_u, f, h, N);
        else
            gaussSeidel(u, new_u, f, h, N);
        iteration++;
        temp = new_u;
        new_u = u;
        u = temp;
        current_residual = residual(u,f,h,N);
        // printf("res: %f, iter: %d\n", residual(u,f,h,N), iteration);
    }

    free(new_u);
    free(f);
}

void printU(double* u, int N){
    for(int i = 0; i < N; i++){
        printf("%f",u[i]);
    }
    printf("\n");
}

int main(int argc, char** argv) {
    int N = strtol(argv[1], NULL, 10);
    double* u = (double*) malloc(N*sizeof(double));
    iteration(N, "jacobi", u);
    // printU(u,N);
   
    iteration(N, "gauss-seidel", u);
    // printU(u,N);
    return 0;
}
