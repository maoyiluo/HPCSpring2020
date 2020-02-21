#include <stdio.h>
#include "utils.h"
#include <cmath>

double residual(double*u, double*f, double h, int N){
    double res = 0;
    double coef = 1.0/(h*h);
    res += (coef*(2 * u[0] - u[1]) - f[0])*(coef*(2 * u[0] - u[1]) - f[0]);
    res += (coef*(2 * u[N-1] - u[N-2]) - f[N-1]) * (coef*(2 * u[N-1] - u[N-2]) - f[N-1]);
    for(int i = 1; i < N-1; i++){
        res += (coef*(-u[i-1] + 2*u[i] - u[i+1]) - f[i])*(coef*(-u[i-1] + 2*u[i] - u[i+1]) - f[i]);
    }
    return sqrt(res);
}

void jacobi(double* u, double* new_u, double* f, double h, int N){
    double coef = 1.0/(h*h);
    new_u[0] = 1.0/(coef*2)*(f[0] - coef*(-u[1]));
    new_u[N-1] = 1.0/(coef*2)*(f[N-1] - coef*(-u[N-2]));
    for(int i = 1; i<N-1; i++){
        new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-u[i-1] - u[i+1]));
    }
}

void gaussSeidel(double* u, double* new_u, double* f, double h, int N){
    double coef = 1.0/(h*h);
    new_u[0] = 1.0/(coef*2)*(f[0] - coef*(-u[1]));
    new_u[N-1] = 1.0/(coef*2)*(f[N-1] - coef*(-new_u[N-2]));
    for(int i = 1; i<N-1; i++){
        new_u[i] = 1.0/(coef*2)*(f[i] - coef*(-new_u[i-1] - u[i+1]));
    }
}

void iteration(int N, std::string method){
    double* new_u = (double*) malloc(N*sizeof(double)); //new_u is used to store
    double* u = (double*) malloc(N*sizeof(double));
    double* f = (double*) malloc(N*sizeof(double));
    double* temp;  //used for swapping u and new_u;
    double h = 1.0/(N+1);

    // Initialized f and u, 
    //let f to be a vector all equals to 1 
    //and u to be the vector all equals to 0.
    for(int i = 0; i<N; i++){
        u[i] = 0;
        f[i] = 1;  
        new_u[i] = 0;
    }

    double initial_residual = residual(u,f,h,N);
    double current_residual = residual(u,f,h,N);
    int iteration = 0;
    while(initial_residual/current_residual <=1e6 && iteration < 100){
        if(method == "jacobi")
            jacobi(u, new_u, f, h, N);
        else
            gaussSeidel(u, new_u, f, h, N);
        iteration++;
        temp = new_u;
        new_u = u;
        u = temp;
        current_residual = residual(u,f,h,N);
        //printf("res: %f, iter: %d, factor: %f\n", residual(u,f,h,N), iteration, initial_residual/current_residual);
    }

    free(u);
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

    Timer t;
    t.tic();
    for(int i = 0; i < 100; i++)
        iteration(N, "jacobi");
    double time = t.toc();
    printf("runtime:%f\n", time);

    t.tic();
    for(int i = 0; i < 100; i++)
        iteration(N, "gauss-seidel");
    time = t.toc();
    printf("runtime:%f\n", time);

    return 0;
}
