#ifndef __PCG_SOLVER_HPP__
#define __PCG_SOLVER_HPP__

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>


class PCGsolver
{
private:
    cublasHandle_t cubHandle;
    cusparseHandle_t cusHandle;
    cusparseMatDescr_t descr_A;
    cusparseMatDescr_t descr_L;

    // device data
    int N;
    int nonzero;
    double *d_ic=0; // nz
    double *d_x=0;  // N
    double *d_y=0;  // N
    double *d_z=0;  // N
    double *d_r=0;  // N
    double *d_rt=0; // N
    double *d_xt=0; // N
    double *d_q=0;  // N
    double *d_p=0;  // N
    double alpha;
    double beta;
    double alpha1 = 1.0;
    double beta0 = 0.0;
    double rTr;
    double pTq;
    double rho;
    double rhot;

    // device data size
    int d_N = 0;
    int d_nz = 0;

    // constant
    int max_iter;
    double tolerance;

public:
    PCGsolver(int max_iter=1000, double tol=1e-12);
    ~PCGsolver();

    void solve_gpumem(int N, int nz,
                    double *d_A, int *d_rowIdx, int *d_colIdx,
                    double *d_b, double *d_guess);

    void convert_coo2csr(const int N, const int nonzero, const int *cooRowIdx, int *csrRowIdx);

    double *get_device_x();

private:
    void check_and_resize();
    void free_memory();

};

#endif
