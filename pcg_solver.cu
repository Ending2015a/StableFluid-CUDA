
#include "pcg_solver.hpp"

#include "../common/error_helper.hpp"


PCGsolver::PCGsolver(int max_iter, double tol)  : max_iter(max_iter), tolerance(tol)
{
    // initialize cuBLAS & cuSPARSE
    error_check(cublasCreate(&cubHandle));
    error_check(cusparseCreate(&cusHandle));

    // create descriptor of matrix A
    error_check(cusparseCreateMatDescr(&descr_A));

    // initialize properties of matrix A
    error_check(cusparseSetMatType(descr_A, CUSPARSE_MATRIX_TYPE_SYMMETRIC));
    error_check(cusparseSetMatFillMode(descr_A, CUSPARSE_FILL_MODE_LOWER));
    error_check(cusparseSetMatDiagType(descr_A, CUSPARSE_DIAG_TYPE_NON_UNIT));
    error_check(cusparseSetMatIndexBase(descr_A, CUSPARSE_INDEX_BASE_ZERO));

    // create descriptor of matrix L
    error_check(cusparseCreateMatDescr(&descr_L));

    // initialize properties of matrix L
    error_check(cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_TRIANGULAR));
    error_check(cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO));
    error_check(cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER));
    error_check(cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_NON_UNIT));
}

PCGsolver::~PCGsolver()
{
    // free data
    free_memory();

    // cusparse
    cusparseDestroyMatDescr(descr_A);
    cusparseDestroyMatDescr(descr_L);
    cusparseDestroy(cusHandle);
    cublasDestroy(cubHandle);
}

void PCGsolver::solve_gpumem(int N, int nz,
                        double *d_A, int *d_rowIdx, int *d_colIdx,
                        double *d_b, double *d_guess)
{
    // check size
    this->N = N;
    this->nonzero = nz;
    free_memory();
    check_and_resize();

    // analyze matrix A (This will be used in incomplete-cholesky factorization)
    cusparseSolveAnalysisInfo_t info_A;
    error_check(cusparseCreateSolveAnalysisInfo(&info_A));
    error_check(cusparseDcsrsv_analysis(cusHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            N, nonzero, descr_A, d_A, d_rowIdx, d_colIdx, info_A));

    // copy matrix A
    error_check(cudaMemcpy(d_ic, d_A, nonzero * sizeof(double), cudaMemcpyDeviceToDevice));

    // compute IC factorization
    error_check(cusparseDcsric0(cusHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            N, descr_A, d_ic, d_rowIdx, d_colIdx, info_A));

    // analyze matrix L & U
    cusparseSolveAnalysisInfo_t info_L;
    error_check(cusparseCreateSolveAnalysisInfo(&info_L));
    error_check(cusparseDcsrsv_analysis(cusHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            N, nonzero, descr_L, d_ic, d_rowIdx, d_colIdx, info_L));

    cusparseSolveAnalysisInfo_t info_U;
    error_check(cusparseCreateSolveAnalysisInfo(&info_U));
    error_check(cusparseDcsrsv_analysis(cusHandle, CUSPARSE_OPERATION_TRANSPOSE,
                            N, nonzero, descr_L, d_ic, d_rowIdx, d_colIdx, info_U));

    // set initial guess
    if(d_guess == NULL)
    {
        error_check(cudaMemset(d_x, 0, N * sizeof(double)));
    }
    else
    {
        error_check(cudaMemcpy(d_x, d_guess, N * sizeof(double), cudaMemcpyDeviceToDevice));
    }
    error_check(cudaMemcpy(d_r, d_b, N * sizeof(double), cudaMemcpyDeviceToDevice));
    
    // solve
    int k;
    for(k=0;k<max_iter;++k)
    {
        error_check(cublasDnrm2(cubHandle, N, d_r, 1, &rTr));
        if(rTr < tolerance)
           break;
        
        error_check(cusparseDcsrsv_solve(cusHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            N, &alpha1, descr_L, d_ic, d_rowIdx, d_colIdx, info_L, d_r, d_y));

        error_check(cusparseDcsrsv_solve(cusHandle, CUSPARSE_OPERATION_TRANSPOSE,
                            N, &alpha1, descr_L, d_ic, d_rowIdx, d_colIdx, info_U, d_y, d_z));

        rhot = rho;  //store last rho
        error_check(cublasDdot(cubHandle, N, d_r, 1, d_z, 1, &rho));  //compute new rho
        
        if(k == 0)
        {
            error_check(cublasDcopy(cubHandle, N, d_z, 1, d_p, 1));
        }
        else
        {
            beta = rho/rhot;
            error_check(cublasDscal(cubHandle, N, &beta, d_p, 1));
            error_check(cublasDaxpy(cubHandle, N, &alpha1, d_z, 1, d_p, 1));
        }
        // Compute q <- Ap
        error_check(cusparseDcsrmv(cusHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        N, N, nonzero, &alpha1, descr_A, d_A, d_rowIdx, d_colIdx, d_p, &beta0, d_q));

        error_check(cublasDdot(cubHandle, N, d_p, 1, d_q, 1, &pTq));
        alpha = rho/pTq;
        error_check(cublasDaxpy(cubHandle, N, &alpha, d_p, 1, d_x, 1));
        double nalpha = -alpha;
        error_check(cublasDaxpy(cubHandle, N, &nalpha, d_q, 1, d_r, 1));
    }

    std::cout << "[PCGsolver] solved in " << k << " iterations, final norm(r) = " 
              << std::scientific << rTr << std::endl;

    error_check(cusparseDestroySolveAnalysisInfo(info_A));
    error_check(cusparseDestroySolveAnalysisInfo(info_L));
    error_check(cusparseDestroySolveAnalysisInfo(info_U));

}

void PCGsolver::convert_coo2csr(const int N, const int nonzero, const int* cooRowIdx, int *csrRowIdx)
{
    error_check(cusparseXcoo2csr(cusHandle, cooRowIdx, nonzero, N, csrRowIdx, CUSPARSE_INDEX_BASE_ZERO));
}

double *PCGsolver::get_device_x()
{
    return d_x;
}

void PCGsolver::free_memory()
{
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_z);
    cudaFree(d_r);
    cudaFree(d_rt);
    cudaFree(d_xt);
    cudaFree(d_q);
    cudaFree(d_p);
    cudaFree(d_ic);
}

void PCGsolver::check_and_resize()
{
    error_check(cudaMalloc(&d_x, N * sizeof(double)));
    error_check(cudaMalloc(&d_y, N * sizeof(double)));
    error_check(cudaMalloc(&d_z, N * sizeof(double)));
    error_check(cudaMalloc(&d_r, N * sizeof(double)));
    error_check(cudaMalloc(&d_rt, N * sizeof(double)));
    error_check(cudaMalloc(&d_xt, N * sizeof(double)));
    error_check(cudaMalloc(&d_q, N * sizeof(double)));
    error_check(cudaMalloc(&d_p, N * sizeof(double)));
    d_N = N;
    error_check(cudaMalloc(&d_ic, nonzero * sizeof(double)));
    d_nz = nonzero;
}
