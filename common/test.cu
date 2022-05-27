#include <iostream>
#include <cuda_runtime.h>
#include "error_helper.hpp"

int main()
{
    // error: invalid value
    double *a = 0, *b = 0;
    error_check(cudaMemcpy(a, b, 10, cudaMemcpyHostToDevice));

    return 0;
}
