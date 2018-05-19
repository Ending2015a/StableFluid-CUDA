#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include <cstring>
#include <cassert>

#include <cuda.h>
#include "pcg_solver.hpp"
#include "../common/error_helper.hpp"

// ===== Constant =====
int nx;
int ny;
double dx;
double dy;

int num;
int steps;
double dt;
double nu;
double rho;
const double g = 9.81;
const int block_size = 16;

int max_iter;
double tol;

__device__ __constant__ int d_nx;
__device__ __constant__ int d_ny;
__device__ __constant__ double d_dx;
__device__ __constant__ double d_dy;
__device__ __constant__ double d_nu;
__device__ __constant__ double d_rho;
__device__ __constant__ double d_g;

// ===== MAC grid =====
double *markers_x;
double *markers_y;

double *d_mx;  // x-coordinate of markers
double *d_my;  // y-coordinate of markers

double *d_u;  // horizontal component of velocity
double *d_v;  // vertical component of velocity
double *d_p;  // pressure
int *d_s;  // status

double *d_bu;  // backup u
double *d_bv;  // backup v

double *d_A;
double *d_rowIdx;
double *d_colIdx;
double *d_b;

// ===== kernel =====

// safe get from 2D matrix
template<typename T>
__device__ inline T _safe_get(const T* const field, const int x, const int y, 
                                const int lim_x, const int lim_y)
{
    if(x < 0 || x >= lim_x || y < 0 || y >= lim_y)
        return 0;
    return field[y*lim_x + x];
}

// safe set to 2D matrix
template<typename T>
__device__ inline void _safe_set(T* const field, const int x, const int y,
                                const int lim_x, const int lim_y, const T value)
{
    if(x < 0 || x >= lim_x || y < 0 || y >= lim_y)
        return;
    field[y*lim_x + x] = value;
}

// bilinear interpolation
__device__ inline double _interpolate(const double* const field, const double x, const double y,
                                     const int lim_x, const int lim_y)
{
    // global coord to index coord
    double idx_x = x/d_dx;
    double idx_y = y/d_dy;
    
    // integer part
    double flx = floor(idx_x);
    double fly = floor(idx_y);

    // lower
    int x0 = (int)flx;
    int y0 = (int)fly;

    //fractional part
    double fcx = idx_x - flx;
    double fcy = idx_y - fly;

    // upper
    int x1 = x0+1;
    int y1 = y0+1;

    // interp x
    double fy1 = (1.-fcx) * _safe_get(field, x0, y0, lim_x, lim_y) 
                    + fcx * _safe_get(field, x1, y0, lim_x, lim_y);
    double fy2 = (1.-fcx) * _safe_get(field, x0, y1, lim_x, lim_y)
                    + fcx * _safe_get(field, x1, y1, lim_x, lim_y);
    // interp y
    double f = (1.-fcy) * fy1 + fcy * fy2;

    return f;
}

// 4th-order Rung-Jutta (ODE solver)
__device__ inline void _RK4(const double* const u, const double* const v, const double dt,
                            double XG, double YG, double *XP, double *YP, 
                            const int lim_ux, const int lim_uy, const int lim_vx, const int lim_vy)
{
    // half index (global coord to field coord)
    const double hx = XG - d_dx/2.;
    const double hy = YG - d_dy/2.;

    // Runge-Kutta 4
    const double xk1 = dt * _interpolate(u, XG, hy, lim_ux, lim_uy);
    const double yk1 = dt * _interpolate(v, hx, YG, lim_vx, lim_vy);

    const double xk2 = dt * _interpolate(u, XG + xk1/2., hy + yk1/2., lim_ux, lim_uy);
    const double yk2 = dt * _interpolate(v, hx + xk1/2., YG + yk1/2., lim_vx, lim_vy);

    const double xk3 = dt * _interpolate(u, XG + xk2/2., hy + yk2/2., lim_ux, lim_uy);
    const double yk3 = dt * _interpolate(v, hx + xk2/2., YG + yk2/2., lim_vx, lim_vy);
    
    const double xk4 = dt * _interpolate(u, XG + xk3, hy + yk3, lim_ux, lim_uy);
    const double yk4 = dt * _interpolate(v, hx + xk3, YG + yk3, lim_vx, lim_vy);

    // advect 
    *XP = XG + (xk1 + 2.*xk2 + 2.*xk3 + xk4)/6.;
    *YP = YG + (yk1 + 2.*yk2 + 2.*yk3 + yk4)/6.;
}

// Update Grid Status (kernel)
template<int block_size>
__global__ void _updatestatus(const int num, const double* const x, const double* const y, int* const st)
{
    // thread index
    const int idx = block_size * blockIdx.x + threadIdx.x;

    // index coord
    const int X = (int)floor(x[idx]/d_dx);
    const int Y = (int)floor(y[idx]/d_dy);

    // i dont care about the race condition anyway...
    const int v = 1;
    _safe_set(st, X, Y, d_nx, d_ny, v);
}

// Advect Velocity (kernel)
template<int block_size>
__global__ void _advect(double* const u, const double* const bu, 
                    const double* const field_u, const double* const field_v, const double dt, 
                    const double off_x, const double off_y, const int lim_x, const int lim_y)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y*lim_x+x;

    // out of boundaries
    if(y >= lim_y || x >= lim_x)return;

    // index coord to global coord
    const double XG = ((double)x + off_x) * d_dx;
    const double YG = ((double)y + off_y) * d_dy;
    double XP;
    double YP;

    // trace backwards
    _RK4(field_u, field_v, -dt, XG, YG, &XP, &YP, d_nx+1, d_ny, d_nx, d_ny+1);

    // update 
    u[idx] = _interpolate(bu, XP-off_x*d_dx, YP-off_y*d_dy, lim_x, lim_y);
}


// Add Force (kernel)
template<int block_size>
__global__ void _addforce(double* const v, const double dt)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y*d_nx + x;

    // out of boundaries
    if(y >= d_ny+1 || x >= d_nx)return;

    // add gravity
    v[idx] = v[idx] - dt * d_g;
}


// Advect markers (kernel)
template<int block_size>
__global__ void _advectmarkers(const int num, double* const x, double* const y,
                            const double* const field_u, const double* const field_v, 
                            const double dt, const double lim_x, const double lim_y)
{
    // thread index
    const int idx = block_size * blockIdx.x + threadIdx.x;

    // out of bounds
    if(idx >= num) return;

    // get current marker position
    const double XP = x[idx];
    const double YP = y[idx];

    double XP1;
    double YP1;

    // trace forwards
    _RK4(field_u, field_v, dt, XP, YP, &XP1, &YP1, d_nx+1, d_ny, d_nx, d_ny+1);

    // clamp into boundaries & set new position
    XP1 = XP1 < 0. ? 0.: XP1 > lim_x ? lim_x : XP1;
    x[idx] = XP1;
    
    YP1 = YP1 < 0. ? 0.: YP1 > lim_y ? lim_y : YP1;
    y[idx] = YP1;
}


// ===== Simulation =====

void updateStatus()
{
    error_check(cudaMemset(d_s, 0, nx*ny*sizeof(int)));

    const dim3 block( num/block_size+1, 1, 1 );
    const dim3 thread( block_size, 1, 1 );

    _updatestatus<block_size><<<block, thread>>>(num, d_mx, d_my, d_s);
}

void advect()
{
    // copy velocity to backup buffers
    error_check(cudaMemcpy(d_bu, d_u, (nx+1)*ny*sizeof(double), cudaMemcpyDeviceToDevice));
    error_check(cudaMemcpy(d_bv, d_v, (ny+1)*nx*sizeof(double), cudaMemcpyDeviceToDevice));

    const dim3 block( nx/block_size + 1, ny/block_size + 1, 1 );
    const dim3 thread( block_size, block_size, 1 );

    // horizontal direction
    _advect<block_size><<<block, thread>>>(d_u, d_bu, d_bu, d_bv, dt, 0, 0.5, nx+1, ny);

    // vertical dircetion
    _advect<block_size><<<block, thread>>>(d_v, d_bv, d_bu, d_bv, dt, 0.5, 0, nx, ny+1);
}

void addForce()
{
    const dim3 block( nx/block_size + 1, ny/block_size + 1, 1 );
    const dim3 thread( block_size, block_size, 1 );

    // vertical direction
    _addforce<block_size><<<block, thread>>>(d_v, dt);
}


void advectMarkers()
{
    const dim3 block( num/block_size+1 , 1, 1 );
    const dim3 thread( block_size, 1, 1 );

    _advectmarkers<block_size><<<block, thread>>>(num, d_mx, d_my, d_u, d_v, dt, nx*dx, ny*dy);
}

// ===== functions =====

void initialize_grid()
{
    // == MEMORY ALLOCATE ==

    // create markers
    error_check(cudaMalloc(&d_mx, num*sizeof(double)));
    error_check(cudaMalloc(&d_my, num*sizeof(double)));

    // create grid cells
    error_check(cudaMalloc(&d_u, (nx+1)*ny*sizeof(double)));
    error_check(cudaMalloc(&d_v, (ny+1)*nx*sizeof(double)));
    error_check(cudaMalloc(&d_p, nx*ny*sizeof(double)));
    error_check(cudaMalloc(&d_s, nx*ny*sizeof(int)));

    // create backup buffers
    error_check(cudaMalloc(&d_bu, (nx+1)*ny*sizeof(double)));
    error_check(cudaMalloc(&d_bv, (ny+1)*nx*sizeof(double)));

    // == INITIALIZE ==

    // initialize  markers
    error_check(cudaMemcpy(d_mx, markers_x, num*sizeof(double), cudaMemcpyHostToDevice));
    error_check(cudaMemcpy(d_my, markers_y, num*sizeof(double), cudaMemcpyHostToDevice));

    // initialize grid cells & backup buffers
    error_check(cudaMemset(d_u, 0, (nx+1)*ny*sizeof(double)));
    error_check(cudaMemset(d_v, 0, (ny+1)*nx*sizeof(double)));
    error_check(cudaMemset(d_p, 0, nx*ny*sizeof(double)));
    error_check(cudaMemset(d_s, 0, nx*ny*sizeof(int)));

    error_check(cudaMemset(d_bu, 0, (nx+1)*ny*sizeof(double)));
    error_check(cudaMemset(d_bv, 0, (ny+1)*nx*sizeof(double)));
}

void initialize_params()
{
    error_check(cudaMemcpyToSymbol(d_nx, &nx, sizeof(int), 0, cudaMemcpyHostToDevice));
    error_check(cudaMemcpyToSymbol(d_ny, &ny, sizeof(int), 0, cudaMemcpyHostToDevice));    
    error_check(cudaMemcpyToSymbol(d_dx, &dx, sizeof(double), 0, cudaMemcpyHostToDevice));    
    error_check(cudaMemcpyToSymbol(d_dy, &dy, sizeof(double), 0, cudaMemcpyHostToDevice)); 

    error_check(cudaMemcpyToSymbol(d_nu, &nu, sizeof(double), 0, cudaMemcpyHostToDevice));    
    error_check(cudaMemcpyToSymbol(d_rho, &rho, sizeof(double), 0, cudaMemcpyHostToDevice));    
    error_check(cudaMemcpyToSymbol(d_g, &g, sizeof(double), 0, cudaMemcpyHostToDevice));
}

void get_result()
{
    error_check(cudaMemcpy(markers_x, d_mx, num*sizeof(double), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(markers_y, d_my, num*sizeof(double), cudaMemcpyDeviceToHost));
}

void finalize_grid()
{
    delete[] markers_x;
    delete[] markers_y;

    cudaFree(d_mx);
    cudaFree(d_my);

    cudaFree(d_u);
    cudaFree(d_v);
    cudaFree(d_p);
    cudaFree(d_s);

    cudaFree(d_bu);
    cudaFree(d_bv);
}


// ===== file IO =====
void write(const char *filename)
{
    std::ofstream fout(filename, std::ios::binary);

    fout.write((char*)&nx, sizeof(int));
    fout.write((char*)&ny, sizeof(int));
    fout.write((char*)&dx, sizeof(double));
    fout.write((char*)&dy, sizeof(double));

    fout.write((char*)&num, sizeof(int));

    for(int i=0;i<num;++i)
    {
        fout.write((char*)&markers_x[i], sizeof(double));
        fout.write((char*)&markers_y[i], sizeof(double));
    }
}

void read(const char *filename)
{
    std::ifstream fin(filename, std::ios::binary);
    
    fin.read((char*)&nx, sizeof(int));
    fin.read((char*)&ny, sizeof(int));
    fin.read((char*)&dx, sizeof(double));
    fin.read((char*)&dy, sizeof(double));

    fin.read((char*)&num, sizeof(int));

    markers_x = new double[num]{};
    markers_y = new double[num]{};

    for(int i=0;i<num;++i)
    {
        fin.read((char*)&markers_x[i], sizeof(double));
        fin.read((char*)&markers_y[i], sizeof(double));
    }
}

// ./exe input output steps dt nu rho max_iter tol
int main(int argc, char **argv)
{
    assert(argc == 9);
    
    // read map
    read(argv[1]);

    // initialize parameters
    steps = atoi(argv[3]);
    dt = atof(argv[4]);
    nu = atof(argv[5]);
    rho = atof(argv[6]);
    max_iter = atoi(argv[7]);
    tol = atof(argv[8]);


    // TODO: create grid
    initialize_grid();
    initialize_params();

    for(int i=0;i<steps;++i)
    {
        // TODO: update status
            updateStatus();
        // TODO: advect
            advect();
        // TODO: addForce
            addForce();
        // TODO: diffuse
        // TODO: project

        // TODO: extrapolate
        // TODO: constrain_velocity
        // TODO: move markers
            advectMarkers();

            get_result();

            char filename[30];
            sprintf(filename, "%s_%03i.sr", argv[2], i);
            write(filename);

    }

    // write map
    //write(argv[2]);

    // TODO: free grid
    finalize_grid();

    return 0;
}

