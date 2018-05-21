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
double rad = 0.08;
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
__device__ __constant__ double d_rad;

int cur_step = 0;

// ===== device =====
double *markers_x;
double *markers_y;

double *d_mx;  // x-coordinate of markers
double *d_my;  // y-coordinate of markers

double *d_u;  // horizontal component of velocity
double *d_v;  // vertical component of velocity
double *d_p;  // pressure
int *d_st;  // status
int *d_pt;
double *d_phi;

int *d_uv;  // valid u
int *d_vv;  // valid v
int *d_buv;  //backup valid u
int *d_bvv;  //backup valid v

double *d_bu;  // backup u
double *d_bv;  // backup v

double *d_A;
int *d_cooRowIdx;
int *d_rowIdx;
int *d_colIdx;
double *d_b;

int *d_neig;  // neighbor count of A for each cell
int *d_idxmap;

int *d_rbuf=0;

// ===== cpu memory =====
int *rbuf=0;
int *st;
int *neig;
int *idxmap;

double *u;
double *v;
double *pressure;

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

template<typename T>
__device__ inline T _safe_get_mr(const T* const field, int x, int y,
                                const int lim_x, const int lim_y)
{
    if(x < 0)
        x = 0;
    else if(x >= lim_x)
        x = lim_x-1;
    if(y < 0)
        y = 0;
    else if(y >= lim_y)
        y = lim_y-1;

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

template<typename T>
__device__ inline void _safe_add(T* const field, const int x, const int y,
                                const int lim_x, const int lim_y, const T value)
{
    if(x < 0 || x >= lim_x || y < 0 || y >= lim_y)
        return;
    atomicAdd(&field[y*lim_x + x], value);
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
    double fy1 = (1.-fcx) * _safe_get_mr(field, x0, y0, lim_x, lim_y) 
                    + fcx * _safe_get_mr(field, x1, y0, lim_x, lim_y);
    double fy2 = (1.-fcx) * _safe_get_mr(field, x0, y1, lim_x, lim_y)
                    + fcx * _safe_get_mr(field, x1, y1, lim_x, lim_y);
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

__device__ inline void _RK2(const double* const u, const double* const v, const double dt,
                            double XG, double YG, double *XP, double *YP,
                            const int lim_ux, const int lim_uy, const int lim_vx, const int lim_vy)
{
    const double hx = XG - d_dx/2.;
    const double hy = YG - d_dy/2.;

    const double xk1 = dt * _interpolate(u, XG, hy, lim_ux, lim_uy);
    const double yk1 = dt * _interpolate(v, hx, YG, lim_vx, lim_vy);

    const double xk2 = dt * _interpolate(u, XG + xk1/2., hy + yk1/2., lim_ux, lim_uy);
    const double yk2 = dt * _interpolate(v, hx + xk1/2., YG + yk1/2., lim_vx, lim_vy);

    *XP = XG + xk2;
    *YP = YG + yk2;
}

template<int block_size>
__global__ void _initialize(double* const grid, const int lim_x, const int lim_y, const double value)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * lim_x + x;

    if(y >= lim_y || x >= lim_x)return;

    grid[idx] = value;
}

// Update Grid Status (kernel)
template<int block_size>
__global__ void _updatestatus(const int num, const double* const x, const double* const y, int* const st, int* const pt)
{
    // thread index
    const int idx = block_size * blockIdx.x + threadIdx.x;

    // out of bounds
    if(idx >= num) return;

    // index coord
    const int X = (int)floor(x[idx]/d_dx);
    const int Y = (int)floor(y[idx]/d_dy);

    const int v = 1;
    _safe_set(st, X, Y, d_nx, d_ny, v);
    _safe_add(pt, X, Y, d_nx, d_ny, v);
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
    _RK2(field_u, field_v, -dt, XG, YG, &XP, &YP, d_nx+1, d_ny, d_nx, d_ny+1);

    // update 
    u[idx] = _interpolate(bu, XP-off_x*d_dx, YP-off_y*d_dy, lim_x, lim_y);
}


// Add Force (kernel)
template<int block_size>
__global__ void _addforce(int* const status, double* const v, const double dt)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y*d_nx + x;

    // out of boundaries
    if(y >= d_ny+1 || x >= d_nx)return;

    if(status[idx] == 0)return;

    // add gravity
    v[idx] = v[idx] - dt * d_g;
}

template<int block_size>
__global__ void _compute_SDF(double* const mx, double* const my, const int num, double* const phi)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    if(y >= d_ny || x >= d_nx)return;

    double ddx = ((double)x + 0.5)/d_dx - mx[num];
    double ddy = ((double)y + 0.5)/d_dy - my[num];
    
    double dist = sqrt(ddx*ddx + ddy*ddy) - d_rad;

    if(dist < phi[idx])
        phi[idx] = dist;
}

template<int block_size>
__global__ void _pre_counting(const int* const status, int* const array)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    // out of bounds
    if(y >= d_ny || x >= d_nx)return;

    // not fluid
    if(status[idx] == 0)return;

    int nb=0;

    if(x > 0 && status[idx-1] != 0)nb++;    
    if(y > 0 && status[idx-d_nx] != 0)nb++;

    array[idx] = nb;
}

template<int block_size>
__global__ void _int_reduce(const int* const input, int* const output, unsigned int n)
{
    extern __shared__ int sdata[];

    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * block_size * 2 + threadIdx.x;
    unsigned int grid_size = block_size * 2 * gridDim.x;

    int sum = 0;
    while(i < n)
    {
        sum += input[i];
        if(i + block_size < n) sum += input[i+block_size];
        i += grid_size;
    }
    sdata[tid] = sum;
    __syncthreads();
    if((block_size >= 512) && (tid < 256)) sdata[tid] = sum = sum + sdata[tid+256];
    __syncthreads();
    if((block_size >= 256) && (tid < 128)) sdata[tid] = sum = sum + sdata[tid+128];
    __syncthreads();
    if((block_size >= 128) && (tid < 64)) sdata[tid] = sum = sum+sdata[tid+64];
    __syncthreads();
    if((block_size >= 64) && (tid < 32)) sdata[tid] = sum = sum+sdata[tid+32];
    __syncthreads();
    if((block_size >= 32) && (tid < 16)) sdata[tid] = sum = sum+sdata[tid+16];
    __syncthreads();
    if((block_size >= 16) && (tid < 8)) sdata[tid] = sum = sum+sdata[tid+8];
    __syncthreads();
    if((block_size >= 8) && (tid < 4)) sdata[tid] = sum = sum+sdata[tid+4];
    __syncthreads();
    if((block_size >= 4) && (tid < 2)) sdata[tid] = sum = sum+sdata[tid+2];
    __syncthreads();
    if((block_size >= 2) && (tid < 1)) sdata[tid] = sum = sum+sdata[tid+1];
    __syncthreads();
    if(tid == 0) output[blockIdx.x] = sum;
}

template<int block_size>
__global__ void _construct_sparse_matrix(int* const status, int* const idxmap, int* const terms,
                                        double* const A, int* const cooRowIdx, int* const colIdx,
                                        const double dt)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    // out of bounds
    if(y >= d_ny || x >= d_nx)return;

    // not fluid
    if(status[idx] == 0)return;

    double scale = dt/d_rho;
    
    int pos = terms[idx];
    int offset = -1;
    double dg = 0;

    double scale_x = scale/(d_dx*d_dx);
    double scale_y = scale/(d_dy*d_dy);
    
    if(x > 0)
    {
        dg += scale_x;
        if(status[idx-1] != 0)
        {
            A[pos+offset] = -scale_x;
            colIdx[pos+offset] = idxmap[idx-1];
            cooRowIdx[pos+offset] = idxmap[idx];
            offset--;
        }
    }
    if(y > 0)
    {
        dg += scale_y;
        if(status[idx-d_nx] != 0)
        {
            A[pos+offset] = -scale_y;
            colIdx[pos+offset] = idxmap[idx-d_nx];
            cooRowIdx[pos+offset] = idxmap[idx];
        }
    }
    if(x < d_nx-1)dg += scale_x;
    if(y < d_ny-1)dg += scale_y;

    A[pos] = dg;
    colIdx[pos] = idxmap[idx];
    cooRowIdx[pos] = idxmap[idx];
}

template<int block_size>
__global__ void _construct_right_hand_side(int* const status, int* const idxmap, 
                                    double* const u, double* const v, double* const b)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    // out of bounds
    if(y >= d_ny || x >= d_nx)return;

    // not fluid
    if(status[idx] == 0)return;

    const int uidx = idx + y;
    const int pos = idxmap[idx];

    double d = 0;
    d += (u[uidx+1]-u[uidx])/d_dx;
    d += (v[idx+d_nx]-v[idx])/d_dy;

    b[pos] = -d;
}

template<int block_size>
__global__ void _update_pressure(int* const status, int* const idxmap,
                            double* const p, double* const xp)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    // out of bounds
    if(y >= d_ny || x >= d_nx)return;

    const int pos = idxmap[idx];
    
    double v;
    // not fluid
    if(status[idx] == 0)
        v = 0;
    else
        v = xp[pos];

    p[idx] = v;
}

template<int block_size>
__global__ void _update_velocity_v_by_pressure(int* const status, int* const valid, int* const pt, double* const p,
                                        double *v, const double dt)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    if(y >= d_ny || x >= d_nx || y < 1)return;
    if(status[idx] == 0 && status[idx-d_nx] == 0)return;

    double pidx = p[idx];
    double pidx_1 = p[idx-d_nx];

    // negative gradient
    const double gd = -dt/d_rho * (pidx-pidx_1)/d_dy;

    // update v
    v[idx] = v[idx] + gd;
    valid[idx] = 1;
}

template<int block_size>
__global__ void _update_velocity_u_by_pressure(int* const status, int* const valid, int* const pt, double* const p,
                                        double *u, const double dt)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    if(y >= d_ny || x >= d_nx || x < 1)return;
    if(status[idx] == 0 && status[idx-1] == 0)return;

    const int uidx = idx + y;

    double pidx = p[idx];// + ((double)pt[idx]/9.);
    double pidx_1 = p[idx-1];// + ((double)pt[idx-1]/9.);

    // negative gradient
    const double gd = -dt/d_rho * (pidx-pidx_1)/d_dx;
    
    // update u
    u[uidx] = u[uidx] + gd;
    valid[uidx] = 1;
}

template<int block_size>
__global__ void _enforce_boundary_v(double* const v)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    if(y != d_ny && y != 0 || x >= d_nx)return;

    if(y == 0 && v[idx] < 0)v[idx] = 0;
    else if(y == d_ny && v[idx] > 0)v[idx] = 0;
}

template<int block_size>
__global__ void _enforce_boundary_u(double* const u)
{
     // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * (d_nx+1) + x;

    if(x != d_nx && x != 0 || y >= d_ny)return;

    if(x == 0 && u[idx] < 0)u[idx] = 0;
    else if(x == d_nx && u[idx] > 0)u[idx] = 0;
}

template<int block_size>
__global__ void _clean_field(int* const status, double* const u, double* const v)
{
     // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    if(x >= d_nx || y >= d_ny)return;
    
    if(_safe_get(status, x, y, d_nx, d_ny) == 0)
    {
        if(_safe_get(status, x-1, y, d_nx, d_ny) == 0)
            u[idx+y] = 0;
        if(_safe_get(status, x, y-1, d_nx, d_ny) == 0)
            v[idx] = 0;

        if(x != d_nx-1 || y != d_ny-1)return;

        if(_safe_get(status, x+1, y, d_nx, d_ny) == 0)
            u[idx+y+1] = 0;
        if(_safe_get(status, x, y+1, d_nx, d_ny) == 0)
            v[idx+y] = 0;
    }
}

template<int block_size>
__global__ void _extrapolate_u(int* const nv, int* const v, double* const nu, double* const u)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * (d_nx+1) + x;

    if(y >= d_ny-1 || y < 1 || x >= d_nx || x < 1)return;
    if(v[idx] != 0)return;

    double sum = 0;
    int count = 0;
    
    if(v[idx+1] != 0)
    {
        sum += u[idx+1];
        count++;
    }

    if(v[idx-1] != 0)
    {
        sum += u[idx-1];
        count++;
    }

    if(v[idx+d_nx] != 0)
    {
        sum += u[idx+d_nx];
        count++;
    }

    if(v[idx-d_nx] != 0)
    {
        sum += u[idx-d_nx];
        count++;
    }

    if(count > 0)
    {
        nu[idx] = sum/(double)count;
        nv[idx] = 1;
    }
}

template<int block_size>
__global__ void _extrapolate_v(int* const nv, int* const v, double* const nu, double* const u)
{
    // thread index
    const int y = block_size * blockIdx.y + threadIdx.y;
    const int x = block_size * blockIdx.x + threadIdx.x;
    const int idx = y * d_nx + x;

    if(y >= d_ny || y < 1 || x >= d_nx-1 || x < 1)return;
    if(v[idx] != 0)return;

    double sum = 0;
    int count = 0;
    
    if(v[idx+1] != 0)
    {
        sum += u[idx+1];
        count++;
    }

    if(v[idx-1] != 0)
    {
        sum += u[idx-1];
        count++;
    }

    if(v[idx+d_nx] != 0)
    {
        sum += u[idx+d_nx];
        count++;
    }

    if(v[idx-d_nx] != 0)
    {
        sum += u[idx-d_nx];
        count++;
    }

    if(count > 0)
    {
        nu[idx] = sum/(double)count;
        nv[idx] = 1;
    }
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
    _RK2(field_u, field_v, dt, XP, YP, &XP1, &YP1, d_nx+1, d_ny, d_nx, d_ny+1);

    // clamp into boundaries & set new position
    XP1 = XP1 < 0. ? 0.: XP1 > lim_x ? lim_x : XP1;
    x[idx] = XP1;
    
    YP1 = YP1 < 0. ? 0.: YP1 > lim_y ? lim_y : YP1;
    y[idx] = YP1;
}

// ===== Simulation =====

void updateStatus()
{
    error_check(cudaMemset(d_st, 0, nx*ny*sizeof(int)));
    error_check(cudaMemset(d_pt, 0, nx*ny*sizeof(int)));

    const dim3 block( num/block_size+1, 1, 1 );
    const dim3 thread( block_size, 1, 1 );

    _updatestatus<block_size><<<block, thread>>>(num, d_mx, d_my, d_st, d_pt);
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
    _addforce<block_size><<<block, thread>>>(d_st, d_v, dt);
}

template<int threads>
int int_reduce(int *d_array, unsigned int size)
{
    const dim3 thread(threads, 1, 1);
    const dim3 block( (size+(threads*2-1))/(threads*2), 1, 1);

    if(rbuf == 0)
        rbuf = new int[block.x]{};
    if(d_rbuf == 0)
        error_check(cudaMalloc(&d_rbuf, block.x*sizeof(int)));

    int sz = (threads <= 32) ? 2*threads*sizeof(int) : threads * sizeof(int);
    _int_reduce<threads><<<block, thread, sz>>>(d_array, d_rbuf, size);

    error_check(cudaMemcpy(rbuf, d_rbuf, block.x*sizeof(int), cudaMemcpyDeviceToHost));

    int sum = 0;
    for(int i=0;i<block.x;++i)
        sum += rbuf[i];
    
    return sum;
}

void compute_SDF()
{
    error_check(cudaMemset(d_phi, 0, nx*ny*sizeof(double)));
    const int sz = 32;
    const dim3 block(nx/sz+1, ny/sz+1, 1);
    const dim3 thread(sz, sz, 1);

    _initialize<sz><<<block, thread>>>(d_phi, nx, ny, rad*3);

    for(int i=0;i<num;++i)
    {
        _compute_SDF<sz><<<block, thread>>>(d_mx, d_my, i, d_phi);
    }
}

void project(PCGsolver &solver)
{
    // counting neighbors
    error_check(cudaMemset(d_neig, 0, nx*ny*sizeof(int)));

    const int sz = 32;
    const dim3 block(nx/sz+1, ny/sz+1, 1);
    const dim3 thread(sz, sz, 1);
    _pre_counting<sz><<<block, thread>>>(d_st, d_neig); //count neighbor

    error_check(cudaMemcpy(st, d_st, nx*ny*sizeof(int), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(neig, d_neig, nx*ny*sizeof(int), cudaMemcpyDeviceToHost));

    int dg = int_reduce<512>(d_st, nx*ny);
    int ng = int_reduce<512>(d_neig, nx*ny);

    // construct sparse matrix A & right-hand-side vector b
    int N = dg;
    int nonzero = dg + ng;

    memset(idxmap, 0, nx*ny*sizeof(int));

    // mapping
    int n = 0;
    int m = 0;
    for(int i=0;i<nx*ny;++i)
    {
        if(st[i] == 1)
        {
            idxmap[i] = n++;
            m += neig[i];
            neig[i] = m;
            m++;
        }
    }

    error_check(cudaMemcpy(d_idxmap, idxmap, nx*ny*sizeof(int), cudaMemcpyHostToDevice));
    error_check(cudaMemcpy(d_neig, neig, nx*ny*sizeof(int), cudaMemcpyHostToDevice));

    error_check(cudaMalloc(&d_A, nonzero*sizeof(double)));
    error_check(cudaMalloc(&d_cooRowIdx, nonzero*sizeof(int)));
    error_check(cudaMalloc(&d_rowIdx, (N+1)*sizeof(int)));
    error_check(cudaMalloc(&d_colIdx, nonzero*sizeof(int)));
    error_check(cudaMalloc(&d_b, N*sizeof(double)));

    // construct sparse matrix A (lower triangular)
    _construct_sparse_matrix<sz><<<block, thread>>>(d_st, d_idxmap, d_neig, d_A, d_cooRowIdx, d_colIdx, dt);
    // construct right hand side vector b
    _construct_right_hand_side<sz><<<block, thread>>>(d_st, d_idxmap, d_u, d_v, d_b);

    // convert from coo to csr format
    solver.convert_coo2csr(N, nonzero, d_cooRowIdx, d_rowIdx);

    // solve linear system
    solver.solve_gpumem(N, nonzero, d_A, d_rowIdx, d_colIdx, d_b, NULL);

    // get answer
    double *d_x = solver.get_device_x();

    // debug.....
    double *A = new double[nonzero]{};
    int *cooRowIdx = new int[nonzero]{};
    int *colIdx = new int[nonzero]{};
    double *b = new double[N]{};
    double *xp = new double[N]{};

    //std::cout << "N: " << N << std::endl;
    //std::cout << "nz: " << nonzero << std::endl;

    error_check(cudaMemcpy(A, d_A, nonzero*sizeof(double), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(cooRowIdx, d_cooRowIdx, (nonzero)*sizeof(int), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(colIdx, d_colIdx, (nonzero)*sizeof(int), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(b, d_b, N*sizeof(double), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(xp, d_x, N*sizeof(double), cudaMemcpyDeviceToHost));


    char name[30];
    sprintf(name, "mat/%s_%03d.in", "step", cur_step);

    std::ofstream fout(name, std::ios::binary);

    fout.write((char*)&N, sizeof(int));
    fout.write((char*)&nonzero, sizeof(int));

    for(int i=0;i<nonzero;++i)
    {
        fout.write((char*)&colIdx[i], sizeof(int));
        fout.write((char*)&cooRowIdx[i], sizeof(int));
        fout.write((char*)&A[i], sizeof(double));
    }

    for(int i=0;i<N;++i)
    {
        fout.write((char*)&b[i], sizeof(double));
    }
    fout.close();
/*
    std::cout << "A: ";
    for(int i=0;i<nonzero;++i)
        std::cout << A[i] << ", ";
    std::cout << std::endl;
    
    std::cout << "row: ";
    for(int i=0;i<nonzero;++i)
        std::cout << cooRowIdx[i] << ", ";
    std::cout << std::endl;

    std::cout << "col: ";
    for(int i=0;i<nonzero;++i)
        std::cout << colIdx[i] << ", ";
    std::cout << std::endl;
 
    std::cout << "b: ";
    for(int i=0;i<N;++i)
        std::cout << b[i] << ", ";
    std::cout << std::endl;
    
    std::cout << "xp: ";
    for(int i=0;i<N;++i)
        std::cout << xp[i] << ", ";
    std::cout << std::endl;
*/

    delete[] A;
    delete[] cooRowIdx;
    delete[] colIdx;
    delete[] b;
    delete[] xp;

    //..... debug



    cudaFree(d_A);
    cudaFree(d_cooRowIdx);
    cudaFree(d_rowIdx);
    cudaFree(d_colIdx);
    cudaFree(d_b);

    error_check(cudaMemset(d_p, 0, nx*ny*sizeof(double)));
    error_check(cudaMemset(d_uv, 0, (nx+1)*ny*sizeof(int)));
    error_check(cudaMemset(d_vv, 0, (ny+1)*nx*sizeof(int)));

    _update_pressure<sz><<<block, thread>>>(d_st, d_idxmap, d_p, d_x);
    _update_velocity_u_by_pressure<sz><<<block, thread>>>(d_st, d_uv, d_pt, d_p, d_u, dt);
    _update_velocity_v_by_pressure<sz><<<block, thread>>>(d_st, d_vv, d_pt, d_p, d_v, dt);
}


void extrapolate_u()
{
    const dim3 block(nx/block_size+1, ny/block_size+1, 1);
    const dim3 thread(block_size, block_size, 1);

    error_check(cudaMemcpy(d_buv, d_uv, (nx+1)*ny*sizeof(int), cudaMemcpyDeviceToDevice));
    error_check(cudaMemcpy(d_bu, d_u, (nx+1)*ny*sizeof(double), cudaMemcpyDeviceToDevice));

    _extrapolate_u<block_size><<<block, thread>>>(d_uv, d_buv, d_u, d_bu);

}

void extrapolate_v()
{
    const dim3 block(nx/block_size+1, ny/block_size+1, 1);
    const dim3 thread(block_size, block_size, 1);

    error_check(cudaMemcpy(d_bvv, d_vv, (ny+1)*nx*sizeof(int), cudaMemcpyDeviceToDevice));
    error_check(cudaMemcpy(d_bv, d_v, (ny+1)*nx*sizeof(double), cudaMemcpyDeviceToDevice));

    _extrapolate_v<block_size><<<block, thread>>>(d_vv, d_bvv, d_v, d_bv);
}

void extrapolate()
{
    for(int i=0;i<4;++i)
    {
        extrapolate_u();
    }

    for(int i=0;i<4;++i)
    {
        extrapolate_v();
    }

}

void enforce_boundary()
{
    const dim3 block( nx/block_size+1, ny/block_size+1, 1);
    const dim3 thread( block_size, block_size, 1);

    _enforce_boundary_u<block_size><<<block, thread>>>(d_u);
    _enforce_boundary_v<block_size><<<block, thread>>>(d_v);
}

void clean_field()
{
    const dim3 block(nx/block_size+1, ny/block_size+1, 1);
    const dim3 thread(block_size, block_size, 1);

    _clean_field<block_size><<<block, thread>>>(d_st, d_u, d_v);
}

void advectMarkers()
{
    const dim3 block( num/block_size+1 , 1, 1 );
    const dim3 thread( block_size, 1, 1 );

    for(int i=0;i<5;++i)
        _advectmarkers<block_size><<<block, thread>>>(num, d_mx, d_my, d_u, d_v, 0.2*dt, nx*dx, ny*dy);
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
    error_check(cudaMalloc(&d_st, nx*ny*sizeof(int)));
    error_check(cudaMalloc(&d_pt, nx*ny*sizeof(int)));
    error_check(cudaMalloc(&d_phi, nx*ny*sizeof(double)));

    error_check(cudaMalloc(&d_uv, (nx+1)*ny*sizeof(int)));
    error_check(cudaMalloc(&d_vv, (ny+1)*nx*sizeof(int)));
    error_check(cudaMalloc(&d_buv, (nx+1)*ny*sizeof(int)));
    error_check(cudaMalloc(&d_bvv, (ny+1)*nx*sizeof(int)));

    // create backup buffers
    error_check(cudaMalloc(&d_bu, (nx+1)*ny*sizeof(double)));
    error_check(cudaMalloc(&d_bv, (ny+1)*nx*sizeof(double)));

    // create other buffers
    error_check(cudaMalloc(&d_neig, nx*ny*sizeof(int)));
    error_check(cudaMalloc(&d_idxmap, nx*ny*sizeof(int)));

    st = new int[nx*ny]{};
    neig = new int[nx*ny]{};
    idxmap = new int[nx*ny]{};

    u = new double[(nx+1)*ny]{};
    v = new double[nx*(ny+1)]{};
    pressure = new double[nx*ny]{};

    // == INITIALIZE ==

    // initialize  markers
    error_check(cudaMemcpy(d_mx, markers_x, num*sizeof(double), cudaMemcpyHostToDevice));
    error_check(cudaMemcpy(d_my, markers_y, num*sizeof(double), cudaMemcpyHostToDevice));

    // initialize grid cells & backup buffers
    error_check(cudaMemset(d_u, 0, (nx+1)*ny*sizeof(double)));
    error_check(cudaMemset(d_v, 0, (ny+1)*nx*sizeof(double)));
    error_check(cudaMemset(d_p, 0, nx*ny*sizeof(double)));
    error_check(cudaMemset(d_st, 0, nx*ny*sizeof(int)));

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
    error_check(cudaMemcpyToSymbol(d_rad, &rad, sizeof(double), 0, cudaMemcpyHostToDevice));
}

void get_result()
{
    error_check(cudaMemcpy(markers_x, d_mx, num*sizeof(double), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(markers_y, d_my, num*sizeof(double), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(u, d_u, (nx+1)*ny*sizeof(double), cudaMemcpyDeviceToHost));
    error_check(cudaMemcpy(v, d_v, ny*(nx+1)*sizeof(double), cudaMemcpyDeviceToHost));

    error_check(cudaMemcpy(pressure, d_p, nx*ny*sizeof(double), cudaMemcpyDeviceToHost));

}

void finalize_grid()
{
    delete[] markers_x;
    delete[] markers_y;
    delete[] rbuf;
    delete[] st;
    delete[] neig;
    delete[] idxmap;

    cudaFree(d_mx);
    cudaFree(d_my);

    cudaFree(d_u);
    cudaFree(d_v);
    cudaFree(d_p);
    cudaFree(d_st);
    cudaFree(d_phi);

    cudaFree(d_pt);

    cudaFree(d_uv);
    cudaFree(d_vv);
    cudaFree(d_buv);
    cudaFree(d_bvv);

    cudaFree(d_bu);
    cudaFree(d_bv);

    cudaFree(d_neig);

    cudaFree(d_rbuf);
    cudaFree(d_idxmap);

    delete[] u;
    delete[] v;
    delete[] pressure;
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

    for(int i=0;i<(nx+1)*ny;++i)
    {
        fout.write((char*)&u[i], sizeof(double));
    }

    for(int i=0;i<nx*(ny+1);++i)
    {
        fout.write((char*)&v[i], sizeof(double));
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

    // TODO: create PCGsolver
    PCGsolver p_solver(max_iter, tol);

    char filename[30];

    for(int i=0;i<steps;++i)
    {
        // TODO: update status
            updateStatus();
        // TODO: advect
            advect();
        // TODO: addForce
            addForce();
        // TODO: enforce boundary
            enforce_boundary();
        // TODO: project
            project(p_solver);
        // TODO: clean field
            clean_field();
        // TODO: extrapolate
            extrapolate();
        // TODO: advect markers
            advectMarkers();

        // Debugging...
            get_result();

            char filename[30];
            sprintf(filename, "%s_all_%03i.sr", argv[2], i);
            write(filename);

    }

    // write map
    //get_result();
    //write(argv[2]);

    // TODO: free grid
    finalize_grid();

    return 0;
}

