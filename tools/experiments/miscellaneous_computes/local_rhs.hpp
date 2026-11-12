#pragma once
// Version 1. Row-major double arrays are borrowed until mc_local_drop_v1.
// centers and means: M x n, projections and tangent frames: M x n x 2,
// coefficients: M x 6 x n in [1,u,v,u²,uv,v²] order, scales: n.
// k+1 <= M; the exact nearest k fitted charts are blended. The native
// specialization covers 2..64 coordinates and 1..64 neighbors.
extern "C" {
void* mc_local_field_v1(int n,int rows,int k,double damping,const double* centers,
    const double* means,const double* projections,const double* frames,
    const double* coefficients,const double* scales);
void mc_local_batch_v1(void* field,int count,const double* states,double* rates);
void mc_local_drop_v1(void* field);
}
