#ifndef RDF_H
#define RDF_H

#include <math.h>

#ifdef __CUDACC__
// GPU device versions (CUDA only)
__device__ static inline double rdf_vol(double r2, double rc2) {
    return (r2 > rc2) ? 0.0 : 1.0;
}

__device__ static inline double rdf_debye(double r2, double kappa) {
    double r = sqrt(r2);
    return exp(-r * kappa) / r;
}

__device__ static inline double rdf_vdw(double r2, double n) {
    double r6 = r2 * r2 * r2;
    return (n < 4.0) ? 1.0/r6 : 1.0/(r6*r6);
}

#else
// CPU versions (regular C/C++)
static inline double rdf_vol(double r2, double rc2) {
    return (r2 > rc2) ? 0.0 : 1.0;
}

static inline double rdf_debye(double r2, double kappa) {
    double r = sqrt(r2);
    return exp(-r * kappa) / r;
}

static inline double rdf_vdw(double r2, double n) {
    double r6 = r2 * r2 * r2;
    return (n < 4.0) ? 1.0/r6 : 1.0/(r6*r6);
}

#endif // __CUDACC__

#endif // RDF_H
