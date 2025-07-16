#ifndef _COMMON_H_
#define _COMMON_H_

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <pthread.h>
#include <assert.h>
#include <limits.h>
#include <fftw3.h>
#include <omp.h>

#define fftw_real double

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double Cut;
    double VLow;
    double kappa;
    double sdie;
    double kBT;
    double escl;
    double vscl;
} SOLV;

#define INFOLEN 31

typedef struct {
    char info[INFOLEN];
        double xyz[3];
        double q;
        double r;
    double Asq; //sqrt(A)
    double Bsq; //sqrt(B)
}ATOM;

// Type definitions
typedef signed char tErn;
#define tErn_MAX SCHAR_MAX
#define tErn_MIN SCHAR_MIN
#define tErn_SCL 8.0

// FFT functions
void fft3d_r2c(int L, int M, int N, fftw_real a[L][M][2*(N/2+1)]);
void fft3d_c2r(int L, int M, int N, fftw_real c[L][M][2*(N/2+1)]);
void fft3d_add_inplace(int L, int M, int N, fftw_real a[L][M][2*(N/2+1)], fftw_real b[L][M][2*(N/2+1)]);

// Grid functions
void volGrdR(int l, fftw_real *grd, int nAtom, ATOM *atoms, double sys[4]);
void debyeGrdR(int l, fftw_real *grd, int nAtom, ATOM *atoms, double sys[4]);
void debyeGrdL2nd(int l, fftw_real *grd, int nAtom, ATOM *atoms, double sys[4]);
void r6GrdR(int l, fftw_real *grd, int nAtom, ATOM *atoms, double sys[4]);
void r6GrdL(int l, fftw_real *grd, int nAtom, ATOM *atoms, double sys[4]);
void r12GrdR(int l, fftw_real *grd, int nAtom, ATOM *atoms, double sys[4]);
void r12GrdL(int l, fftw_real *grd, int nAtom, ATOM *atoms, double sys[4]);
int init_lattice_data();

// Score functions
void SetBool(int l, bool vol[l][l][l], bool label);
void SetZero(int l, fftw_real grd[l][l][2*(l/2+1)]);
void vRep(int l, bool vol[l][l][l], double cutoff, fftw_real grd[l][l][2*(l/2+1)], double scl, const char* str, double kBT);
void AddTo1filt(int l, fftw_real grd1[l][l][2*(l/2+1)], const double scl1, fftw_real grd2[l][l][2*(l/2+1)], const double scl2, bool vol[l][l][l]);
void toBool(int l, fftw_real *grid, bool *gridbool, fftw_real cutoff);
void showscorenet(const char* func, double sum, double frac, double u, double u1);

// Grid generation
void softgrd(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], int l, double dx, SOLV* sol, tErn softfb[l][l][l]);

// Cleanup functions
void fft3d_cleanup(void);
void softgrd_cleanup(void);

#ifdef __cplusplus
}
#endif

#endif  /* _COMMON_H_ */
