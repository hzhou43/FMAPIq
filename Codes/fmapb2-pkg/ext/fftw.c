#include "common.h"
#include <omp.h>

void fft3d_r2c(int L, int M, int N, fftw_real a[L][M][2*(N/2+1)]){
	fftw_complex *A;
	fftw_plan p;
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

	A = (fftw_complex*) &a[0][0][0];
	p=fftw_plan_dft_r2c_3d(L, M, N,&a[0][0][0], A,
                                  FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
    fftw_cleanup_threads();
	//fftw_cleanup();
}

void fft3d_c2r(int L, int M, int N, fftw_real c[L][M][2*(N/2+1)]){
        fftw_complex *C;
	fftw_plan p;
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    C = (fftw_complex*) &c[0][0][0];
    p=fftw_plan_dft_c2r_3d(L, M, N, C, &c[0][0][0],
                                  FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_cleanup_threads();
	//fftw_cleanup();
}

/* cross-correlation, complex conjugate f.g*:
(a0+a1i)*(b0-b1i)=a0b0+a1b1+(a1*b0-a0*b1)i */
void fft3d_add_inplace(int L, int M, int N,fftw_real a[L][M][2*(N/2+1)], fftw_real b[L][M][2*(N/2+1)]){
	fftw_real scale = 1.0 / ((fftw_real)L*(fftw_real) M * (fftw_real)N);
	fftw_complex *A, *B;
	A = (fftw_complex*) &a[0][0][0];
	B = (fftw_complex*) &b[0][0][0];
	int i, j, k;
        fftw_complex tmpC;
	#pragma omp parallel for private (j,k,tmpC)	
	for (i = 0; i < L; ++i)
	    for (j = 0; j < M; ++j)
	        for (k = 0; k < N/2+1; ++k){
                int ijk = k+ (N/2+1)*(i*M + j);
                tmpC[0]=(A[ijk][0] * B[ijk][0] + A[ijk][1] * B[ijk][1])* scale;
                tmpC[1]=(A[ijk][1] * B[ijk][0] - A[ijk][0] * B[ijk][1])* scale;
		B[ijk][0]=tmpC[0];
                B[ijk][1]=tmpC[1];
                }
}
