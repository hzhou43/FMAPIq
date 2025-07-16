#include "common.h"


/**
 * Safe exponential function that prevents overflow
 * @param x: Input value for exponential
 * @return: exp(x) or 0.0 if x would cause overflow
 */
double myEXP(double x) {
    if (x > FLT_MAX_EXP) {
        return 0.0;
    } else {
        return exp(x);
    }
}

/**
 * Display scoring network results in formatted output
 * @param func: Function name/identifier
 * @param sum: Sum value
 * @param frac: Fraction/occupancy value
 * @param u: Average value
 * @param u1: Net value
 */
void showscorenet(const char* func, double sum, double frac, double u, double u1) {
    printf("%s\t", func);
    printf("Sum: %16e Occ:%16.12f Ave:%16.8f Net:%16.8f\n", sum, frac, u, u1);
}

/**
 * Convert real grid values to boolean grid based on cutoff
 * @param l: Grid dimension
 * @param grid: Input real grid data
 * @param gridbool: Output boolean grid
 * @param cutoff: Threshold value for conversion
 */
void toBool(int l, fftw_real *grid, bool *gridbool, fftw_real cutoff) {
    fftw_real (*grd)[l][2*(l/2+1)] = (void*)grid;
    bool (*gbl)[l][l] = (void*)gridbool;
    
    int i, j, k;
    #pragma omp parallel for private(j,k)
    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < l; k++) {
                gbl[i][j][k] = (grd[i][j][k] > cutoff);
            }
        }
    }
}

/**
 * Set all values in boolean volume to specified label
 * @param l: Grid dimension
 * @param vol: 3D boolean array
 * @param label: Boolean value to set
 */
void SetBool(int l, bool vol[l][l][l], bool label) {
    int i, j, k;
    #pragma omp parallel for private(j,k)
    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < l; k++) {
                vol[i][j][k] = label;
            }
        }
    }
}

/**
 * Initialize all grid values to zero
 * @param l: Grid dimension
 * @param grd: 3D real grid array
 */
void SetZero(int l, fftw_real grd[l][l][2*(l/2+1)]) {
    int i, j, k;
    #pragma omp parallel for private(j,k)
    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < 2*(l/2+1); k++) {
                grd[i][j][k] = 0.0;
            }
        }
    }
}

/**
 * Volume reporting function with thermodynamic calculations
 * @param l: Grid dimension
 * @param vol: Boolean volume mask
 * @param cutoff: Cutoff value for volume determination
 * @param grd: Real grid data
 * @param scl: Scaling factor
 * @param str: String identifier for output
 * @param kBT: Boltzmann constant * temperature
 */
void vRep(int l, bool vol[l][l][l], double cutoff, fftw_real grd[l][l][2*(l/2+1)], 
          double scl, const char* str, double kBT) {
    double sum = 0.0;
    double count = 0.0;
    double dl3 = (double)l * (double)l * (double)l;
    
    // Initialize volume if string starts with "vol"
    if (!strncmp(str, "vol", 3)) {
        toBool(l, &grd[0][0][0], &vol[0][0][0], cutoff);
    }
    
    // Calculate sum and count - parallelize the most expensive part
    int i, j, k;
    #pragma omp parallel for private(j,k) reduction(+:sum,count)
    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < l; k++) {
                if (vol[i][j][k]) {
                    count += 1.0;
                } else {
                    sum += myEXP(grd[i][j][k] * scl / (-kBT));
                }
            }
        }
    }
    
    // Calculate thermodynamic quantities
    double u, un;
    if (!strncmp(str, "vol", 3)) {
        u = -kBT * log(1.0 - count/dl3);
        un = 0.0;
    } else {
        u = -kBT * log(sum / dl3);
        un = -kBT * log(sum / (dl3 - count));
    }
    
    showscorenet(str, sum, count/dl3, u, un);
}

/**
 * Add scaled values from one grid to another
 * @param l: Grid dimension
 * @param grd1: First grid (modified in place)
 * @param scl1: Scaling factor for first grid
 * @param grd2: Second grid (read only)
 * @param scl2: Scaling factor for second grid
 */
void AddTo1filt(int l, fftw_real grd1[l][l][2*(l/2+1)], const double scl1, fftw_real grd2[l][l][2*(l/2+1)], const double scl2,bool vol[l][l][l] ){
        int i,j,k;
        #pragma omp parallel for private(j,k)
        for (i=0;i<l;i++){
                for (j=0;j<l;j++){
                        for (k=0;k<l;k++){
                    if (vol[i][j][k]){
                        grd1[i][j][k]=grd1[i][j][k]*scl1+FLT_MAX_EXP;
                    }else{
                                            grd1[i][j][k]=grd1[i][j][k]*scl1+grd2[i][j][k]*scl2;
                    }
                        }
                }
        }
}
