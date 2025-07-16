#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "fft.h"

// Create FFT grid from real-space grid
fftwf_complex* calculate_fft(fftwf_complex* fft_grid, int grid_size) {
    // Create FFT plan
    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    
    fftwf_plan plan = fftwf_plan_dft_3d(grid_size, grid_size, grid_size,
                                        fft_grid, fft_grid, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!plan) {
        fprintf(stderr, "Error: Failed to create FFT plan\n");
        return NULL;
    }
    
    // Execute FFT
    fftwf_execute(plan);
    
    // Clean up plan
    fftwf_destroy_plan(plan);
    fftwf_cleanup_threads();
    
    return fft_grid;
}

// Create q-vectors for binning
QVectors* create_q_vectors(int grid_size, double grid_spacing) {
    QVectors* qvecs = malloc(sizeof(QVectors));
    if (!qvecs) {
        fprintf(stderr, "Error: Failed to allocate memory for q-vectors\n");
        return NULL;
    }
    
    // Set grid size
    qvecs->box[0] = qvecs->box[1] = qvecs->box[2] = grid_size;
    qvecs->spacing = grid_spacing;
    qvecs->total = (size_t)grid_size * (size_t)grid_size * (size_t)grid_size;
    
    // Allocate arrays for q-vectors
    qvecs->qx = malloc(grid_size * sizeof(double));
    qvecs->qy = malloc(grid_size * sizeof(double));
    qvecs->qz = malloc(grid_size * sizeof(double));
    
    if (!qvecs->qx || !qvecs->qy || !qvecs->qz) {
        fprintf(stderr, "Error: Failed to allocate memory for q-vector arrays\n");
        free_q_vectors(qvecs);
        return NULL;
    }
    
    // Calculate q-vectors
    const double two_pi = 2.0 * (double)M_PI;
    const double factor = two_pi / (grid_size * grid_spacing);
    const int half_size = grid_size / 2;
    
    // Calculate frequency components - parallelized first part
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            for (int i = 0; i <= half_size; i++) {
                qvecs->qx[i] = i * factor;
            }
            for (int i = half_size + 1; i < grid_size; i++) {
                qvecs->qx[i] = (i - grid_size) * factor;
            }
        }
        
        #pragma omp section
        {
            for (int i = 0; i <= half_size; i++) {
                qvecs->qy[i] = i * factor;
            }
            for (int i = half_size + 1; i < grid_size; i++) {
                qvecs->qy[i] = (i - grid_size) * factor;
            }
        }
        
        #pragma omp section
        {
            for (int i = 0; i <= half_size; i++) {
                qvecs->qz[i] = i * factor;
            }
            for (int i = half_size + 1; i < grid_size; i++) {
                qvecs->qz[i] = (i - grid_size) * factor;
            }
        }
    }
    
    return qvecs;
}

// Free q-vectors
void free_q_vectors(QVectors* qvecs) {
    if (qvecs) {
        free(qvecs->qx);
        free(qvecs->qy);
        free(qvecs->qz);
        free(qvecs);
    }
}