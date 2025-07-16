#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>

#include "saxs_calculator.h"
#include "args.h"
#include "pqr_reader.h"
#include "grid_mapping.h"
#include "rdf.h"
#include "fft.h"
#include "binning.h"
#include "cromer_mann_4gauss.h"
#include "saxs_functions.h"

// Add values from source grid to target grid with scaling
void add2grid(float* source, float* target, int grid_size, double scale) {
    size_t total_size = (size_t)grid_size * (size_t)grid_size * (size_t)grid_size;
    
    #pragma omp parallel for simd
    for (size_t i = 0; i < total_size; i++) {
        target[i] += source[i] * scale;
    }
}

// Add uniform values to target grid for grid points above threshold
void addU2grid(float* source, float* target, int grid_size, double grid_spacing, double scale) {
    const double dV = grid_spacing * grid_spacing * grid_spacing;
    size_t total_size = (size_t)grid_size * (size_t)grid_size * (size_t)grid_size;
    double total_volume = 0.0;
    
    #pragma omp parallel for reduction(+:total_volume)
    for (size_t i = 0; i < total_size; i++) {
        if (source[i] > 0.5) {
            target[i] += scale*dV;
            total_volume += dV;
        }
    }

    printf("Volume: %.12e\n", total_volume);
}

// Calculate B-factor from displacement
double u2B(double u) {
    return 8.0 * M_PI * M_PI * u * u;
}

// Get elapsed time in seconds since start time
double get_elapsed_seconds(struct timeval start_time) {
    struct timeval current_time;
    gettimeofday(&current_time, NULL);
    
    return (double)(current_time.tv_sec - start_time.tv_sec) + 
           (double)(current_time.tv_usec - start_time.tv_usec) / 1000000.0;
}

// Multi-step calculation of SAXS data
QIData calculate_saxs_data_multi_step(
    Structure* structure,
    int grid_size,
    double grid_spacing,
    int steps,
    double cutoff_offset,
    double se_weight,
    double ou_weight,
    double in_weight,
    double b_factor
) {
    QIData result = {0};
    struct timeval start_time;
    gettimeofday(&start_time, NULL);

    // Allocate grids - one for vacuum, one for other contributions, and one work grid
    size_t grid_total = (size_t)grid_size * (size_t)grid_size * (size_t)grid_size;
    float* vacuum_grid = (float*)malloc(grid_total * sizeof(float));
    float* other_grid = (float*)malloc(grid_total * sizeof(float));
    float* work_grid = (float*)malloc(grid_total * sizeof(float));

    if (!vacuum_grid || !other_grid || !work_grid) {
        fprintf(stderr, "Error: Failed to allocate memory for grids\n");
        if (vacuum_grid) free(vacuum_grid);
        if (other_grid) free(other_grid);
        if (work_grid) free(work_grid);
        return result;
    }

    // Initialize grids to zeros
    memset(vacuum_grid, 0, grid_total * sizeof(float));
    memset(other_grid, 0, grid_total * sizeof(float));
    memset(work_grid, 0, grid_total * sizeof(float));

    printf("[%.3f s] Mapping atoms to grid...\n", get_elapsed_seconds(start_time));

    // Step 1: Map atoms to grid with electron count normalization (vacuum contribution)
    RdfParams rdf_params = {
        .global_b = b_factor,
        .rho0 = 0.334,
        .norm_mode = NORM_ELECTRON_COUNT
    };

    map_atoms_to_grid(
        vacuum_grid,
        grid_size,
        grid_spacing,
        structure,
        0.0, // Default cutoff
        &rdf_params,
        get_element_type
    );

    // Step 2: Map atoms to grid with volume normalization for solvent excluded volume
    if (steps > 1) {
        printf("[%.3f s] Mapping solvent excluded volume...\n", get_elapsed_seconds(start_time));

        memset(work_grid, 0, grid_total * sizeof(float));

        rdf_params.norm_mode = NORM_VOLUME;

        map_atoms_to_grid(
            work_grid,
            grid_size,
            grid_spacing,
            structure,
            0.0, // Default cutoff
            &rdf_params,
            get_element_type
        );

        add2grid(work_grid, other_grid, grid_size, se_weight);
    }

    // Step 3: Map outer shell with no normalization (goes to other grid)
    if (steps > 2) {
        printf("[%.3f s] Mapping outer shell...\n", get_elapsed_seconds(start_time));

        memset(work_grid, 0, grid_total * sizeof(float));

        rdf_params.norm_mode = NORM_NONE;

        map_atoms_to_grid(
            work_grid,
            grid_size,
            grid_spacing,
            structure,
            cutoff_offset, // Use cutoff offset for outer shell
            &rdf_params,
            get_element_type
        );

        addU2grid(work_grid, other_grid, grid_size, grid_spacing, ou_weight);
    }

    // Step 4: Map inner shell with no normalization (goes to other grid)
    if (steps > 3) {
        printf("[%.3f s] Mapping inner shell...\n", get_elapsed_seconds(start_time));

        memset(work_grid, 0, grid_total * sizeof(float));

        rdf_params.norm_mode = NORM_NONE;

        map_atoms_to_grid(
            work_grid,
            grid_size,
            grid_spacing,
            structure,
            0.0, // No offset for inner shell
            &rdf_params,
            get_element_type
        );

        addU2grid(work_grid, other_grid, grid_size, grid_spacing, in_weight);
    }

    free(work_grid);

    printf("[%.3f s] Calculating FFTs...\n", get_elapsed_seconds(start_time));

    // Create q-vectors for binning
    QVectors* qvecs = create_q_vectors(grid_size, grid_spacing);
    if (!qvecs) {
        fprintf(stderr, "Error: Failed to create q-vectors\n");
        free(vacuum_grid);
        free(other_grid);
        return result;
    }

    // Allocate FFT grids
    fftwf_complex* fft_vacuum = fftwf_alloc_complex(grid_total);
    fftwf_complex* fft_other = NULL;
    
    if (!fft_vacuum) {
        fprintf(stderr, "Error: Failed to allocate memory for FFT grids\n");
        free(vacuum_grid);
        free(other_grid);
        free_q_vectors(qvecs);
        return result;
    }
    
    // Copy vacuum grid to complex grid
    #pragma omp parallel for simd
    for (size_t i = 0; i < grid_total; i++) {
        fft_vacuum[i][0] = vacuum_grid[i];
        fft_vacuum[i][1] = 0.0f;
    }
    
    free(vacuum_grid);
    
    // Perform FFT on vacuum grid
    if (!calculate_fft(fft_vacuum, grid_size)) {
        fprintf(stderr, "Error: Failed to calculate vacuum FFT\n");
        fftwf_free(fft_vacuum);
        free(other_grid);
        free_q_vectors(qvecs);
        return result;
    }
    
    // Only allocate and calculate FFT for other grid if we have steps > 1
    if (steps > 1) {
        fft_other = fftwf_alloc_complex(grid_total);
        if (!fft_other) {
            fprintf(stderr, "Error: Failed to allocate memory for other FFT grid\n");
            fftwf_free(fft_vacuum);
            free(other_grid);
            free_q_vectors(qvecs);
            return result;
        }
        
        // Copy other grid to complex grid
        #pragma omp parallel for simd
        for (size_t i = 0; i < grid_total; i++) {
            fft_other[i][0] = other_grid[i];
            fft_other[i][1] = 0.0f;
        }
        
        free(other_grid);
        
        // Perform FFT on other grid
        if (!calculate_fft(fft_other, grid_size)) {
            fprintf(stderr, "Error: Failed to calculate other FFT\n");
            fftwf_free(fft_vacuum);
            fftwf_free(fft_other);
            free_q_vectors(qvecs);
            return result;
        }
    } else {
        free(other_grid);
    }

    printf("[%.3f s] Binning intensity...\n", get_elapsed_seconds(start_time));

    // Create binning structure
    BinData* bins = create_bin_data(qvecs);
    if (!bins) {
        free_q_vectors(qvecs);
        fftwf_free(fft_vacuum);
        if (fft_other) fftwf_free(fft_other);
        fprintf(stderr, "Error: Failed to create binning structure\n");
        return result;
    }

    // Accumulate intensity in bins
    accumulate_intensity(fft_vacuum, fft_other, qvecs, bins, b_factor);

    // Free FFT results
    fftwf_free(fft_vacuum);
    if (fft_other) fftwf_free(fft_other);

    // Create final I(q) data
    result = create_final_data(bins);

    printf("[%.3f s] SAXS calculation complete.\n", get_elapsed_seconds(start_time));

    // Free bin data
    free_bin_data(bins);
    free_q_vectors(qvecs);

    return result;
}