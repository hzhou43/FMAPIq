#ifndef SAXS_FUNCTIONS_H
#define SAXS_FUNCTIONS_H

#include <sys/time.h>
#include "saxs_calculator.h"
#include "pqr_reader.h"

// Add source grid values to target grid with scaling
void add2grid(float* source, float* target, int grid_size, double scale);

// Add uniform values to target grid for grid points above threshold
void addU2grid(float* source, float* target, int grid_size, double grid_spacing, double scale);

// Calculate B-factor from displacement
double u2B(double u);

// Get elapsed time in seconds
double get_elapsed_seconds(struct timeval start_time);

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
);

#endif // SAXS_FUNCTIONS_H