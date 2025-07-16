#ifndef ARGS_H
#define ARGS_H

#include "rdf.h"
#include <stdbool.h>

// Structure to store command line arguments
typedef struct {
    const char *input_file;    // Input PQR file
    const char *output_file;   // Output data file
    const char *traj_file;     // Trajectory file (optional)
    int grid_size;             // Size of grid in each dimension
    double grid_spacing;       // Spacing between grid points (Å)
    double cutoff_offset;      // Offset for cutoff radius (Å)
    double b_factor;           // Global B-factor
    int steps;                 // Number of calculation steps to perform
    double se_weight;          // Weight for solvent excluded volume step
    double ou_weight;          // Weight for outer shell step
    double in_weight;          // Weight for inner shell step
    bool has_traj;             // Flag to indicate if trajectory file is provided
    bool traj_is_txt;          // Flag to indicate if trajectory file is in text format
} CommandArgs;

// Parse command line arguments
int parse_args(int argc, char *argv[], CommandArgs *args);

#endif /* ARGS_H */