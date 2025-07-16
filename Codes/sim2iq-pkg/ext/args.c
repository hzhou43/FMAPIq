#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "args.h"

// Parse command line arguments
int parse_args(int argc, char *argv[], CommandArgs *args) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <input.pqr> <grid_size> <grid_spacing> <output.dat> [options]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -b <value>        B-factor (default: 0.0)\n");
        fprintf(stderr, "  -offset <value>   Cutoff offset for outer shell (default: 2.8)\n");
        fprintf(stderr, "  -steps <value>    Number of calculation steps to perform (default: 1)\n");
        fprintf(stderr, "  -se <value>       Solvent excluded volume weight (default: -1.0)\n");
        fprintf(stderr, "  -ou <value>       Outer shell weight (default: 0.011)\n");
        fprintf(stderr, "  -in <value>       Inner shell weight (default: -0.011)\n");
        fprintf(stderr, "  -traj <file>      PD6D trajectory file (optional)\n");
        fprintf(stderr, "  -txt              Indicate trajectory file is in text format (default: PDB format)\n");
        return 0;
    }
    
    // Required arguments
    args->input_file = argv[1];
    args->grid_size = atoi(argv[2]);
    args->grid_spacing = atof(argv[3]);
    args->output_file = argv[4];
    
    // Set defaults for optional arguments
    args->b_factor = 0.0;
    args->cutoff_offset = 2.8;
    args->steps = 1;
    args->se_weight = -1.0;
    args->ou_weight = 0.011;
    args->in_weight = -0.011;
    args->has_traj = false;
    args->traj_is_txt = false;
    args->traj_file = NULL;
    
    // Parse optional arguments
    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            args->b_factor = atof(argv[++i]);
        } else if (strcmp(argv[i], "-offset") == 0 && i + 1 < argc) {
            args->cutoff_offset = atof(argv[++i]);
        } else if (strcmp(argv[i], "-steps") == 0 && i + 1 < argc) {
            args->steps = atoi(argv[++i]);
            if (args->steps < 1 || args->steps > 4) {
                fprintf(stderr, "Error: Steps must be between 1 and 4\n");
                return 0;
            }
        } else if (strcmp(argv[i], "-se") == 0 && i + 1 < argc) {
            args->se_weight = atof(argv[++i]);
        } else if (strcmp(argv[i], "-ou") == 0 && i + 1 < argc) {
            args->ou_weight = atof(argv[++i]);
        } else if (strcmp(argv[i], "-in") == 0 && i + 1 < argc) {
            args->in_weight = atof(argv[++i]);
        } else if (strcmp(argv[i], "-traj") == 0 && i + 1 < argc) {
            args->traj_file = argv[++i];
            args->has_traj = true;
        } else if (strcmp(argv[i], "-txt") == 0) {
            args->traj_is_txt = true;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 0;
        }
    }
    
    // Validate arguments
    if (args->grid_size <= 0) {
        fprintf(stderr, "Error: Grid size must be positive\n");
        return 0;
    }
    
    if (args->grid_spacing <= 0.0) {
        fprintf(stderr, "Error: Grid spacing must be positive\n");
        return 0;
    }
    
    if (args->has_traj && args->traj_file == NULL) {
        fprintf(stderr, "Error: Trajectory file must be specified with -traj option\n");
        return 0;
    }
    
    // Print parameters
    printf("Input file: %s\n", args->input_file);
    printf("Grid size: %d\n", args->grid_size);
    printf("Grid spacing: %.6f Å\n", args->grid_spacing);
    printf("B-factor: %.6f\n", args->b_factor);
    printf("Cutoff offset: %.6f Å\n", args->cutoff_offset);
    printf("Number of steps: %d\n", args->steps);
    printf("Solvent excluded volume weight: %.6f\n", args->se_weight);
    printf("Outer shell weight: %.6f\n", args->ou_weight);
    printf("Inner shell weight: %.6f\n", args->in_weight);
    printf("Output file: %s\n", args->output_file);
    
    if (args->has_traj) {
        printf("Trajectory file: %s\n", args->traj_file);
        printf("Trajectory format: %s\n", args->traj_is_txt ? "Text" : "PDB");
    }
    
    return 1;
}