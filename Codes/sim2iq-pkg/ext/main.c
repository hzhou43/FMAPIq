#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <stdbool.h>

#include "saxs_calculator.h"
#include "args.h"
#include "pqr_reader.h"
#include "saxs_functions.h"
#include "pd6d.h"
#include "binning.h"

// Apply trajectory transformation to a structure
Structure process_trajectory(Structure* original_structure, const char* traj_file, bool is_txt) {
    // Count atoms in the trajectory file
    int natoms = PD6DCountAtoms((char*)traj_file, is_txt);
    
    // Allocate memory for the transformation and rotation data
    double (*tr)[3] = malloc(natoms * sizeof(*tr));
    double (*xr)[3] = malloc(natoms * sizeof(*xr));
    
    if (!tr || !xr) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        if (tr) free(tr);
        if (xr) free(xr);
        exit(EXIT_FAILURE);
    }
    
    // Read the transformation data
    FILE* file = fopen(traj_file, "r");
    int read_atoms = PD6DReadxyz(file, natoms, tr, xr, is_txt);
    fclose(file);
    
    if (read_atoms != natoms) {
        fprintf(stderr, "Warning: Read %d proteins from trajectory file, expected %d\n", read_atoms, natoms);
    }
    
    printf("Successfully read %d proteins from trajectory file\n", read_atoms);
    
    // Create a new structure with multiple copies of the original
    int orig_num_atoms = original_structure->count;
    int total_atoms = orig_num_atoms * natoms;
    
    Structure new_structure;
    new_structure.count = total_atoms;
    new_structure.atoms = malloc(total_atoms * sizeof(Atom));
    
    if (!new_structure.atoms) {
        fprintf(stderr, "Error: Memory allocation failed for new structure\n");
        free(tr);
        free(xr);
        exit(EXIT_FAILURE);
    }
    
    // Process each transformation separately in parallel
    #pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < natoms; t++) {
        // Thread-local arrays for coordinates
        double (*xyz_old)[3] = malloc(orig_num_atoms * sizeof(*xyz_old));
        double (*xyz_new)[3] = malloc(orig_num_atoms * sizeof(*xyz_new));
        
        if (!xyz_old || !xyz_new) {
            if (xyz_old) free(xyz_old);
            if (xyz_new) free(xyz_new);
            continue; // Skip this iteration but continue with others
        }
        
        // Fill old coordinates array with original structure
        for (int i = 0; i < orig_num_atoms; i++) {
            xyz_old[i][0] = original_structure->atoms[i].x;
            xyz_old[i][1] = original_structure->atoms[i].y;
            xyz_old[i][2] = original_structure->atoms[i].z;
        }
        
        // Apply transformation to the copy - thread-safe operation
        pd6dmmTR(tr[t], xr[t], orig_num_atoms, xyz_old, xyz_new);
        
        // Copy the transformed coordinates to the new structure at the appropriate offset
        int offset = t * orig_num_atoms;
        for (int i = 0; i < orig_num_atoms; i++) {
            // Copy atom properties from original
            new_structure.atoms[offset + i] = original_structure->atoms[i];
            
            // Update with new coordinates
            new_structure.atoms[offset + i].x = xyz_new[i][0];
            new_structure.atoms[offset + i].y = xyz_new[i][1];
            new_structure.atoms[offset + i].z = xyz_new[i][2];
        }
        
        // Clean up arrays for this transformation
        free(xyz_old);
        free(xyz_new);
    }
    
    // Clean up
    free(tr);
    free(xr);
    
    printf("Trajectory processing complete - created structure with %ld atoms\n", new_structure.count);
    return new_structure;
}

int main(int argc, char *argv[]) {
    // Step 0: Process command-line arguments
    CommandArgs args;
    if (!parse_args(argc, argv, &args)) {
        return EXIT_FAILURE;
    }

    // Step 1: Read PQR file
    Structure original_structure = read_pqr_file(args.input_file);

    // Step 2: Process trajectory if provided
    Structure* structure_to_use = &original_structure;
    Structure transformed_structure;

    if (args.has_traj) {
        printf("Processing trajectory file: %s\n", args.traj_file);
        transformed_structure = process_trajectory(&original_structure, args.traj_file, args.traj_is_txt);
        structure_to_use = &transformed_structure;
        printf("Using transformed structure for calculations\n");
    }

    // Calculate B-factor based on grid spacing
    double b_factor = u2B(0.25 * args.grid_spacing);
    if (args.b_factor > 0.0) {
        b_factor = args.b_factor;
    }
    printf("Actual B-factor: %.6f\n", b_factor);

    // Calculate SAXS data
    QIData result = calculate_saxs_data_multi_step(
        structure_to_use,
        args.grid_size,
        args.grid_spacing,
        args.steps,
        args.cutoff_offset,
        args.se_weight,
        args.ou_weight,
        args.in_weight,
        b_factor
    );

    // Free structures as they're no longer needed
    if (args.has_traj) {
        structure_free(&transformed_structure);
    }
    structure_free(&original_structure);

    // Check if result is valid
    if (result.count == 0) {
        fprintf(stderr, "Error: No data points generated\n");
        qidata_free(&result);
        return EXIT_FAILURE;
    }

    if (!write_output(&result, args.output_file)) {
        fprintf(stderr, "Error: Failed to write output file\n");
        qidata_free(&result);
        return EXIT_FAILURE;
    }

    // Free result data
    qidata_free(&result);

    return EXIT_SUCCESS;
}