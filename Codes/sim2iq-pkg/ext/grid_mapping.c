#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "grid_mapping.h"
#include "rdf.h"
#include "cromer_mann_4gauss.h"

// Calculate the radius of a sphere from its volume
double sphere_radius_from_volume(double volume) {
    return pow((3.0 * volume) / (4.0 * M_PI), 1.0 / 3.0);
}

// Map a single atom to a small grid covering only the cutoff region
int grid_for_one_atom(
    size_t* full_indices,
    double* r2_values,
    int small_size,
    int offsets[3],
    int grid_size,
    double grid_spacing,
    const Atom* atom,
    double cutoff_radius
) {
    double cutoff_radius_squared = cutoff_radius * cutoff_radius;
    
    // Calculate atom position in grid coordinates (relative to full grid)
    double atom_x = atom->x / grid_spacing;
    double atom_y = atom->y / grid_spacing;
    double atom_z = atom->z / grid_spacing;
    
    // Calculate atom position relative to small grid
    double rel_x = atom_x - offsets[0];
    double rel_y = atom_y - offsets[1];
    double rel_z = atom_z - offsets[2];
    
    int valid_count = 0;
    
    // Calculate valid points and their distances
    for (int x = 0; x < small_size; x++) {
        double dx = (x - rel_x) * grid_spacing;
        
        for (int y = 0; y < small_size; y++) {
            double dy = (y - rel_y) * grid_spacing;
            
            for (int z = 0; z < small_size; z++) {
                double dz = (z - rel_z) * grid_spacing;
                
                // Calculate squared distance
                double r2 = dx*dx + dy*dy + dz*dz;
                
                // Only process points within cutoff
                if (r2 <= cutoff_radius_squared) {
                    // Map to full grid with PBC
                    int full_x = (offsets[0] + x + grid_size) % grid_size;
                    int full_y = (offsets[1] + y + grid_size) % grid_size;
                    int full_z = (offsets[2] + z + grid_size) % grid_size;
                    
                    // Calculate full grid index
                    full_indices[valid_count] = (size_t)full_x + (size_t)full_y * grid_size + 
                                               (size_t)full_z * (size_t)grid_size * grid_size;
                    
                    // Store squared distance
                    r2_values[valid_count] = r2;
                    
                    valid_count++;
                }
            }
        }
    }
    
    return valid_count;
}

// Map all atoms to the grid with normalization
void map_atoms_to_grid(
    float* grid,
    int grid_size,
    double grid_spacing,
    const Structure* structure,
    double cutoff_offset,
    const RdfParams* rdf_params,
    ElementTypeFunc get_element_type
) { 
    // Initialize the grid to zeros
    size_t total_size = (size_t)grid_size * (size_t)grid_size * (size_t)grid_size;
    #pragma omp parallel for simd
    for (size_t i = 0; i < total_size; i++) {
        grid[i] = 0.0f;
    }
    
    // Process atoms in parallel
    #pragma omp parallel for schedule(guided, 16)
    for (size_t i = 0; i < structure->count; i++) {
        const Atom* atom = &structure->atoms[i];
        
        // Get element type and parameters
        int element_type = get_element_type(atom->element);
        const CromerMann4Gauss* cm = get_cromer_mann_params(element_type);
        
        // Determine cutoff radius based on normalization mode
        double cutoff_radius;
        switch (rdf_params->norm_mode) {
            case NORM_ELECTRON_COUNT:
                cutoff_radius = 3.0;
                break;
            case NORM_VOLUME:
                cutoff_radius = 2.0 * cm->vdw_radius; // 2 * pdb.vdW
                break;
            case NORM_DEFAULT:
                cutoff_radius = 3.0;
                break;
            default:
                cutoff_radius = cutoff_offset + cm->vdw_radius;
        }
        
        // Calculate small grid size based on cutoff
        int small_size = 2 * (int)ceil(cutoff_radius / grid_spacing) + 1;
        
        // Allocate memory for this atom
        int max_points = small_size * small_size * small_size;
        size_t* full_indices = (size_t*)malloc(max_points * sizeof(size_t));
        double* r2_values = (double*)malloc(max_points * sizeof(double));
        double* contributions = (double*)malloc(max_points * sizeof(double));
        
        if (!full_indices || !r2_values || !contributions) {
            fprintf(stderr, "Error: Failed to allocate memory for atom %zu\n", i);
            if (full_indices) free(full_indices);
            if (r2_values) free(r2_values);
            if (contributions) free(contributions);
            continue;  // Skip this atom
        }
        
        // Calculate atom position in centered grid coordinates
        int atom_x = (int)floor(atom->x / grid_spacing);
        int atom_y = (int)floor(atom->y / grid_spacing);
        int atom_z = (int)floor(atom->z / grid_spacing);
        
        // Calculate offsets (relative to small grid)
        int offsets[3] = {
            atom_x - small_size/2,
            atom_y - small_size/2,
            atom_z - small_size/2
        };
        
        // Get valid grid points and their distances
        int num_valid = grid_for_one_atom(
            full_indices,
            r2_values,
            small_size,
            offsets,
            grid_size,
            grid_spacing,
            atom,
            cutoff_radius
        );
        
        if (num_valid > 0) {
            // Calculate contributions for all valid points (vectorized)
            double total_contribution = 0.0;
            
            // Vectorized calculation of contributions
            #pragma omp simd reduction(+:total_contribution)
            for (int j = 0; j < num_valid; j++) {
                double r2 = r2_values[j];
                double contribution;
                
                // Calculate RDF contribution
                switch (rdf_params->norm_mode) {
                    case NORM_VOLUME:
                        contribution = simple_gauss_rdf(r2, atom->vol, rdf_params->rho0);
                        break;
                    case NORM_NONE:
                        contribution = 1.0;
                        break;
                    default:
                        contribution = multigauss_rdf(r2, element_type, rdf_params->global_b);
                        break;
                }
                
                // Store contribution
                contributions[j] = contribution;
                total_contribution += contribution;
            }
            
            // Determine normalization factor based on mode
            double normalization_factor = 1.0;
            if (total_contribution > 0.0) {
                switch (rdf_params->norm_mode) {
                    case NORM_ELECTRON_COUNT:
                        if (cm && cm->electrons > 0) {
                            normalization_factor = cm->electrons / total_contribution;
                        } else {
                            normalization_factor = 1.0 / total_contribution;
                        }
                        break;
                    case NORM_VOLUME:
                        normalization_factor = atom->vol * rdf_params->rho0 / total_contribution;
                        break;
                    case NORM_DEFAULT:
                        normalization_factor = 1.0 / total_contribution;
                        break;
                    default:
                        normalization_factor = 1.0;
                }
                
                // Apply normalization and add to main grid
                for (int j = 0; j < num_valid; j++) {
                    double normalized_contribution = contributions[j] * normalization_factor;
                    
                    if (normalized_contribution > 0.0) {
                        // Add to full grid with atomic operation
                        #pragma omp atomic
                        grid[full_indices[j]] += (float)normalized_contribution;
                    }
                }
            }
        }
        
        // Free memory for this atom
        free(full_indices);
        free(r2_values);
        free(contributions);
    } 
}