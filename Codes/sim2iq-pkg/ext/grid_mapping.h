#ifndef GRID_MAPPING_H
#define GRID_MAPPING_H

#include "saxs_calculator.h"
#include "rdf.h"

// Calculate the radius of a sphere from its volume
double sphere_radius_from_volume(double volume);

typedef int (*ElementTypeFunc)(const char* element);

// Map all atoms to the grid with normalization
void map_atoms_to_grid(
    float* grid,
    int grid_size,
    double grid_spacing,
    const Structure* structure,
    double cutoff_offset,
    const RdfParams* rdf_params,
    ElementTypeFunc get_element_type
);

#endif // GRID_MAPPING_H