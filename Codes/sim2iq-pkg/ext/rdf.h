#ifndef RDF_H
#define RDF_H

#include <stdbool.h>

// Normalization modes for RDF density mapping
typedef enum {
    NORM_DEFAULT = 0,       // Default normalization (sum to 1.0)
    NORM_ELECTRON_COUNT = 1, // Normalize by electron count
    NORM_VOLUME = 2,         // Normalize by volume
    NORM_NONE = 3            // No normalization (use RDF values directly)
} NormalizationMode;

// Parameters for RDF calculations
typedef struct {
    double global_b;        // Global B-factor
    double vol;             // Volume for simple Gaussian model
    double rho0;            // Electron density for simple Gaussian model
    NormalizationMode norm_mode; // Normalization mode for grid mapping
} RdfParams;

// Calculate Multi-Gaussian RDF based on Cromer-Mann parameters
double multigauss_rdf(double r2, int element_type, double B);

// Calculate Simple Gaussian RDF
double simple_gauss_rdf(double r2, double V, double rho0);

#endif /* RDF_H */