#include <math.h>
#include <stdbool.h>
#include "cromer_mann_4gauss.h"
#include "rdf.h"

// Calculate Multi-Gaussian RDF based on Cromer-Mann parameters
double multigauss_rdf(double r2, int element_type, double B) {
    const CromerMann4Gauss* cm = get_cromer_mann_params(element_type);
    if (!cm) return 0.0;  // Unknown element
    
    double ff = 0.0;
    
    for (int i = 0; i < 4; i++) {
        double B_term = cm->b[i] + B;
        double prefactor = pow(4 * M_PI / B_term, 3.0/2.0) * cm->a[i];
        ff += prefactor * exp(-4 * M_PI * M_PI * r2 / B_term);
    }
    
    return ff;
}

// Calculate Simple Gaussian RDF
double simple_gauss_rdf(double r2, double V, double rho0) {
    
    double B_term = 4 * M_PI * pow(V, 2.0/3.0);
    double prefactor = 8 * rho0 * pow(M_PI, 3.0/2.0) * V / pow(B_term, 3.0/2.0);
    return prefactor * exp(-4 * M_PI * M_PI * r2 / B_term);
}