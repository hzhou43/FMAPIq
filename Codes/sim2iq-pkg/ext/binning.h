#ifndef BINNING_H
#define BINNING_H

#include <fftw3.h>
#include "saxs_calculator.h"

// Create binning structure
BinData* create_bin_data(const QVectors* qvecs);

// Free binning structure
void free_bin_data(BinData* bins);

// Calculate intensity at q=0
void calculate_I0(const fftwf_complex* ft_out, BinData* bins);

// Accumulate intensity in bins with separate vacuum and other contributions
void accumulate_intensity(
    const fftwf_complex* ft_vac, 
    const fftwf_complex* ft_other, 
    const QVectors* qvecs, 
    BinData* bins,
    double b_factor
);

// Create final I(q) data from bins
QIData create_final_data(const BinData* bins);

// Write I(q) data to file
int write_output(const QIData* data, const char* filename);

#endif /* BINNING_H */