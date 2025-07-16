#ifndef FFT_H
#define FFT_H

#include <fftw3.h>
#include "saxs_calculator.h"

// Create FFT grid from real-space grid
fftwf_complex* calculate_fft(fftwf_complex* fft_grid, int grid_size);

#endif /* FFT_H */