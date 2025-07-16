#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "binning.h"

// Create binning structure
BinData* create_bin_data(const QVectors* qvecs) {
    // Calculate qmax based on maximum frequency
    const double qmax = qvecs->qx[qvecs->box[0]/2-1] - 1e-8; 
    
    // Calculate q_bin_width
    const double q_bin_width = 2.0 * (double)M_PI / (qvecs->box[0] * qvecs->spacing) - 1e-8;
    
    // Calculate max_bins
    const size_t max_bins = (size_t)(qmax / q_bin_width) + 1;
    
    // Allocate bin data
    BinData* bins = malloc(sizeof(BinData));
    if (!bins) {
        fprintf(stderr, "Error: Failed to allocate memory for bin data\n");
        return NULL;
    }
    
    // Initialize bin data
    bins->sums = calloc(max_bins, sizeof(double));
    bins->counts = calloc(max_bins, sizeof(size_t));
    bins->q_sums = calloc(max_bins, sizeof(double));
    bins->Iq0 = 0.0;
    
    if (!bins->sums || !bins->counts || !bins->q_sums) {
        fprintf(stderr, "Error: Failed to allocate memory for bin arrays\n");
        free_bin_data(bins);
        return NULL;
    }
    
    bins->max_bins = max_bins;
    bins->q_bin_width = q_bin_width;
    bins->qmax = qmax;
    
    return bins;
}

// Free binning structure
void free_bin_data(BinData* bins) {
    if (bins) {
        free(bins->sums);
        free(bins->counts);
        free(bins->q_sums);
        free(bins);
    }
}

// Calculate intensity at q=0
void calculate_I0(const fftwf_complex* ft_out, BinData* bins) {
    // Calculate intensity at q=0 (origin)
    const double re = ft_out[0][0];
    const double im = ft_out[0][1];
    bins->Iq0 = re*re + im*im;
}

// Accumulate intensity in bins with separate vacuum and other contributions
void accumulate_intensity(const fftwf_complex* ft_vac, const fftwf_complex* ft_other, const QVectors* qvecs, BinData* bins, double b_factor) {
    // Calculate I(q=0) from the sum of vacuum and other contributions
    fftwf_complex sum_at_origin;
    sum_at_origin[0] = ft_vac[0][0] + (ft_other ? ft_other[0][0] : 0.0f);
    sum_at_origin[1] = ft_vac[0][1] + (ft_other ? ft_other[0][1] : 0.0f);
    const double re = sum_at_origin[0];
    const double im = sum_at_origin[1];
    bins->Iq0 = re*re + im*im;
    
    // For each point in reciprocal space
    #pragma omp parallel
    {
        // Create thread-local bins to avoid contention
        double* local_sums = (double*)calloc(bins->max_bins, sizeof(double));
        double* local_q_sums = (double*)calloc(bins->max_bins, sizeof(double));
        size_t* local_counts = (size_t*)calloc(bins->max_bins, sizeof(size_t));
        
        if (local_sums && local_q_sums && local_counts) {
            #pragma omp for
            for (int k = 0; k < qvecs->box[2]; k++) {
                for (int j = 0; j < qvecs->box[1]; j++) {
                    for (int i = 0; i < qvecs->box[0]; i++) {
                        // Skip origin (already processed)
                        if (i == 0 && j == 0 && k == 0) continue;
                        
                        // Calculate q magnitude
                        const double q = sqrt(
                            qvecs->qx[i] * qvecs->qx[i] +
                            qvecs->qy[j] * qvecs->qy[j] +
                            qvecs->qz[k] * qvecs->qz[k]
                        );
                        
                        // Only process q values within range
                        if (q > 0 && q < bins->qmax) {
                            // Determine bin index
                            const size_t bin_idx = (size_t)(q / bins->q_bin_width);
                            
                            // Make sure bin index is valid
                            if (bin_idx < bins->max_bins) {
                                // Calculate index in FFT result
                                const size_t idx = i + (size_t)j * qvecs->box[0] + (size_t)k * qvecs->box[0] * qvecs->box[1];
                                
                                // Calculate B-factor sharpening (only for vacuum)
                                double bsharp3D = 1.0;
                                if (b_factor > 0.0) {
                                    bsharp3D = exp(b_factor * pow((q / (4.0 * M_PI)), 2));
                                }
                                
                                // Apply B-factor to vacuum term
                                double ft_vac_bsharp[2];
                                ft_vac_bsharp[0] = ft_vac[idx][0] * bsharp3D;
                                ft_vac_bsharp[1] = ft_vac[idx][1] * bsharp3D;
                                
                                // Sum vacuum and other (if present)
                                double re, im;
                                if (ft_other) {
                                    re = ft_vac_bsharp[0] + ft_other[idx][0];
                                    im = ft_vac_bsharp[1] + ft_other[idx][1];
                                } else {
                                    re = ft_vac_bsharp[0];
                                    im = ft_vac_bsharp[1];
                                }
                                
                                // Calculate intensity
                                const double intensity = re*re + im*im;
                                
                                // Add to local bins
                                local_sums[bin_idx] += intensity;
                                local_q_sums[bin_idx] += q;
                                local_counts[bin_idx]++;
                            }
                        }
                    }
                }
            }
            
            // Define reduction variables with direct assignment
            double* sums = bins->sums;
            double* q_sums = bins->q_sums;
            size_t* counts = bins->counts;
            size_t max_bins = bins->max_bins;

            // Use the reduction on these variables
            #pragma omp parallel for reduction(+:sums[:max_bins], q_sums[:max_bins], counts[:max_bins])
            for (size_t i = 0; i < max_bins; i++) {
                sums[i] += local_sums[i];
                q_sums[i] += local_q_sums[i];
                counts[i] += local_counts[i];
            }
        }
        
        // Free thread-local memory
        free(local_sums);
        free(local_q_sums);
        free(local_counts);
    }
}

// Create final I(q) data from bins
QIData create_final_data(const BinData* bins) {
    // Add an extra entry for q=0
    QIData result = qi_data_create(bins->max_bins + 1);
    if (!result.data) {
        fprintf(stderr, "Error: Failed to allocate memory for QI data\n");
        return result;
    }
    
    // Add I(q=0) point
    result.data[0].q = 0.0;
    result.data[0].intensity = bins->Iq0;
    result.count = 1;
    
    // Add binned points
    for (size_t i = 0; i < bins->max_bins; i++) {
        // Only include bins with data
        if (bins->counts[i] > 0) {
            // Calculate average intensity for this bin
            double avg_intensity = bins->sums[i] / bins->counts[i];
            
            // Calculate mean q value (qbinsc) for this bin
            double q_mean;
            if (i == 0) {
                // Special case for the first bin (q=0)
                q_mean = 0.0;
            } else {
                // Calculate mean q from accumulated q values
                q_mean = bins->q_sums[i] / bins->counts[i];
            }
            
            // Add to result
            result.data[result.count].q = q_mean;
            result.data[result.count].intensity = avg_intensity;
            result.count++;
        }
    }
    
    return result;
}

// Write I(q) data to file
int write_output(const QIData* data, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: Could not open output file: %s\n", filename);
        return 0;
    }

    // Write data
    for (size_t i = 0; i < data->count; i++) {
        fprintf(file, "%.12e\t%.12e\n",
                data->data[i].q,
                data->data[i].intensity);
    }

    fclose(file);
    printf("Wrote %zu data points to %s\n", data->count, filename);
    return 1;
}