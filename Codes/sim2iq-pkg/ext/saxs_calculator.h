#ifndef SAXS_CALCULATOR_H
#define SAXS_CALCULATOR_H

#include <stddef.h>
#include <string.h>

// Use double precision FFTW
#define FFTW_FLOAT
#include <fftw3.h>

#define MAX_ELEMENT_LENGTH 8
#define INITIAL_CAPACITY 1024
#define MAX_LINE_LENGTH 256

#ifndef M_PI
#   define M_PI 3.1415926535897932384626433832
#endif

// Basic structures for atomic data
typedef struct {
    char element[MAX_ELEMENT_LENGTH];
    double x, y, z, vol;
} Atom;

typedef struct {
    Atom* atoms;
    size_t count;
    size_t capacity;
} Structure;

typedef struct {
    double q;
    double intensity;
} QIPoint;

typedef struct {
    QIPoint* data;
    size_t count;
    size_t capacity;
} QIData;

// Q-vector management structures
typedef struct {
    double* qx;
    double* qy;
    double* qz;
    int box[3];
    double spacing;
    size_t total;
} QVectors;

// Binning data structure
typedef struct {
    double* sums;           // Sum of intensities in each bin
    size_t* counts;        // Count of points in each bin
    double* q_sums;         // Sum of q values in each bin
    size_t max_bins;
    double q_bin_width;
    double qmax;
    double Iq0;
} BinData;

// Memory management utilities
void* check_malloc(size_t size);
void qidata_free(QIData* data);
QIData qi_data_create(size_t initial_capacity);
Structure structure_create(void);
void structure_ensure_capacity(Structure* structure);
void structure_free(Structure* structure);

// Q-vector management functions
QVectors* create_q_vectors(int grid_size, double spacing);
void free_q_vectors(QVectors* qvecs);

#endif // SAXS_CALCULATOR_H