// fmap_energy.c - Self-contained implementation of FMAP energy functions with consolidated conversions
#include "fmap_energy.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
// System V IPC headers (may need adjustment based on platform)
#ifdef __unix__
#include <sys/ipc.h>
#include <sys/shm.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Type definitions
typedef signed char tErn;
#define tErn_MAX SCHAR_MAX
#define tErn_MIN SCHAR_MIN
#define tErn_SCL 8.0

typedef short tAng;

// FMAP context structure
struct FMAP_Context {
    int nrot;           // Number of rotations
    int blen;           // Box length
    int nbin;           // Number of bins
    double scl;         // Scale factor
    double rcut2;       // Squared cutoff radius

    // Data pointers (from shared memory)
    tErn *energy_data;
    tAng *angle_indices;
    tAng *angle_mapping;
    double (*rotation_matrices)[3][3];
};

// ============================================================================
// COORDINATE CONVERSION FUNCTIONS (consolidated from other files)
// ============================================================================

// Convert quaternion to ZXZ Euler angles
void fmap_quat_to_euler(const double quat[4], double euler[3]) {
    const double qw = quat[0], qx = quat[1], qy = quat[2], qz = quat[3];
    
    double cos_middle = 2.0 * (qw*qw + qz*qz) - 1.0;
    cos_middle = fmax(-1.0, fmin(1.0, cos_middle)); // clamp
    
    const double epsilon = 2.2204460492503131e-15;
    
    if (fabs(cos_middle - 1.0) < epsilon) {
        euler[0] = 2.0 * atan2(qz, qw);
        euler[1] = 0.0;
        euler[2] = 0.0;
    }
    else if (fabs(cos_middle + 1.0) < epsilon) {
        euler[0] = 2.0 * atan2(qy, qx);
        euler[1] = M_PI;
        euler[2] = 0.0;
    }
    else {
        euler[0] = atan2(2.0 * (qw*qy + qx*qz), 2.0 * (qw*qx - qy*qz));
        euler[1] = acos(cos_middle);
        euler[2] = -atan2(2.0 * (qw*qy - qx*qz), 2.0 * (qw*qx + qy*qz));
    }
}

// Convert ZXZ Euler angles to quaternion
void fmap_euler_to_quat(const double euler[3], double quat[4]) {
    const double half_x = euler[0] * 0.5;
    const double half_y = euler[1] * 0.5;
    const double half_z = euler[2] * 0.5;
    
    const double cos_half_y = cos(half_y);
    const double sin_half_y = sin(half_y);
    
    const double sum_xz = half_x + half_z;
    const double diff_xz = half_x - half_z;
    
    const double cos_sum = cos(sum_xz);
    const double sin_sum = sin(sum_xz);
    const double cos_diff = cos(diff_xz);
    const double sin_diff = sin(diff_xz);
    
    quat[0] = cos_half_y * cos_sum;    // w component
    quat[1] = sin_half_y * cos_diff;   // x component
    quat[2] = sin_half_y * sin_diff;   // y component
    quat[3] = cos_half_y * sin_sum;    // z component
}

// Normalize quaternion
int fmap_quat_normalize(double quat[4]) {
    const double epsilon = 2.2204460492503131e-15;
    double mag = sqrt(quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3]);
    
    if (mag < epsilon) {
        return -1;  // Error: cannot normalize zero quaternion
    }
    
    const double inv_mag = 1.0 / mag;
    quat[0] *= inv_mag;
    quat[1] *= inv_mag;
    quat[2] *= inv_mag;
    quat[3] *= inv_mag;
    
    return 0;  // Success
}

// Convert Euler angles to bin indices
void fmap_euler_to_bins(const double euler[3], double bins[3], int nbin) {
    const double TWO_PI = 2.0 * M_PI;
    
    // Alpha
    double angle = fmod(euler[0] + M_PI, TWO_PI);
    if (angle < 0) angle += TWO_PI;
    bins[0] = (angle * nbin) / TWO_PI;
    
    // Beta (clamped to [0, PI])
    double beta = fmax(0.0, fmin(M_PI, euler[1]));
    bins[1] = ((1.0 - cos(beta)) * nbin) / 2.0;
    
    // Gamma
    angle = fmod(euler[2] + M_PI, TWO_PI);
    if (angle < 0) angle += TWO_PI;
    bins[2] = (angle * nbin) / TWO_PI;
}

// Convert bin indices to Euler angles
void fmap_bins_to_euler(const double bins[3], double euler[3], int nbin) {
    const double TWO_PI = 2.0 * M_PI;
    
    // Alpha
    euler[0] = -M_PI + (bins[0] / nbin) * TWO_PI;
    
    // Beta
    euler[1] = acos(1.0 - bins[1] / nbin * 2.0);
    
    // Gamma
    euler[2] = -M_PI + (bins[2] / nbin) * TWO_PI;
}

// Direct quaternion to bins conversion (optimized path)
void fmap_quat_to_bins(const double quat[4], double bins[3], int nbin) {
    double euler[3];
    fmap_quat_to_euler(quat, euler);
    fmap_euler_to_bins(euler, bins, nbin);
}

// Direct bins to quaternion conversion (optimized path)
void fmap_bins_to_quat(const double bins[3], double quat[4], int nbin) {
    double euler[3];
    fmap_bins_to_euler(bins, euler, nbin);
    fmap_euler_to_quat(euler, quat);
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

static inline double apply_pbc(double x, double box_length) {
    double offset = fmod(x, box_length);
    return (offset > box_length/2) ? offset - box_length : 
           (offset < -box_length/2) ? offset + box_length : offset;
}

static inline int to_grid_index(double x, int grid_size) {
    int idx = (int)round(fmod(x, grid_size));
    return (idx < 0) ? idx + grid_size : (idx == grid_size) ? 0 : idx;
}

static inline int orientation_to_index(const double orient[3], int nbin) {
    int i = to_grid_index(orient[0], nbin);
    int j = to_grid_index(orient[1], nbin);
    int k = to_grid_index(orient[2], nbin);
    return i * nbin * nbin + j * nbin + k;
}

static inline void apply_rotation(const double rotation[3][3], const double coords[3], double result[3]) {
    result[0] = rotation[0][0] * coords[0] + rotation[0][1] * coords[1] + rotation[0][2] * coords[2];
    result[1] = rotation[1][0] * coords[0] + rotation[1][1] * coords[1] + rotation[1][2] * coords[2];
    result[2] = rotation[2][0] * coords[0] + rotation[2][1] * coords[1] + rotation[2][2] * coords[2];
}

// ============================================================================
// KEY PARSING AND SHARED MEMORY ACCESS
// ============================================================================

#ifdef __unix__
// Helper function to parse key string to key_t
static key_t parse_key_string(const char* key_string) {
    if (!key_string) {
        fprintf(stderr, "Error: NULL key string provided\n");
        return -1;
    }
    
    char *endptr;
    long parsed_key = strtol(key_string, &endptr, 16);
    
    if (*endptr != '\0') {
        fprintf(stderr, "Error: Invalid key format: %s (contains non-hex characters)\n", key_string);
        return -1;
    }
    
    if (parsed_key < 0) {
        fprintf(stderr, "Error: Invalid key format: %s (negative value)\n", key_string);
        return -1;
    }
    
    // Apply the shift that was previously done in simulation.c
    key_t keyp = (key_t)(parsed_key << 4);
    
    return keyp;
}
#endif

// ============================================================================
// ENERGY CALCULATION FUNCTIONS
// ============================================================================

// Step 1: Calculate distance and relative coordinates (returns early if beyond cutoff)
static bool calculate_distance_and_coords(const FMAP_Context* ctx,
                                         const double pos1[3], const double pos2[3], 
                                         const double box[3], double rel_coords[3]) {
    // Calculate relative coordinates with PBC
    for (int i = 0; i < 3; i++) {
        rel_coords[i] = apply_pbc(pos1[i] - pos2[i], box[i]);
    }
    
    // Check distance cutoff
    double d2 = rel_coords[0]*rel_coords[0] + rel_coords[1]*rel_coords[1] + rel_coords[2]*rel_coords[2];
    if (d2 > ctx->rcut2) {
        return false; // Beyond cutoff
    }
    
    return true; // Within cutoff
}

// Step 2: Calculate energy from relative coordinates and orientations
static double calculate_energy_from_orientations(const FMAP_Context* ctx,
                                               const double rel_coords[3],
                                               const double orient1[3], const double orient2[3]) {
    // Get angular indices
    int a1 = ctx->angle_indices[orientation_to_index(orient1, ctx->nbin)];
    int a2 = ctx->angle_indices[orientation_to_index(orient2, ctx->nbin)];
    int a3 = ctx->angle_mapping[a1 * ctx->nrot + a2];
    
    // Apply rotation using first orientation
    double rotated_coords[3];
    apply_rotation((const double (*)[3])ctx->rotation_matrices[a1], rel_coords, rotated_coords);
    
    // Calculate grid indices
    int h = ctx->blen / 2;
    int i = (int)round(rotated_coords[0] * ctx->scl) + h;
    int j = (int)round(rotated_coords[1] * ctx->scl) + h;
    int k = (int)round(rotated_coords[2] * ctx->scl) + h;
    
    // Bounds check
    if (i < 0 || i >= ctx->blen || j < 0 || j >= ctx->blen || k < 0 || k >= ctx->blen) {
        return 0.0;
    }
    
    // Get energy value
    size_t idx = (size_t)a3 * ctx->blen * ctx->blen * ctx->blen + 
                 i * ctx->blen * ctx->blen + j * ctx->blen + k;
    tErn energy_raw = ctx->energy_data[idx];
    
    return (energy_raw == tErn_MAX) ? 9999.0 : (double)energy_raw / tErn_SCL;
}

// ============================================================================
// PUBLIC API FUNCTIONS
// ============================================================================

// Create FMAP context from shared memory key string
FMAP_Context* fmap_create(const char* key_string) {
#ifdef __unix__
    // Parse the key string
    key_t keyp = parse_key_string(key_string);
    if (keyp == -1) {
        return NULL;
    }
    
    //printf("FMAP: Parsed key string '%s' to key 0x%x\n", key_string, keyp);
    
    FMAP_Context* ctx = malloc(sizeof(FMAP_Context));
    if (!ctx) {
        fprintf(stderr, "Error: Failed to allocate FMAP context\n");
        return NULL;
    }
    
    // Access parameters shared memory
    int shmid = shmget(keyp, sizeof(double) * 10, 0666);
    if (shmid < 0) {
        fprintf(stderr, "Error: Failed to access parameters shared memory with key 0x%x\n", keyp);
        free(ctx);
        return NULL;
    }
    
    double *params = (double*)shmat(shmid, NULL, 0);
    if (params == (void*)-1) {
        fprintf(stderr, "Error: Failed to attach parameters shared memory\n");
        free(ctx);
        return NULL;
    }
    
    // Read parameters
    ctx->nrot = (int)params[0];
    ctx->blen = (int)params[1];
    ctx->nbin = (int)params[2];
    double scld = 1.0; //params[3];
    double spacing = params[4];
    int lsel = (int)params[5];
    
    ctx->scl = scld / spacing;
    double rcut = lsel / 2.0 / ctx->scl;
    ctx->rcut2 = rcut * rcut;
    //fprintf(stderr,"scl (%lf) rcut (%lf)\n",ctx->scl,rcut);
    
    shmdt(params);
    
    // Access energy data
    shmid = shmget(keyp + 1, sizeof(tErn) * (size_t)ctx->blen * ctx->blen * ctx->blen * ctx->nrot, 0666);
    if (shmid < 0) { 
        fprintf(stderr, "Error: Failed to access energy data shared memory\n");
        free(ctx); 
        return NULL; 
    }
    ctx->energy_data = (tErn*)shmat(shmid, NULL, 0);
    if (ctx->energy_data == (void*)-1) { 
        fprintf(stderr, "Error: Failed to attach energy data shared memory\n");
        free(ctx); 
        return NULL; 
    }
    
    // Access angle indices
    shmid = shmget(keyp + 2, sizeof(tAng) * (size_t)ctx->nbin * ctx->nbin * ctx->nbin, 0666);
    if (shmid < 0) { 
        fprintf(stderr, "Error: Failed to access angle indices shared memory\n");
        free(ctx); 
        return NULL; 
    }
    ctx->angle_indices = (tAng*)shmat(shmid, NULL, 0);
    if (ctx->angle_indices == (void*)-1) { 
        fprintf(stderr, "Error: Failed to attach angle indices shared memory\n");
        free(ctx); 
        return NULL; 
    }
    
    // Access angle mapping
    shmid = shmget(keyp + 3, sizeof(tAng) * (size_t)ctx->nrot * ctx->nrot, 0666);
    if (shmid < 0) { 
        fprintf(stderr, "Error: Failed to access angle mapping shared memory\n");
        free(ctx); 
        return NULL; 
    }
    ctx->angle_mapping = (tAng*)shmat(shmid, NULL, 0);
    if (ctx->angle_mapping == (void*)-1) { 
        fprintf(stderr, "Error: Failed to attach angle mapping shared memory\n");
        free(ctx); 
        return NULL; 
    }
    
    // Access rotation matrices
    shmid = shmget(keyp + 5, sizeof(double) * (size_t)ctx->nrot * 9, 0666);
    if (shmid < 0) { 
        fprintf(stderr, "Error: Failed to access rotation matrices shared memory\n");
        free(ctx); 
        return NULL; 
    }
    ctx->rotation_matrices = (double (*)[3][3])shmat(shmid, NULL, 0);
    if (ctx->rotation_matrices == (void*)-1) { 
        fprintf(stderr, "Error: Failed to attach rotation matrices shared memory\n");
        free(ctx); 
        return NULL; 
    }
   
    printf("FMAP initialized successfully with key string '%s'\n", key_string);
    printf("FMAP parameters: nrot=%d, blen=%d, nbin=%d\n",
           ctx->nrot, ctx->blen, ctx->nbin); 
    return ctx;
#else
    fprintf(stderr, "Error: Shared memory not supported on this platform\n");
    return NULL;
#endif
}

// Cleanup FMAP context
void fmap_cleanup(FMAP_Context* ctx) {
    if (!ctx) return;
    
#ifdef __unix__
    if (ctx->energy_data != (void*)-1) shmdt(ctx->energy_data);
    if (ctx->angle_indices != (void*)-1) shmdt(ctx->angle_indices);
    if (ctx->angle_mapping != (void*)-1) shmdt(ctx->angle_mapping);
    if (ctx->rotation_matrices != (void*)-1) shmdt(ctx->rotation_matrices);
#endif
    
    free(ctx);
}

// Calculate energy using direct bin indices
double fmap_energy_bins(const FMAP_Context* ctx,
                       const double pos1[3], const double bins1[3],
                       const double pos2[3], const double bins2[3],
                       const double box[3]) {
    if (!ctx) return 0.0;
    
    double rel_coords[3];
    if (!calculate_distance_and_coords(ctx, pos1, pos2, box, rel_coords)) {
        return 0.0; // Beyond cutoff
    }
    return calculate_energy_from_orientations(ctx, rel_coords, bins1, bins2);
}

// Calculate energy using Euler angles
double fmap_energy_euler(const FMAP_Context* ctx,
                        const double pos1[3], const double euler1[3],
                        const double pos2[3], const double euler2[3],
                        const double box[3]) {
    if (!ctx) return 0.0;
    
    double rel_coords[3];
    if (!calculate_distance_and_coords(ctx, pos1, pos2, box, rel_coords)) {
        return 0.0; // Beyond cutoff
    }
    
    // Convert orientations to bins only when needed
    double bins1[3], bins2[3];
    fmap_euler_to_bins(euler1, bins1, ctx->nbin);
    fmap_euler_to_bins(euler2, bins2, ctx->nbin);
    
    return calculate_energy_from_orientations(ctx, rel_coords, bins1, bins2);
}

// Calculate energy using quaternions
double fmap_energy_quaternion(const FMAP_Context* ctx,
                             const double pos1[3], const double quat1[4],
                             const double pos2[3], const double quat2[4],
                             const double box[3]) {
    if (!ctx) return 0.0;
    
    double rel_coords[3];
    if (!calculate_distance_and_coords(ctx, pos1, pos2, box, rel_coords)) {
        return 0.0; // Beyond cutoff
    }
    
    // Convert orientations to bins only when needed
    double bins1[3], bins2[3];
    fmap_quat_to_bins(quat1, bins1, ctx->nbin);
    fmap_quat_to_bins(quat2, bins2, ctx->nbin);
    
    return calculate_energy_from_orientations(ctx, rel_coords, bins1, bins2);
}
