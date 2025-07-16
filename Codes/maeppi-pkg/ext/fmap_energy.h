#ifndef FMAP_ENERGY_H
#define FMAP_ENERGY_H

#ifdef __cplusplus
extern "C" {
#endif

// Forward declaration
typedef struct FMAP_Context FMAP_Context;

// ============================================================================
// CORE FMAP FUNCTIONS
// ============================================================================

// Create/destroy FMAP context
FMAP_Context* fmap_create(const char* key_string);
void fmap_cleanup(FMAP_Context* ctx);

// Energy calculation functions
double fmap_energy_bins(const FMAP_Context* ctx,
                       const double pos1[3], const double bins1[3],
                       const double pos2[3], const double bins2[3],
                       const double box[3]);

double fmap_energy_euler(const FMAP_Context* ctx,
                        const double pos1[3], const double euler1[3],
                        const double pos2[3], const double euler2[3],
                        const double box[3]);

double fmap_energy_quaternion(const FMAP_Context* ctx,
                             const double pos1[3], const double quat1[4],
                             const double pos2[3], const double quat2[4],
                             const double box[3]);

// ============================================================================
// COORDINATE CONVERSION FUNCTIONS
// ============================================================================

// Quaternion <-> Euler conversions
void fmap_quat_to_euler(const double quat[4], double euler[3]);
void fmap_euler_to_quat(const double euler[3], double quat[4]);

// Euler <-> Bins conversions
void fmap_euler_to_bins(const double euler[3], double bins[3], int nbin);
void fmap_bins_to_euler(const double bins[3], double euler[3], int nbin);

// Direct Quaternion <-> Bins conversions (optimized paths)
void fmap_quat_to_bins(const double quat[4], double bins[3], int nbin);
void fmap_bins_to_quat(const double bins[3], double quat[4], int nbin);

// Quaternion utilities
int fmap_quat_normalize(double quat[4]);

#ifdef __cplusplus
}
#endif

#endif // FMAP_ENERGY_H