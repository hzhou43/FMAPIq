// simulation_cffi.h - CFFI interface for Monte Carlo simulation
#ifndef SIMULATION_CFFI_H
#define SIMULATION_CFFI_H

#ifdef __cplusplus
extern "C" {
#endif

// Forward declaration of simulation box
typedef struct Box Box;

// ============================================================================
// SIMULATION SETUP AND MANAGEMENT
// ============================================================================

// Create and initialize simulation box
Box* sim_create_box(int n_particles, double box_length, double temperature, 
                    double max_trans_param, double max_rot_param, 
                    unsigned int seed, const char *key_string, double zratio);

// Destroy simulation box
void sim_destroy_box(Box *box);

// Initialize with random configuration
int sim_initialize_random(Box *box);

// Load configuration from file
int sim_load_configuration(Box *box, const char *filename);

// Save configuration to file
int sim_save_configuration(const Box *box, const char *filename);

// ============================================================================
// SIMULATION EXECUTION
// ============================================================================

// Run simulation for specified number of steps
void sim_run_steps(Box *box, long n_steps);

// ============================================================================
// DATA ACCESS AND TRAJECTORY
// ============================================================================

// Get current simulation statistics
void sim_get_stats(const Box *box, double *total_energy, long *n_moves, 
                   long *n_accepted, double *acceptance_rate);

// Get particle coordinates (positions and orientations as Euler angles)
// coords[0:3] = position [x, y, z]
// coords[3:6] = orientation [euler_x, euler_y, euler_z]
int sim_get_particle_coords(const Box *box, int particle_id, double coords[6]);

// Set particle coordinates (positions and orientations as Euler angles)
// coords[0:3] = position [x, y, z]
// coords[3:6] = orientation [euler_x, euler_y, euler_z]
int sim_set_particle_coords(Box *box, int particle_id, const double coords[6]);

// Get box dimensions
void sim_get_box_dimensions(const Box *box, double box_len[3]);

// Get simulation parameters
void sim_get_parameters(const Box *box, int *n_particles, double *temperature,
                       double *max_trans_param, double *max_rot_param);

#ifdef __cplusplus
}
#endif

#endif // SIMULATION_CFFI_H
