// simulation.c - Generic Monte Carlo simulation engine with CFFI interface
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "particle.h"
#include "fmap_energy.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Linear Congruential Generator state
typedef struct {
    unsigned long state;
} LCG_State;

// Simulation box structure
typedef struct {
    int n_particles;
    double box_len[3];
    double temperature;
    
    Particle **particles;         // Array of particle pointers
    double *particle_energies;   // Individual particle energies
    double **energy_matrix;      // Pair energy matrix [i][j] for i < j
    double total_energy;
    
    // FMAP contexts array (for future multi-type support)
    FMAP_Context **fmap_contexts;
    int n_fmap_contexts;
    
    // Move parameters
    double max_translation_param;       
    double max_rotation_param;
    
    // Statistics
    long n_moves;
    long n_accepted;
    
    // Random number generator state and function
    void *rng_state;              // Generic pointer to RNG state
    RandomFunc random_func;       // Function pointer to random generator
    
    // Temporary storage for trial move energies
    double *trial_energies;
    
    // Reusable trial particle - allocated once, reused for all moves
    Particle *trial_particle;
    
    // Density scale for move scaling
    double density_scale;
} Box;

// Simple linear congruential generator
static double lcg_random(void *state) {
    LCG_State *lcg = (LCG_State*)state;
    lcg->state = (lcg->state * 1103515245UL + 12345UL) & 0x7fffffffUL;
    return (double)(lcg->state) / (double)0x7fffffffUL;
}

// Allocate energy matrix (upper triangular)
static double** allocate_energy_matrix(int n) {
    double **matrix = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        matrix[i] = malloc(n * sizeof(double));
        memset(matrix[i], 0, n * sizeof(double));
    }
    return matrix;
}

// Free energy matrix
static void free_energy_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Initialize FMAP contexts
static int initialize_fmap_contexts(Box *box, const char *key_string) {
    box->n_fmap_contexts = 1;
    box->fmap_contexts = malloc(box->n_fmap_contexts * sizeof(FMAP_Context*));
    
    if (!box->fmap_contexts) {
        printf("Error: Failed to allocate FMAP contexts array\n");
        return 0;
    }
    
    box->fmap_contexts[0] = fmap_create(key_string);
    if (!box->fmap_contexts[0]) {
        printf("Error: Failed to create FMAP context with key string '%s'\n", key_string);
        free(box->fmap_contexts);
        box->fmap_contexts = NULL;
        return 0;
    }
    
    printf("FMAP context 0 initialized successfully with key string '%s'\n", key_string);
    return 1;
}

// Cleanup FMAP contexts
static void cleanup_fmap_contexts(Box *box) {
    if (box->fmap_contexts) {
        for (int i = 0; i < box->n_fmap_contexts; i++) {
            if (box->fmap_contexts[i]) {
                fmap_cleanup(box->fmap_contexts[i]);
                box->fmap_contexts[i] = NULL;
            }
        }
        free(box->fmap_contexts);
        box->fmap_contexts = NULL;
    }
    box->n_fmap_contexts = 0;
}

// Calculate energies for loaded configuration
static void calculate_loaded_energies(Box *box) {
    printf("Calculating energies for loaded configuration...\n");
    
    // Initialize energy arrays
    memset(box->particle_energies, 0, box->n_particles * sizeof(double));
    
    // Calculate all pairwise energies and build energy matrix
    #pragma omp parallel for schedule(guided,16)
    for (int i = 0; i < box->n_particles; i++) {
        for (int j = i + 1; j < box->n_particles; j++) {
            double energy = PARTICLE_ENERGY(box->particles[i], box->particles[j], box->box_len, box->fmap_contexts[0]);
            box->energy_matrix[i][j] = energy;
            
            #pragma omp atomic
            box->particle_energies[i] += energy;
            #pragma omp atomic
            box->particle_energies[j] += energy;
        }
    }
    
    printf("Energy calculation completed for loaded configuration\n");
}

// Initialize particles with random positions and orientations
static int initialize_particles_random(Box *box) {
    printf("Initializing %d particles randomly...\n", box->n_particles);
    
    for (int i = 0; i < box->n_particles; i++) {
        int max_attempts = 1000;
        int attempts = 0;
        double insertion_energy;
        
        while (attempts < max_attempts) {
            PARTICLE_RANDOM(box->particles[i], box->box_len, box->random_func, box->rng_state);
            
            insertion_energy = 0.0;
            for (int j = 0; j < i; j++) {
                double energy = PARTICLE_ENERGY(box->particles[j], box->particles[i], box->box_len, box->fmap_contexts[0]);
                insertion_energy += energy;
                box->energy_matrix[j][i] = energy;
            }
            
            if (insertion_energy < 999.0) {
                box->particle_energies[i] = insertion_energy;
                
                for (int j = 0; j < i; j++) {
                    box->particle_energies[j] += box->energy_matrix[j][i];
                }
                break;
            }
            
            attempts++;
        }
        
        if (attempts >= max_attempts) {
            printf("Failed to place particle %d after %d attempts\n", i, max_attempts);
            return 0;
        }
        
        if ((i + 1) % (box->n_particles / 2) == 0 || i == box->n_particles - 1) {
            printf("Placed %d/%d particles\n", i + 1, box->n_particles);
        }
    }
    
    return 1;
}

// Calculate initial total energy
static void calculate_initial_energy(Box *box) {
    box->total_energy = 0.0;
    
    for (int i = 0; i < box->n_particles; i++) {
        for (int j = i + 1; j < box->n_particles; j++) {
            box->total_energy += box->energy_matrix[i][j];
        }
    }
}

// Calculate energy change for a trial move
static double calculate_energy_change(Box *box, int particle_id, const Particle *trial) {
    double old_energy = box->particle_energies[particle_id];
    double new_energy = 0.0;
    
    #pragma omp parallel for reduction(+:new_energy) schedule(guided, 16)
    for (int j = 0; j < box->n_particles; j++) {
        if (j == particle_id) {
            box->trial_energies[j] = 0.0;
            continue;
        }
        
        double energy;
        if (j < particle_id) {
            energy = PARTICLE_ENERGY(box->particles[j], trial, box->box_len, box->fmap_contexts[0]);
        } else {
            energy = PARTICLE_ENERGY(trial, box->particles[j], box->box_len, box->fmap_contexts[0]);
        }
        
        box->trial_energies[j] = energy;
        new_energy += energy;
    }
    
    return new_energy - old_energy;
}

// Accept move and update data structures
static void accept_move(Box *box, int particle_id, const Particle *trial, double delta_energy) {
    PARTICLE_COPY(trial, box->particles[particle_id]);
    
    double new_particle_energy = 0.0;
    
    for (int j = 0; j < box->n_particles; j++) {
        if (j == particle_id) continue;
        
        double new_energy = box->trial_energies[j];
        double old_energy;
        
        if (j < particle_id) {
            old_energy = box->energy_matrix[j][particle_id];
            box->energy_matrix[j][particle_id] = new_energy;
        } else {
            old_energy = box->energy_matrix[particle_id][j];
            box->energy_matrix[particle_id][j] = new_energy;
        }
        
        box->particle_energies[j] += (new_energy - old_energy);
        new_particle_energy += new_energy;
    }
    
    box->particle_energies[particle_id] = new_particle_energy;
    box->total_energy += delta_energy;
}

// Perform one Monte Carlo step
static void monte_carlo_step(Box *box) {
    int particle_id = (int)(ceil(box->random_func(box->rng_state) * box->n_particles))-1;
    if (particle_id==-1){particle_id=0;}
    
    PARTICLE_TRIAL_MOVE(box->particles[particle_id], box->trial_particle, 
                        box->max_translation_param, box->max_rotation_param,
                        box->box_len, box->density_scale, box->random_func, box->rng_state);
    
    double delta_energy = calculate_energy_change(box, particle_id, box->trial_particle);
    
    int accept = 0;
    if (delta_energy <= 0) {
        accept = 1;
    } else {
        double probability = exp(-delta_energy / box->temperature);
        if (box->random_func(box->rng_state) < probability) {
            accept = 1;
        }
    }
    
    if (accept) {
        accept_move(box, particle_id, box->trial_particle, delta_energy);
        box->n_accepted++;
    }
    
    box->n_moves++;
}

// Cleanup box
static void destroy_box(Box *box) {
    if (box) {
        for (int i = 0; i < box->n_particles; i++) {
            PARTICLE_DESTROY(box->particles[i]);
        }
        free(box->particles);
        if (box->particle_energies) free(box->particle_energies);
        if (box->energy_matrix) free_energy_matrix(box->energy_matrix, box->n_particles);
        if (box->trial_energies) free(box->trial_energies);
        if (box->trial_particle) PARTICLE_DESTROY(box->trial_particle);
        cleanup_fmap_contexts(box);
        if (box->rng_state) free(box->rng_state);
        free(box);
    }
}

// ============================================================================
// CFFI INTERFACE FUNCTIONS
// ============================================================================

// Create and initialize simulation box
Box* sim_create_box(int n_particles, double box_length, double temperature, 
                    double max_trans_param, double max_rot_param, 
                    unsigned int seed, const char *key_string, double zratio) {
    Box *box = malloc(sizeof(Box));
    if (!box) return NULL;
    
    box->n_particles = n_particles;
    box->box_len[0] = box->box_len[1] = box_length;
    box->box_len[2] = box_length * zratio;
    box->temperature = temperature;
    
    // Initialize random number generator
    LCG_State *lcg_state = malloc(sizeof(LCG_State));
    if (!lcg_state) {
        free(box);
        return NULL;
    }
    lcg_state->state = seed;
    box->rng_state = lcg_state;
    box->random_func = lcg_random;
    
    // Initialize FMAP contexts
    box->fmap_contexts = NULL;
    box->n_fmap_contexts = 0;
    
    // Set move parameters
    const double max_trans_default = 1.62;
    box->max_translation_param = max_trans_param * max_trans_default;
    box->max_rotation_param = max_rot_param;
    
    box->total_energy = 0.0;
    box->n_moves = 0;
    box->n_accepted = 0;
    
    // Calculate density scale
    const double sigma = 40.5;
    double volume = box->box_len[0] * box->box_len[1] * box->box_len[2];
    double density = box->n_particles * (sigma * sigma * sigma) / volume;
    box->density_scale = 0.3 / density;
    
    if (!initialize_fmap_contexts(box, key_string)) {
        free(lcg_state);
        free(box);
        return NULL;
    }
    
    // Allocate particle array
    box->particles = malloc(n_particles * sizeof(Particle*));
    for (int i = 0; i < n_particles; i++) {
        box->particles[i] = PARTICLE_CREATE();
        if (!box->particles[i]) {
            for (int j = 0; j < i; j++) {
                PARTICLE_DESTROY(box->particles[j]);
            }
            free(box->particles);
            cleanup_fmap_contexts(box);
            free(lcg_state);
            free(box);
            return NULL;
        }
    }
    
    // Allocate reusable trial particle
    box->trial_particle = PARTICLE_CREATE();
    if (!box->trial_particle) {
        printf("Error: Failed to create trial particle\n");
        destroy_box(box);
        return NULL;
    }
    
    box->particle_energies = malloc(n_particles * sizeof(double));
    box->energy_matrix = allocate_energy_matrix(n_particles);
    box->trial_energies = malloc(n_particles * sizeof(double));
    
    if (!box->particle_energies || !box->energy_matrix || !box->trial_energies) {
        printf("Error: Memory allocation failed\n");
        destroy_box(box);
        return NULL;
    }
    
    memset(box->particle_energies, 0, n_particles * sizeof(double));
    
    return box;
}

// Destroy simulation box
void sim_destroy_box(Box *box) {
    destroy_box(box);
}

// Initialize with random configuration
int sim_initialize_random(Box *box) {
    if (!box) return 0;
    
    if (!initialize_particles_random(box)) {
        return 0;
    }
    calculate_initial_energy(box);
    return 1;
}

// Load configuration from file
int sim_load_configuration(Box *box, const char *filename) {
    if (!box) return 0;
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error: Cannot open %s for reading\n", filename);
        return 0;
    }
    
    char line[256];
    int particle_count = 0;
    
    printf("Loading configuration from %s...\n", filename);
    
    while (fgets(line, sizeof(line), file) && particle_count < box->n_particles) {
        if (line[0] == '#') continue;
        
        if (PARTICLE_LOAD(box->particles[particle_count], line)) {
            particle_count++;
        }
    }
    
    fclose(file);
    
    if (particle_count != box->n_particles) {
        printf("Error: Expected %d particles, but loaded %d from %s\n", 
               box->n_particles, particle_count, filename);
        return 0;
    }
    
    printf("Successfully loaded %d particles from %s\n", particle_count, filename);
    calculate_loaded_energies(box);
    calculate_initial_energy(box);
    return 1;
}

// Save configuration to file
int sim_save_configuration(const Box *box, const char *filename) {
    if (!box) return 0;
    
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Cannot open %s for writing\n", filename);
        return 0;
    }
    
    fprintf(file, "# Box state with %d particles\n", box->n_particles);
    fprintf(file, "# Total energy: %.6f\n", box->total_energy);
    
    for (int i = 0; i < box->n_particles; i++) {
        PARTICLE_SAVE(box->particles[i], file, i + 1);
    }
    
    fclose(file);
    printf("Configuration saved to %s\n", filename);
    return 1;
}

// Run simulation for specified number of steps
void sim_run_steps(Box *box, long n_steps) {
    if (!box) return;
    
    for (long step = 0; step < n_steps; step++) {
        monte_carlo_step(box);
    }
}

// Get current simulation statistics
void sim_get_stats(const Box *box, double *total_energy, long *n_moves, 
                   long *n_accepted, double *acceptance_rate) {
    if (!box) return;
    
    if (total_energy) *total_energy = box->total_energy;
    if (n_moves) *n_moves = box->n_moves;
    if (n_accepted) *n_accepted = box->n_accepted;
    if (acceptance_rate) {
        *acceptance_rate = (box->n_moves > 0) ? (double)box->n_accepted / box->n_moves : 0.0;
    }
}

// Get particle coordinates (positions and orientations as Euler angles)
int sim_get_particle_coords(const Box *box, int particle_id, double coords[6]) {
    if (!box || particle_id < 0 || particle_id >= box->n_particles) return 0;
    
    // Get position
    PARTICLE_GET_POS(box->particles[particle_id], coords);
    
    // For simple particles, we need to convert bins to Euler angles
    // This assumes we're using the simple particle implementation
    double euler[3];
    double bins[3];
    
    // Get orientation bins from particle (assuming simple particle structure)
    // This is a bit of a hack - in a real implementation, we'd have a proper interface
    const Particle *p = box->particles[particle_id];
    // Access the coordinate array directly (this assumes simple particle structure)
    const double *particle_coords = (const double*)p + 3;  // Skip position, get orientation
    bins[0] = particle_coords[0];
    bins[1] = particle_coords[1]; 
    bins[2] = particle_coords[2];
    
    fmap_bins_to_euler(bins, euler, 120);
    coords[3] = euler[0];
    coords[4] = euler[1];
    coords[5] = euler[2];
    
    return 1;
}

// Set particle coordinates (positions and orientations as Euler angles)
int sim_set_particle_coords(Box *box, int particle_id, const double coords[6]) {
    if (!box || particle_id < 0 || particle_id >= box->n_particles) return 0;
    
    // Set position
    PARTICLE_SET_POS(box->particles[particle_id], coords);
    
    // Convert Euler angles to bins and set orientation
    double euler[3] = {coords[3], coords[4], coords[5]};
    double bins[3];
    fmap_euler_to_bins(euler, bins, 120);
    
    // Set orientation bins (this assumes simple particle structure)
    Particle *p = box->particles[particle_id];
    double *particle_coords = (double*)p + 3;  // Skip position, get orientation
    particle_coords[0] = bins[0];
    particle_coords[1] = bins[1];
    particle_coords[2] = bins[2];
    
    return 1;
}

// Get box dimensions
void sim_get_box_dimensions(const Box *box, double box_len[3]) {
    if (!box || !box_len) return;
    
    box_len[0] = box->box_len[0];
    box_len[1] = box->box_len[1];
    box_len[2] = box->box_len[2];
}

// Get simulation parameters
void sim_get_parameters(const Box *box, int *n_particles, double *temperature,
                       double *max_trans_param, double *max_rot_param) {
    if (!box) return;
    
    if (n_particles) *n_particles = box->n_particles;
    if (temperature) *temperature = box->temperature;
    if (max_trans_param) *max_trans_param = box->max_translation_param;
    if (max_rot_param) *max_rot_param = box->max_rotation_param;
}

// ============================================================================
// STANDALONE MAIN (for backward compatibility)
// ============================================================================

int main(int argc, char *argv[]) {
    if (argc < 10) {
        printf("Usage: %s <n_particles> <box_length> <temperature> <n_steps> <seed> <start.dat> <stop.dat> <sample_freq> <keyp> [tx] [rx] [zx]\n", argv[0]);
        printf("  start.dat: '-' for random initialization, or filename to read initial configuration\n");
        printf("  stop.dat:  filename to save final configuration\n");
        printf("  keyp:      shared memory key for FMAP data (hex string, e.g., 'abc123')\n");
        return 1;
    }
    
    int n_particles = atoi(argv[1]);
    double box_length = atof(argv[2]);
    double temperature = atof(argv[3]);
    long n_steps = atol(argv[4]);
    unsigned int seed = (unsigned int)atoi(argv[5]);
    char *start_file = argv[6];
    char *stop_file = argv[7];
    long sample_freq = atol(argv[8]);
    char *key_string = argv[9];
    
    double max_trans_param = (argc > 10) ? atof(argv[10]) : 1.0;
    double max_rot_param = (argc > 11) ? atof(argv[11]) : 1.0;
    double zratio = (argc > 12) ? atof(argv[12]) : 1.0;
    
    printf("Starting NVT Monte Carlo simulation\n");
    
    Box *box = sim_create_box(n_particles, box_length, temperature, 
                             max_trans_param, max_rot_param, seed, key_string, zratio);
    if (!box) {
        printf("Error: Failed to create simulation box\n");
        return 1;
    }
    
    if (strcmp(start_file, "-") == 0) {
        printf("Using random initialization for configuration\n");
        if (!sim_initialize_random(box)) {
            printf("Error: Failed to initialize particles\n");
            sim_destroy_box(box);
            return 1;
        }
    } else {
        if (!sim_load_configuration(box, start_file)) {
            printf("Error: Failed to load configuration from %s\n", start_file);
            sim_destroy_box(box);
            return 1;
        }
    }
    
    printf("\nRunning %ld MC steps...\n", n_steps);
    printf("Temperature: %.3f\n", box->temperature);
    printf("Number of particles: %d\n", box->n_particles);
    printf("Box dimensions: %.3f x %.3f x %.3f\n", 
           box->box_len[0], box->box_len[1], box->box_len[2]);
    
    double initial_energy;
    sim_get_stats(box, &initial_energy, NULL, NULL, NULL);
    printf("Initial energy: %.6f\n", initial_energy);
    
    #ifdef _OPENMP
    printf("OpenMP threads: %d\n", omp_get_max_threads());
    #endif
    printf("----------------------------------------\n");
    
    clock_t start_time = clock();
    
    for (long step = 0; step < n_steps; step++) {
        sim_run_steps(box, 1);
        
        if ((step + 1) % sample_freq == 0) {
            double energy;
            long moves, accepted;
            double acceptance_rate;
            sim_get_stats(box, &energy, &moves, &accepted, &acceptance_rate);
            
            double elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
            double steps_per_sec = (step + 1) / elapsed;
            
            printf("Step %8ld: Energy = %12.6f, Accept = %.3f, Speed = %.0f steps/s\n",
                   step + 1, energy, acceptance_rate, steps_per_sec);
        }
    }
    
    double final_energy;
    long final_moves, final_accepted;
    double final_acceptance;
    sim_get_stats(box, &final_energy, &final_moves, &final_accepted, &final_acceptance);
    
    double total_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    
    printf("----------------------------------------\n");
    printf("Simulation completed!\n");
    printf("Final energy: %.6f\n", final_energy);
    printf("Overall acceptance rate: %.3f\n", final_acceptance);
    printf("Total time: %.2f seconds\n", total_time);
    printf("Average speed: %.0f steps/second\n", n_steps / total_time);
    
    sim_save_configuration(box, stop_file);
    sim_destroy_box(box);
    
    return 0;
}
