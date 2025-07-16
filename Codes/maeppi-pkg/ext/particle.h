// particle.h - Generic particle interface with focused RNG interface
#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include "fmap_energy.h"

// Forward declarations
typedef struct Particle Particle;

// Random number generator function pointer type
typedef double (*RandomFunc)(void *rng_state);

// Function pointers for particle operations
typedef struct {
    // Memory management
    Particle* (*create)(void);
    void (*destroy)(Particle *p);
    
    // Coordinate access
    void (*get_position)(const Particle *p, double pos[3]);
    void (*set_position)(Particle *p, const double pos[3]);
    
    // File I/O
    int (*load_from_line)(Particle *p, const char *line);
    void (*save_to_file)(const Particle *p, FILE *file, int particle_id);
    
    // Monte Carlo moves - now uses focused RNG interface
    void (*generate_random)(Particle *p, const double box_len[3], RandomFunc rng_func, void *rng_state);
    void (*generate_trial_move)(const Particle *original, Particle *trial, 
                               const double max_trans_param, double max_rot_param,
                               const double box_len[3], double density_scale, 
                               RandomFunc rng_func, void *rng_state);
    
    // Energy calculation (direct FMAP interface)
    double (*pair_energy)(const Particle *p1, const Particle *p2, const double box_len[3], const FMAP_Context *fmap_ctx);
    
    // Utility
    void (*copy)(const Particle *src, Particle *dst);
    void (*print_info)(const Particle *p);
    
} ParticleOps;

// Global particle operations (set by particle implementation)
extern ParticleOps particle_ops;

// Convenience macros for cleaner code
#define PARTICLE_CREATE() particle_ops.create()
#define PARTICLE_DESTROY(p) particle_ops.destroy(p)
#define PARTICLE_GET_POS(p, pos) particle_ops.get_position(p, pos)
#define PARTICLE_SET_POS(p, pos) particle_ops.set_position(p, pos)
#define PARTICLE_LOAD(p, line) particle_ops.load_from_line(p, line)
#define PARTICLE_SAVE(p, file, id) particle_ops.save_to_file(p, file, id)
#define PARTICLE_RANDOM(p, box_len, rng_func, rng_state) particle_ops.generate_random(p, box_len, rng_func, rng_state)
#define PARTICLE_TRIAL_MOVE(orig, trial, trans, rot, box_len, scale, rng_func, rng_state) \
    particle_ops.generate_trial_move(orig, trial, trans, rot, box_len, scale, rng_func, rng_state)
#define PARTICLE_ENERGY(p1, p2, box_len, fmap_ctx) particle_ops.pair_energy(p1, p2, box_len, fmap_ctx)
#define PARTICLE_COPY(src, dst) particle_ops.copy(src, dst)
#define PARTICLE_PRINT(p) particle_ops.print_info(p)

#endif // PARTICLE_H
