// particle_simple.c - Simple 6-coordinate particle implementation with focused RNG interface
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "particle.h"
#include "fmap_energy.h"

static const double rot_default = 3.0; // rotation displacement in bin units

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Simple particle structure
struct Particle {
    double coord[6];  // [x, y, z, theta1, theta2, theta3]
};

// Apply periodic boundary conditions
static void apply_pbc(double coord[6], const double box_len[3]) {
    // Position coordinates
    for (int i = 0; i < 3; i++) {
        coord[i] = fmod(coord[i], box_len[i]);
        if (coord[i] < 0) coord[i] += box_len[i];
    }
    
    // Orientation coordinates (0-120 degrees)
    for (int i = 3; i < 6; i++) {
        coord[i] = fmod(coord[i], 120.0);
        if (coord[i] < 0) coord[i] += 120.0;
    }
}

// Implementation of particle operations
static Particle* create_particle(void) {
    Particle *p = malloc(sizeof(Particle));
    if (p) {
        memset(p->coord, 0, 6 * sizeof(double));
    }
    return p;
}

static void destroy_particle(Particle *p) {
    free(p);
}

static void get_position(const Particle *p, double pos[3]) {
    memcpy(pos, p->coord, 3 * sizeof(double));
}

static void set_position(Particle *p, const double pos[3]) {
    memcpy(p->coord, pos, 3 * sizeof(double));
}

static int load_from_line(Particle *p, const char *line) {
    double coord[6];
    double euler[3];
    int parsed = sscanf(line, "%*d %lf %lf %lf %lf %lf %lf",
                       &coord[0], &coord[1], &coord[2],
                       &euler[0], &euler[1], &euler[2]);
    
    if (parsed == 6) {
        memcpy(p->coord, coord, 3 * sizeof(double));
        fmap_euler_to_bins(euler, &(p->coord[3]), 120);  // Use fmap function
        return 1;
    }
    return 0;
}

static void save_to_file(const Particle *p, FILE *file, int particle_id) {
    double euler[3];
    fmap_bins_to_euler(&(p->coord[3]), euler, 120);  // Use fmap function
    fprintf(file, "%16d%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n",
            particle_id,
            p->coord[0], p->coord[1], p->coord[2],
            euler[0], euler[1], euler[2]);
}

static void generate_random(Particle *p, const double box_len[3], RandomFunc rng_func, void *rng_state) {
    // Generate random position and orientation
    for (int i = 0; i < 3; i++) {
        p->coord[i] = rng_func(rng_state) * box_len[i];
        p->coord[i+3] = rng_func(rng_state) * 120.0;
    }
}

static void generate_trial_move(const Particle *original, Particle *trial, 
                               const double max_trans_param, double max_rot_param,
                               const double box_len[3], double density_scale, 
                               RandomFunc rng_func, void *rng_state) {
    // Copy original coordinates
    memcpy(trial->coord, original->coord, 6 * sizeof(double));

    // Apply random displacement
    for (int i = 0; i < 3; i++) {
        trial->coord[i] += (rng_func(rng_state) - 0.5) * 
                          max_trans_param * density_scale;
        trial->coord[i+3] += (rng_func(rng_state) - 0.5) * 
                            max_rot_param * rot_default;
    }
    
    // Apply periodic boundary conditions
    apply_pbc(trial->coord, box_len);
}

static double pair_energy(const Particle *p1, const Particle *p2, const double box_len[3], const FMAP_Context *fmap_ctx) {
    return fmap_energy_bins(fmap_ctx, p1->coord, &p1->coord[3], p2->coord, &p2->coord[3], box_len);
}

static void copy_particle(const Particle *src, Particle *dst) {
    memcpy(dst, src, sizeof(Particle));
}

static void print_info(const Particle *p) {
    printf("Particle: coord=[%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n",
           p->coord[0], p->coord[1], p->coord[2],
           p->coord[3], p->coord[4], p->coord[5]);
}

// Initialize the particle operations structure
ParticleOps particle_ops = {
    .create = create_particle,
    .destroy = destroy_particle,
    .get_position = get_position,
    .set_position = set_position,
    .load_from_line = load_from_line,
    .save_to_file = save_to_file,
    .generate_random = generate_random,
    .generate_trial_move = generate_trial_move,
    .pair_energy = pair_energy,
    .copy = copy_particle,
    .print_info = print_info
};
