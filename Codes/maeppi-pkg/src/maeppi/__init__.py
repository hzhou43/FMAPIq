"""
MAEPPI - Monte Carlo Energy Analysis for Protein-Protein Interactions

A Python package for FMAP energy calculations and Monte Carlo simulations.
"""

import os
import sys
import platform
import glob
from cffi import FFI

# Version info
__version__ = "0.1.0"
__author__ = "Your Name"

# Determine library extension based on platform
if platform.system() == "Darwin":
    LIB_EXT = "dylib"
elif platform.system() == "Windows":
    LIB_EXT = "dll"
else:
    LIB_EXT = "so"

# Get package directory
_package_dir = os.path.dirname(os.path.abspath(__file__))

# Find library files using glob patterns
def _find_library(pattern):
    """Find library file using glob pattern"""
    lib_pattern = os.path.join(_package_dir, pattern)
    matches = glob.glob(lib_pattern)
    if not matches:
        raise ImportError(f"Library not found matching pattern: {pattern}")
    return matches[0]

# Library paths - find the actual filenames
_fmap_lib_path = _find_library(f"libmaeppiern*.{LIB_EXT}")
_sim_lib_path = _find_library(f"libmaeppisim*.{LIB_EXT}")

# Initialize CFFI for FMAP energy library
fmap_ffi = FFI()
fmap_ffi.cdef("""
    // FMAP Context
    typedef struct FMAP_Context FMAP_Context;
    
    // Core FMAP functions
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
    
    // Coordinate conversion functions
    void fmap_quat_to_euler(const double quat[4], double euler[3]);
    void fmap_euler_to_quat(const double euler[3], double quat[4]);
    void fmap_euler_to_bins(const double euler[3], double bins[3], int nbin);
    void fmap_bins_to_euler(const double bins[3], double euler[3], int nbin);
    void fmap_quat_to_bins(const double quat[4], double bins[3], int nbin);
    void fmap_bins_to_quat(const double bins[3], double quat[4], int nbin);
    int fmap_quat_normalize(double quat[4]);
""")

# Initialize CFFI for simulation library
sim_ffi = FFI()
sim_ffi.cdef("""
    // Simulation Box
    typedef struct Box Box;
    
    // Simulation setup and management
    Box* sim_create_box(int n_particles, double box_length, double temperature, 
                        double max_trans_param, double max_rot_param, 
                        unsigned int seed, const char *key_string, double zratio);
    void sim_destroy_box(Box *box);
    int sim_initialize_random(Box *box);
    int sim_load_configuration(Box *box, const char *filename);
    int sim_save_configuration(const Box *box, const char *filename);
    
    // Simulation execution
    void sim_run_steps(Box *box, long n_steps);
    
    // Data access and trajectory
    void sim_get_stats(const Box *box, double *total_energy, long *n_moves, 
                       long *n_accepted, double *acceptance_rate);
    int sim_get_particle_coords(const Box *box, int particle_id, double coords[6]);
    int sim_set_particle_coords(Box *box, int particle_id, const double coords[6]);
    void sim_get_box_dimensions(const Box *box, double box_len[3]);
    void sim_get_parameters(const Box *box, int *n_particles, double *temperature,
                           double *max_trans_param, double *max_rot_param);
""")

# Load the libraries
try:
    ern = fmap_ffi.dlopen(_fmap_lib_path)
    sim = sim_ffi.dlopen(_sim_lib_path)
except Exception as e:
    raise ImportError(f"Failed to load MAEPPI libraries: {e}")

# Helper functions to hide FFI complexity
def _to_c_array(data, size):
    """Convert Python list/array to C array"""
    return fmap_ffi.new(f"double[{size}]", data)

# Export the main interfaces
__all__ = [
    'ern', 'sim', '__version__'
]

# Convenience wrapper classes
class Ern:
    """Wrapper for FMAP energy context with automatic cleanup."""
    
    def __init__(self, key_string):
        self._ctx = ern.fmap_create(key_string.encode('utf-8'))
        if self._ctx == fmap_ffi.NULL:
            raise RuntimeError("Failed to create FMAP context")
    
    def __del__(self):
        if hasattr(self, '_ctx') and self._ctx != fmap_ffi.NULL:
            ern.fmap_cleanup(self._ctx)
    
    def energy(self, pos1, orient1, pos2, orient2, box, kind='euler'):
        """Calculate energy using different orientation representations
        
        Args:
            pos1: Position of particle 1 [x, y, z]
            orient1: Orientation of particle 1 (format depends on kind)
            pos2: Position of particle 2 [x, y, z]
            orient2: Orientation of particle 2 (format depends on kind)
            box: Box dimensions [lx, ly, lz]
            kind: Type of orientation representation ('euler', 'quat', 'bins')
                - 'euler': Euler angles [α, β, γ]
                - 'quat': Quaternions [w, x, y, z]  
                - 'bins': Binned representation [bin1, bin2, bin3]
        
        Returns:
            Energy value as float
        """
        pos1_c = _to_c_array(pos1, 3)
        pos2_c = _to_c_array(pos2, 3)
        box_c = _to_c_array(box, 3)
        
        if kind == 'euler':
            orient1_c = _to_c_array(orient1, 3)
            orient2_c = _to_c_array(orient2, 3)
            return ern.fmap_energy_euler(self._ctx, pos1_c, orient1_c, pos2_c, orient2_c, box_c)
        elif kind == 'quat':
            orient1_c = _to_c_array(orient1, 4)
            orient2_c = _to_c_array(orient2, 4)
            return ern.fmap_energy_quaternion(self._ctx, pos1_c, orient1_c, pos2_c, orient2_c, box_c)
        elif kind == 'bins':
            orient1_c = _to_c_array(orient1, 3)
            orient2_c = _to_c_array(orient2, 3)
            return ern.fmap_energy_bins(self._ctx, pos1_c, orient1_c, pos2_c, orient2_c, box_c)
        else:
            raise ValueError(f"Unknown orientation kind: {kind}. Must be 'euler', 'quat', or 'bins'")
    

class Sim:
    """Wrapper for simulation box with automatic cleanup."""
    
    def __init__(self, n_particles, box_length, temperature, 
                 max_trans_param, max_rot_param, seed, key_string, zratio):
        self._box = sim.sim_create_box(
            n_particles, box_length, temperature,
            max_trans_param, max_rot_param, seed,
            key_string.encode('utf-8'), zratio
        )
        if self._box == sim_ffi.NULL:
            raise RuntimeError("Failed to create simulation box")
    
    def __del__(self):
        if hasattr(self, '_box') and self._box != sim_ffi.NULL:
            sim.sim_destroy_box(self._box)
    
    def run_steps(self, n_steps):
        """Run simulation steps"""
        sim.sim_run_steps(self._box, n_steps)
    
    def get_stats(self):
        """Get simulation statistics"""
        energy = sim_ffi.new("double *")
        moves = sim_ffi.new("long *")
        accepted = sim_ffi.new("long *")
        acceptance = sim_ffi.new("double *")
        sim.sim_get_stats(self._box, energy, moves, accepted, acceptance)
        return energy[0], moves[0], accepted[0], acceptance[0]

# Add convenience wrapper classes to exports
__all__.extend(['Ern', 'Sim'])

#print(f"MAEPPI v{__version__} initialized successfully")
#print(f"✓ FMAP library loaded: {_fmap_lib_path}")
#print(f"✓ Simulation library loaded: {_sim_lib_path}")
