"""
sim2iq - Simple SAXS calculator
"""

import os
import glob
import cffi
import numpy as np

_package_dir = os.path.dirname(__file__)

def find_library(pattern):
    """Find library file using glob pattern"""
    lib_pattern = os.path.join(_package_dir, pattern)
    matches = glob.glob(lib_pattern)
    if not matches:
        raise ImportError(f"Library not found matching pattern: {pattern}")
    return matches[0]

# Create CFFI instance
ffi = cffi.FFI()

# Define the C interface
ffi.cdef("""
    typedef struct {
        char element[8];
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

    Structure read_pqr_file(const char* filename);
    void structure_free(Structure* structure);
    void qidata_free(QIData* data);
    
    QIData calculate_saxs_data_multi_step(
        Structure* structure,
        int grid_size,
        double grid_spacing,
        int steps,
        double cutoff_offset,
        double se_weight,
        double ou_weight,
        double in_weight,
        double b_factor
    );
    
    Structure process_trajectory(Structure* original_structure, 
                               const char* traj_file, 
                               bool is_txt);
    
    double u2B(double u);
""")

class Calc:
    """SAXS Calculator"""
    
    def __init__(self):
        lib_path = find_library("libsim2iq.*")
        self.lib = ffi.dlopen(lib_path)
    
    def __call__(self, pqr_file, grid_size=64, grid_spacing=1.0, 
                 steps=1, cutoff_offset=2.8, se_weight=-1.0, 
                 ou_weight=0.011, in_weight=-0.011, b_factor=0.0,
                 traj_file=None, traj_is_txt=False):
        """Calculate SAXS profile"""
        
        structure = self.lib.read_pqr_file(pqr_file.encode('utf-8'))
        
        try:
            # Process trajectory if provided
            if traj_file:
                transformed_structure = self.lib.process_trajectory(
                    ffi.addressof(structure), 
                    traj_file.encode('utf-8'), 
                    traj_is_txt
                )
                structure_to_use = transformed_structure
            else:
                structure_to_use = structure
            
            if b_factor == 0.0:
                b_factor = self.lib.u2B(0.25 * grid_spacing)
            
            result = self.lib.calculate_saxs_data_multi_step(
                ffi.addressof(structure_to_use), grid_size, grid_spacing, 
                steps, cutoff_offset, se_weight, ou_weight, 
                in_weight, b_factor
            )
            
            q = np.array([result.data[i].q for i in range(result.count)])
            intensity = np.array([result.data[i].intensity for i in range(result.count)])
            
            self.lib.qidata_free(ffi.addressof(result))
            return q, intensity
            
        finally:
            if traj_file:
                self.lib.structure_free(ffi.addressof(transformed_structure))
            self.lib.structure_free(ffi.addressof(structure))
