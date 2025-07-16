import os
import sys
import glob
from cffi import FFI

# Set data path
os.environ["FMAPB2_DATA_PATH"] = os.path.join(os.path.dirname(__file__), "data")

# Create FFI instance
ffi = FFI()

# Define C structures and functions
ffi.cdef("""
    typedef struct {
        char info[31];
        double xyz[3];
        double q, r, Asq, Bsq;
    } ATOM;
    
    typedef struct {
        double Cut, VLow, kappa, sdie, kBT, escl, vscl;
    } SOLV;
    
    typedef signed char tErn;
    
    void softgrd(int nCrd, ATOM Crds[], int nPro, ATOM Pros[], 
                 int l, double dx, SOLV* sol, tErn softfb[]);
    void softgrd_cleanup(void);
""")

# Load the shared library
def _load_library():
    """Load the compiled shared library"""
    lib_dir = os.path.dirname(__file__)
    
    # Try different library extensions
    if sys.platform == "win32":
        lib_names = ["libfmapb2.dll", "fmapb2.dll"]
    elif sys.platform == "darwin":
        lib_names = ["libfmapb2.dylib"]
    else:
        lib_names = ["libfmapb2.so"]
    
    for lib_name in lib_names:
        lib_path = os.path.join(lib_dir, lib_name)
        if os.path.exists(lib_path):
            try:
                return ffi.dlopen(lib_path)
            except OSError as e:
                print(f"Warning: Failed to load {lib_path}: {e}")
                continue
    
    # If no local library found, try system library
    try:
        return ffi.dlopen(glob.glob(os.path.join(os.path.dirname(__file__), "libfmapb2*"))[0])
    except IndexError:
        raise ImportError(
        "Could not find libfmapb2 shared library. Try make."
        )

# Load library and expose interface
lib = _load_library()

# Export the same interface as before for compatibility
__all__ = ['ffi', 'lib']
