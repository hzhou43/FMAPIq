#!/usr/bin/env python3
"""
Python replacement for shm_manager.c using CFFI
Maintains full compatibility with C clients using System V IPC
"""

import os
import sys
import time
import signal
from cffi import FFI
from pathlib import Path

# Initialize CFFI
ffi = FFI()

# Define C structures and functions we need
ffi.cdef("""
    // System V IPC
    int shmget(int key, size_t size, int shmflg);
    void* shmat(int shmid, const void *shmaddr, int shmflg);
    int shmdt(const void *shmaddr);
    int shmctl(int shmid, int cmd, struct shmid_ds *buf);
    
    // File operations for gzip
    typedef struct gzFile_s *gzFile;
    gzFile gzopen(const char *path, const char *mode);
    int gzread(gzFile file, void *buf, unsigned int len);
    int gzwrite(gzFile file, const void *buf, unsigned int len);
    int gzclose(gzFile file);
    
    // Standard C
    void* malloc(size_t size);
    void free(void *ptr);
    FILE* fopen(const char *filename, const char *mode);
    int fclose(FILE *stream);
    size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
    size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
    
    // Constants
    #define IPC_CREAT 01000
    #define IPC_RMID 0
    #define NULL ...
""")

# Load system libraries
try:
    libc = ffi.dlopen(None)  # Standard C library
    libz = ffi.dlopen("libz.so.1")  # zlib for gzip
except OSError as e:
    print(f"Error loading libraries: {e}")
    sys.exit(1)

class SimpleINIParser:
    """
    Simple INI parser that supports:
    - [section] headers
    - key = value pairs
    - Comments starting with ';'
    - Whitespace trimming
    """
    
    def __init__(self):
        self.data = {}
    
    def parse(self, filename):
        """Parse INI file and return nested dictionary"""
        self.data = {}
        current_section = None
        
        try:
            with open(filename, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Skip comments
                    if line.startswith(';'):
                        continue
                    
                    # Remove inline comments
                    if ';' in line:
                        line = line.split(';')[0].strip()
                    
                    # Section header
                    if line.startswith('[') and line.endswith(']'):
                        current_section = line[1:-1].strip()
                        if current_section not in self.data:
                            self.data[current_section] = {}
                        continue
                    
                    # Key-value pair
                    if '=' in line and current_section is not None:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        self.data[current_section][key] = value
                        continue
                    
                    # Invalid line
                    if line:  # Only warn for non-empty lines
                        print(f"Warning: Invalid line {line_num} in {filename}: {line}")
        
        except FileNotFoundError:
            print(f"Error: Config file not found: {filename}")
            sys.exit(1)
        except Exception as e:
            print(f"Error parsing {filename}: {e}")
            sys.exit(1)
        
        return self.data
    
    def get(self, section, key, default=None):
        """Get value from section.key, return default if not found"""
        return self.data.get(section, {}).get(key, default)

class Configuration:
    """Configuration container matching the C structure"""
    
    def __init__(self):
        self.version = 1
        self.keyp = 0
        self.nrot = 0
        self.blen = 0
        self.nbin = 0
        self.fpat_ern = None
        self.fpat_rot = None
        self.spacing = 0.0
        self.scld = 0.0
        self.lsel = 0
    
    def load_from_ini(self, filename):
        """Load configuration from INI or TXT file"""
        if filename.endswith('.txt'):
            self.load_from_txt(filename)
        else:
            # Original INI loading code
            parser = SimpleINIParser()
            data = parser.parse(filename)
            
            # Parse protocol section
            if 'protocol' in data:
                self.version = int(data['protocol'].get('version', 1))
            
            # Parse FMAP section
            if 'FMAP' in data:
                fmap = data['FMAP']
                
                # Parse keyp as hex value with << 4 shift
                if 'keyp' in fmap:
                    keyp_hex = int(fmap['keyp'], 16)
                    self.keyp = keyp_hex << 4
                
                # Parse integer values
                self.nrot = int(fmap.get('nrot', 0))
                self.blen = int(fmap.get('blen', 0))
                self.nbin = int(fmap.get('nbin', 0))
                self.lsel = int(fmap.get('lsel', 0))
                
                # Parse float values
                self.spacing = float(fmap.get('spacing', 0.0))
                self.scld = float(fmap.get('scld', 0.0))
                
                # Parse string values
                self.fpat_ern = fmap.get('fpat_ern')
                self.fpat_rot = fmap.get('fpat_rot')
            
            # Validate required fields
            if self.keyp == 0:
                print("Error: keyp not specified in config")
                sys.exit(1)
            if self.nrot <= 0 or self.blen <= 0 or self.nbin <= 0:
                print("Error: nrot, blen, nbin must be positive")
                sys.exit(1)
    
    def load_from_txt(self, filename):
        """Load configuration from TXT file like parms.txt"""
        import hashlib
        
        # Generate keyp from absolute path using simple hash
        abs_path = os.path.abspath(filename)
        hash_obj = hashlib.md5(abs_path.encode())
        keyp_hex = int(hash_obj.hexdigest()[:6], 16)
        self.keyp = keyp_hex << 4
        
        print(f"Generated keyp: 0x{keyp_hex:06x}")
        
        # Parse txt file for blen and nrot
        try:
            with open(filename, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('blen'):
                        parts = line.split(':')
                        if len(parts) >= 2:
                            self.blen = int(parts[1].strip())
                    elif line.startswith('nrot'):
                        parts = line.split(':')
                        if len(parts) >= 2:
                            self.nrot = int(parts[1].strip())
                    elif line.startswith('spacing'):
                        parts = line.split(':')
                        if len(parts) >= 2:
                            self.spacing = float(parts[1].strip())
        except FileNotFoundError:
            print(f"Error: Config file not found: {filename}")
            sys.exit(1)
        except Exception as e:
            print(f"Error parsing {filename}: {e}")
            sys.exit(1)
        
        # Set defaults for missing values
        if self.nrot == 0:
            self.nrot = 72
        self.nbin = 120
        self.lsel = self.blen
        self.scld = 1.0
        if self.spacing == 0.0:
            self.spacing = 0.6
        
        # Set file patterns based on script path
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.fpat_ern = "ve/tErn.%05d.mat.gz"
        
        # Try script-relative path first
        fpat_rot = os.path.join(script_dir, "../data/sspNF%05d" % self.nrot)
        if os.path.exists(fpat_rot):
            self.fpat_rot = fpat_rot
        else:
            # Try maeppi package path
            try:
                import maeppi
                maeppi_path = maeppi.__path__[0]
                fpat_rot = os.path.join(maeppi_path, "data/sspNF%05d" % self.nrot)
                if os.path.exists(fpat_rot):
                    self.fpat_rot = fpat_rot
                else:
                    print(f"Error: Could not find data directory sspNF{self.nrot:05d}")
                    print(f"Tried: {os.path.join(script_dir, '../data/sspNF%05d' % self.nrot)}")
                    print(f"Tried: {fpat_rot}")
                    sys.exit(1)
            except ImportError:
                print(f"Error: Could not find data directory sspNF{self.nrot:05d}")
                print(f"Tried: {os.path.join(script_dir, '../data/sspNF%05d' % self.nrot)}")
                print("maeppi package not available")
                sys.exit(1)
        
        self.lsel = self.blen
        
        # Validate
        if self.nrot <= 0 or self.blen <= 0:
            print("Error: nrot and blen must be positive")
            sys.exit(1)
        
        print(f"Loaded from txt: blen={self.blen}, nrot={self.nrot}, spacing={self.spacing}")

class SharedMemoryManager:
    """Manages System V shared memory segments using CFFI"""
    
    def __init__(self):
        self.segments = {}  # Track created segments
        self.killed = False
        
    def create_segment(self, key, size, init_data=None):
        """Create or attach to shared memory segment"""
        # Create shared memory segment
        shmid = libc.shmget(key, size, 0o666 | 0o1000)  # IPC_CREAT | 0666
        if shmid < 0:
            print(f"Error: shmget failed for key 0x{key:x}")
            sys.exit(1)
        
        # Attach to segment
        ptr = libc.shmat(shmid, ffi.NULL, 0)
        if ptr == ffi.cast("void*", -1):
            print(f"Error: shmat failed for key 0x{key:x}")
            sys.exit(1)
        
        # Store segment info
        self.segments[key] = {
            'shmid': shmid,
            'ptr': ptr,
            'size': size
        }
        
        # Initialize data if provided
        if init_data is not None:
            ffi.memmove(ptr, init_data, len(init_data))
        
        return ptr
    
    def attach_segment(self, key, size):
        """Attach to existing shared memory segment"""
        shmid = libc.shmget(key, size, 0o666)
        if shmid < 0:
            return None  # Segment doesn't exist
        
        ptr = libc.shmat(shmid, ffi.NULL, 0)
        if ptr == ffi.cast("void*", -1):
            return None
        
        self.segments[key] = {
            'shmid': shmid,
            'ptr': ptr,
            'size': size
        }
        
        return ptr
    
    def remove_segment(self, key):
        """Remove shared memory segment"""
        if key not in self.segments:
            # Try to get segment ID for cleanup
            shmid = libc.shmget(key, 0, 0o666)
            if shmid >= 0:
                libc.shmctl(shmid, 0, ffi.NULL)  # IPC_RMID
            return
        
        segment = self.segments[key]
        libc.shmdt(segment['ptr'])
        libc.shmctl(segment['shmid'], 0, ffi.NULL)  # IPC_RMID
        del self.segments[key]
    
    def cleanup_all(self):
        """Clean up all segments"""
        for key in list(self.segments.keys()):
            self.remove_segment(key)

def load_file_data(filename, element_size, count):
    """Load binary data from file (with gzip support)"""
    if filename.endswith('.gz'):
        # Handle gzip files
        gz_file = libz.gzopen(filename.encode('utf-8'), b"rb")
        if gz_file == ffi.NULL:
            print(f"Error: Cannot open gzip file {filename}")
            return None
        
        buffer = ffi.new("char[]", count * element_size)
        bytes_read = libz.gzread(gz_file, buffer, count * element_size)
        libz.gzclose(gz_file)
        
        if bytes_read != count * element_size:
            print(f"Error: Expected {count * element_size} bytes, read {bytes_read}")
            return None
        
        return buffer
    else:
        # Handle regular files
        file_ptr = libc.fopen(filename.encode('utf-8'), b"rb")
        if file_ptr == ffi.NULL:
            print(f"Error: Cannot open file {filename}")
            return None
        
        buffer = ffi.new("char[]", count * element_size)
        items_read = libc.fread(buffer, element_size, count, file_ptr)
        libc.fclose(file_ptr)
        
        if items_read != count:
            print(f"Error: Expected {count} items, read {items_read}")
            return None
        
        return buffer

def run_server(config_file, load_only=False):
    """Run the shared memory server"""
    mode_str = "load-only" if load_only else "persistent"
    print(f"Starting server ({mode_str}) with config: {config_file}")
    
    # Load configuration
    config = Configuration()
    config.load_from_ini(config_file)
    
    # Create shared memory manager
    shm_manager = SharedMemoryManager()
    
    # Set up signal handler for graceful shutdown
    def signal_handler(signum, frame):
        print(f"\nReceived signal {signum}, shutting down...")
        shm_manager.killed = True
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    try:
        # Create parameter storage segment (10 doubles)
        param_size = 10 * 8  # 10 doubles * 8 bytes each
        param_ptr = shm_manager.create_segment(config.keyp, param_size)
        
        # Initialize parameters as double array
        params = ffi.cast("double*", param_ptr)
        params[0] = float(config.nrot)   # nrot
        params[1] = float(config.blen)   # blen  
        params[2] = float(config.nbin)   # nbin
        params[3] = config.scld          # scld
        params[4] = config.spacing       # spacing
        params[5] = float(config.lsel)   # lsel
        params[6] = float(config.keyp)   # keyp
        params[7] = 0.0                  # reserved
        params[8] = 0.0                  # reserved
        params[9] = 0.0                  # reserved
        
        # Create data segments with correct type sizes
        # typedef signed char tErn;  -> 1 byte
        # typedef short tAng;        -> 2 bytes  
        # typedef double;            -> 8 bytes
        segments_info = [
            (config.keyp + 1, "ern", config.blen**3 * config.nrot, 1),  # tErn = signed char = 1 byte
            (config.keyp + 2, "ijk", config.nbin**3, 2),                # tAng = short = 2 bytes  
            (config.keyp + 3, "ij2k", config.nrot**2, 2),               # tAng = short = 2 bytes
            (config.keyp + 4, "rot", config.nrot * 3, 8),               # double = 8 bytes
            (config.keyp + 5, "rotm", config.nrot * 9, 8),              # double = 8 bytes
        ]
        
        for key, name, count, elem_size in segments_info:
            size = count * elem_size
            ptr = shm_manager.create_segment(key, size)
            print(f"Created {name} segment: key=0x{key:x}, size={size}")
            
            # Load data if files are specified
            if name == "ern" and config.fpat_ern:
                # Load multiple files for ern data
                for i in range(config.nrot):
                    filename = config.fpat_ern % i
                    if os.path.exists(filename):
                        file_data = load_file_data(filename, elem_size, config.blen**3)
                        if file_data:
                            offset = i * config.blen**3 * elem_size
                            ffi.memmove(ffi.cast("char*", ptr) + offset, file_data, config.blen**3 * elem_size)
                            print(f"Loaded {filename}")
            
            elif config.fpat_rot:
                # Load single files for other data types
                file_map = {
                    "ijk": "angijk.dat.gz",
                    "ij2k": "angij2k.dat.gz", 
                    "rot": "rot.dat.gz",
                    "rotm": "rotm.dat.gz"
                }
                
                if name in file_map:
                    filename = f"{config.fpat_rot}/{file_map[name]}"
                    if os.path.exists(filename):
                        file_data = load_file_data(filename, elem_size, count)
                        if file_data:
                            ffi.memmove(ptr, file_data, size)
                            print(f"Loaded {filename}")
        
        print(f"Server running: key=0x{config.keyp:x}, nrot={config.nrot}, blen={config.blen}, nbin={config.nbin}")
        print(f"Additional params: scld={config.scld:.6f}, spacing={config.spacing:.6f}, lsel={config.lsel}")
        
        if load_only:
            print("Data loaded successfully. Server exiting.")
            print(f"Shared memory segments created with base key: 0x{config.keyp:x}")
            print("Use 'ipcrm -M <key>' to remove segments manually if needed")
            print("Or use 'python shm_manager.py cleanup config.ini' for automatic cleanup")
            return 0
        
        print("Press Ctrl+C to stop")
        
        # Main server loop - wait for client to set params[0] to 0 or negative
        while params[0] > 0 and not shm_manager.killed:
            time.sleep(1)
        
        print("Server shutdown requested")
        
    except Exception as e:
        print(f"Server error: {e}")
    finally:
        # Cleanup only in persistent mode or on error
        if not load_only or shm_manager.killed:
            print("Cleaning up shared memory segments...")
            shm_manager.cleanup_all()
            print("Server stopped")

def stop_server(config_file):
    """Stop the running server"""
    config = Configuration()
    config.load_from_ini(config_file)
    
    shm_manager = SharedMemoryManager()
    
    # Try to attach to parameter segment
    param_ptr = shm_manager.attach_segment(config.keyp, 80)
    if param_ptr is None:
        print(f"No server running with key 0x{config.keyp:x}")
        return 1
    
    # Send stop signal by setting params[0] to -1
    params = ffi.cast("double*", param_ptr)
    params[0] = -1.0
    print("Stop signal sent to server")
    return 0

def check_status(config_file):
    """Check status of shared memory segments"""
    config = Configuration()
    config.load_from_ini(config_file)
    
    shm_manager = SharedMemoryManager()
    
    # Try to attach to parameter segment
    param_ptr = shm_manager.attach_segment(config.keyp, 80)  # 10 doubles
    if param_ptr is None:
        print(f"Status: No server running with key 0x{config.keyp:x}")
        return 1
    
    # Read parameters
    params = ffi.cast("double*", param_ptr)
    print(f"Status: key=0x{config.keyp:x} nrot={params[0]:.0f} blen={params[1]:.0f} nbin={params[2]:.0f}")
    print(f"Additional params: scld={params[3]:.6f} spacing={params[4]:.6f} lsel={params[5]:.0f} keyp=0x{int(params[6]):x}")
    
    # Check other segments with correct sizes
    segments_info = [
        (config.keyp + 1, "shm_ern", config.blen**3 * config.nrot * 1),  # tErn = 1 byte
        (config.keyp + 2, "shm_ijk", config.nbin**3 * 2),                # tAng = 2 bytes
        (config.keyp + 3, "shm_ij2k", config.nrot**2 * 2),               # tAng = 2 bytes
        (config.keyp + 4, "shm_rot", config.nrot * 3 * 8),               # double = 8 bytes
        (config.keyp + 5, "shm_rotm", config.nrot * 9 * 8),              # double = 8 bytes
    ]
    
    for key, name, size in segments_info:
        ptr = shm_manager.attach_segment(key, size)
        print(f"{name}: {ptr if ptr else 'NULL'}")
    
    return 0

def cleanup_segments(config_file):
    """Manually cleanup shared memory segments"""
    config = Configuration()
    config.load_from_ini(config_file)
    
    shm_manager = SharedMemoryManager()
    
    # List of all segment keys
    keys = [
        config.keyp,        # parameters
        config.keyp + 1,    # ern
        config.keyp + 2,    # ijk  
        config.keyp + 3,    # ij2k
        config.keyp + 4,    # rot
        config.keyp + 5,    # rotm
    ]
    
    print(f"Cleaning up shared memory segments for key base: 0x{config.keyp:x}")
    
    cleaned = 0
    for key in keys:
        try:
            shm_manager.remove_segment(key)
            print(f"Removed segment: 0x{key:x}")
            cleaned += 1
        except:
            # Segment probably doesn't exist, which is fine
            pass
    
    print(f"Cleanup complete. Removed {cleaned} segments.")
    return 0

def usage(prog):
    """Print usage information"""
    print(f"Usage: {prog} <command> <config.ini|config.txt> [options]")
    print("Commands:")
    print("  start       - Start the server")
    print("    --load-only : Load data and exit (don't keep server running)")
    print("  status      - Check status of shared memory segments") 
    print("  stop        - Stop the running server")
    print("  cleanup     - Remove all shared memory segments")
    print()
    print("Config files:")
    print("  .ini files  - Use INI format with [FMAP] section")
    print("  .txt files  - Use parameter format (blen:, nrot:, spacing:)")
    print()
    print("Examples:")
    print(f"  {prog} start parms.txt            # Persistent server")
    print(f"  {prog} start parms.txt --load-only # Load data from txt and exit")
    print(f"  {prog} status parms.txt           # Check what's running")
    print(f"  {prog} stop parms.txt             # Stop persistent server")
    print(f"  {prog} cleanup parms.txt          # Force cleanup segments")

def main():
    if len(sys.argv) < 3:
        usage(sys.argv[0])
        sys.exit(1)
    
    command = sys.argv[1]
    config_file = sys.argv[2]
    
    if command == "start":
        # Check for --load-only flag
        load_only = "--load-only" in sys.argv
        return run_server(config_file, load_only)
    elif command == "status":
        return check_status(config_file)
    elif command == "stop":
        return stop_server(config_file)
    elif command == "cleanup":
        return cleanup_segments(config_file)
    else:
        print(f"Unknown command: {command}")
        usage(sys.argv[0])
        return 1

if __name__ == "__main__":
    sys.exit(main())
