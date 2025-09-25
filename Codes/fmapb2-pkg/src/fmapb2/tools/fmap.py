#!/usr/bin/env python3
"""
Python CFFI replacement for fmap.cen.c with parameter file input
"""

import sys
import os
import zlib
import math
import traceback
from copy import deepcopy

try:
    from fmapb2 import ffi, lib
except ModuleNotFoundError:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    from fmapb2 import ffi, lib

class Atom:
    """Python representation of ATOM structure"""
    def __init__(self, info="", xyz=None, q=0.0, r=0.0, Asq=0.0, Bsq=0.0):
        self.info = info[:30]
        self.xyz = xyz if xyz is not None else [0.0, 0.0, 0.0]
        self.q = q
        self.r = r
        self.Asq = Asq
        self.Bsq = Bsq

class Solvent:
    """Python representation of SOLV structure"""
    def __init__(self, Cut=12.0, VLow=0.0, kappa=0.1, sdie=80.0, escl=1.0, vscl=0.16, kBT=0.591):
        self.Cut = Cut
        self.VLow = VLow
        self.kappa = kappa
        self.sdie = sdie
        self.escl = escl
        self.vscl = vscl
        self.kBT = kBT
    
    def to_c_struct(self):
        """Convert to CFFI C structure"""
        c_solv = ffi.new("SOLV *")
        c_solv.Cut = self.Cut
        c_solv.VLow = self.VLow
        c_solv.kappa = self.kappa
        c_solv.sdie = self.sdie
        c_solv.escl = self.escl
        c_solv.vscl = self.vscl
        c_solv.kBT = self.kBT
        return c_solv

def read_parameters(filename):
    """Read parameters from file"""
    params = {}
    # Only read the parameters we actually need
    needed_params = ['kbt', 'qscl', 'sdie', 'kappa', 'rscl', 'escl', 'vscl', 'rcut', 'spacing', 'blen']
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if ':' in line and not line.startswith('#'):
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Only process needed parameters
                    if key in needed_params:
                        # Extract numeric value (remove units)
                        value = value.split()[0]
                        try:
                            if key == 'blen':
                                params[key] = int(value)
                            else:
                                params[key] = float(value)
                        except ValueError:
                            continue
    except FileNotFoundError:
        print(f"Error: Parameter file {filename} not found")
        sys.exit(1)
    return params

def read_vdw(filename):
    """Read PQR file and return list of Atom objects"""
    atoms = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    info = line[:30]
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        q = float(line[54:62])
                        r = float(line[62:70])
                        
                        start_pos = 70 + 8 + 16 + 16
                        Asq = float(line[start_pos:start_pos+16])
                        
                        Bsq = float(line[start_pos+16:start_pos+32])
                        
                        atoms.append(Atom(info=info, xyz=[x, y, z], q=q, r=r, Asq=Asq, Bsq=Bsq))
                    
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse line: {line.strip()}")
                        continue
                        
    except FileNotFoundError:
        print(f"Error: File {filename} not found")
        sys.exit(1)
    return atoms

def read_angles(filename):
    """Read angles file"""
    angles = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    angles.append([float(parts[0]), float(parts[1]), float(parts[2])])
    except FileNotFoundError:
        print(f"Error: File {filename} not found")
        sys.exit(1)
    return angles

def euler_to_rotmat(psi, theta, phi):
    """Convert Euler angles to rotation matrix"""
    cos_psi, sin_psi = math.cos(psi), math.sin(psi)
    cos_theta, sin_theta = math.cos(theta), math.sin(theta)
    cos_phi, sin_phi = math.cos(phi), math.sin(phi)
    
    r11 = cos_psi * cos_phi - sin_psi * cos_theta * sin_phi
    r12 = -cos_psi * sin_phi - sin_psi * cos_theta * cos_phi
    r13 = sin_psi * sin_theta
    
    r21 = sin_psi * cos_phi + cos_psi * cos_theta * sin_phi
    r22 = -sin_psi * sin_phi + cos_psi * cos_theta * cos_phi
    r23 = -cos_psi * sin_theta
    
    r31 = sin_theta * sin_phi
    r32 = sin_theta * cos_phi
    r33 = cos_theta
    
    return [[r11, r12, r13],
            [r21, r22, r23],
            [r31, r32, r33]]

def rotate_atoms(atoms, rotmat):
    """Apply rotation matrix to atoms"""
    rotated_atoms = []
    for atom in atoms:
        new_x = rotmat[0][0] * atom.xyz[0] + rotmat[0][1] * atom.xyz[1] + rotmat[0][2] * atom.xyz[2]
        new_y = rotmat[1][0] * atom.xyz[0] + rotmat[1][1] * atom.xyz[1] + rotmat[1][2] * atom.xyz[2]
        new_z = rotmat[2][0] * atom.xyz[0] + rotmat[2][1] * atom.xyz[1] + rotmat[2][2] * atom.xyz[2]
        
        new_atom = Atom(info=atom.info, xyz=[new_x, new_y, new_z], 
                       q=atom.q, r=atom.r, Asq=atom.Asq, Bsq=atom.Bsq)
        rotated_atoms.append(new_atom)
    return rotated_atoms

def scale_coordinates(atoms, dx):
    """Scale coordinates by dx"""
    scaled_atoms = []
    for atom in atoms:
        new_xyz = [atom.xyz[i] / dx for i in range(3)]
        new_atom = Atom(info=atom.info, xyz=new_xyz, 
                       q=atom.q, r=atom.r, Asq=atom.Asq, Bsq=atom.Bsq)
        scaled_atoms.append(new_atom)
    return scaled_atoms

def atoms_to_c_array(atoms):
    """Convert Python atoms list to C array"""
    c_atoms = ffi.new("ATOM[]", len(atoms))
    for i, atom in enumerate(atoms):
        c_atoms[i].info = atom.info.encode('utf-8')
        for j in range(3):
            c_atoms[i].xyz[j] = atom.xyz[j]
        c_atoms[i].q = atom.q
        c_atoms[i].r = atom.r
        c_atoms[i].Asq = atom.Asq
        c_atoms[i].Bsq = atom.Bsq
    return c_atoms

def save_grid_data(filename, softfb_array, l):
    """Save grid data to file with gzip compression"""
    buffer_data = ffi.buffer(softfb_array, l * l * l * ffi.sizeof("tErn"))
    raw_data = bytes(buffer_data)
    
    # Use zlib compressobj to compress in gzip format
    compressor = zlib.compressobj(wbits=16+15)
    compressed_data = compressor.compress(raw_data)
    compressed_data += compressor.flush()
    
    with open(filename, 'wb') as f:
        f.write(compressed_data)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]  # Get command line args, excluding script name
    
    if len(argv) < 3:
        print("Usage: python %s subA.vdw parms.txt ang.dat [fnamepat] [subB.vdw]" % (sys.argv[0] if argv is None else "fmapb2"), file=sys.stderr)
        sys.exit(1)
    
    # Parse command line arguments
    subA_file = argv[0]
    parms_file = argv[1]
    ang_file = argv[2]
    fnamepat = "ve/tErn.{:05d}.mat.gz" if len(argv) < 4 else argv[3] 
    subB_file = argv[4] if len(argv) > 4 else None
    
    # Validate format pattern
    try:
        fnamepat.format(0)
    except Exception as e:
        print(f"Invalid filename pattern: {fnamepat}",file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
    
    # Read parameters
    params = read_parameters(parms_file)
    #print(params) 
    
    # Extract needed parameters
    try:
        kbt = params['kbt']
        kappa = params['kappa']
        sdie = params['sdie']
        rscl = params['rscl']
        qscl = params['qscl']
        escl = params['escl']
        vscl = params['vscl']
        rcut = params['rcut']
        dx = params['spacing']  # Use spacing from parameter file
        blen = int(params['blen'])
    except KeyError as e:
        print(f"Missing required parameter: {e}",file=sys.stderr)
        sys.exit(1)
    print(f"Grid: {blen}x{blen}x{blen}, spacing: {dx} A, rcut: {rcut:.0f} A",file=sys.stderr)
    
    # Create solvent object
    solvent = Solvent(
        Cut=rcut,
        VLow=0.0,
        kappa=kappa,
        sdie=sdie,
        escl=escl,
        vscl=vscl,
        kBT=kbt,
    )
    c_solvent = solvent.to_c_struct()
    
    # Read coordinate files
    crds = read_vdw(subA_file)  # A is always crds
    
    if subB_file:
        pros = read_vdw(subB_file)  # B is pros when provided
    else:
        pros = deepcopy(crds)  
    
    # Scale radii and charges
    for atom in crds:
        atom.r *= rscl
        atom.q *= qscl
    for atom in pros:
        atom.r *= rscl
    
    # Scale coordinates by dx
    crds = scale_coordinates(crds, dx)
    pros = scale_coordinates(pros, dx)
    
    # Read angles
    angles = read_angles(ang_file)
    
    print(f"Loaded: {len(crds)} crds, {len(pros)} pros, {len(angles)} rotations",file=sys.stderr)
    
    # Create output directory if it doesn't exist
    path = os.path.dirname(fnamepat)
    if path:
        os.makedirs(path, exist_ok=True)
    
    # Convert to C arrays
    c_crds = atoms_to_c_array(crds)
    
    # Process each orientation
    print("Processing orientations...(check output file for progress)", flush=True,file=sys.stderr)
    for i, angle in enumerate(angles):
        
        # Rotate probe atoms
        rotmat = euler_to_rotmat(angle[0], angle[1], angle[2])
        rotated_pros = rotate_atoms(pros, rotmat)
        c_pros = atoms_to_c_array(rotated_pros)
        
        # Allocate softfb array
        softfb_array = ffi.new("tErn[]", blen * blen * blen)
        
        # Call the main computation function
        lib.softgrd(len(crds), c_crds, len(rotated_pros), c_pros, blen, dx, c_solvent, softfb_array)
        
        # Save the result
        filename = fnamepat.format(i)
        save_grid_data(filename, softfb_array, blen)
    
    print()  # New line after progress
    
    lib.softgrd_cleanup()
    print("Completed successfully!",file=sys.stderr)

if __name__ == "__main__":
    main()
