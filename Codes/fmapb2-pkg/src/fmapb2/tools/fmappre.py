#!/usr/bin/env python

import sys
import math
from operator import add
from math import sqrt

VDW_DATA = {
    'H': [1.78179743628, 0.0200],
    'O': [2.85087589805, 0.2000],
    'C': [3.56359487256, 0.1500],
    'N': [3.11814551349, 0.1600],
    'S': [3.56359487256, 0.2000],
    'P': [3.74177461619, 0.2],
    'IM': [4.40103966761, 0.1],
    'Li': [2.02590368505, 0.0183],
    'IP': [3.32839761097, 0.00277],
    'Na': [3.32839761097, 0.00277],
    'K': [4.73601758563, 0.000328],
    'Rb': [5.26699322165, 0.00017],
    'Cs': [6.04920229617, 8.06e-05],
    'MG': [1.412252648, 0.8947],
    'C0': [3.05239718809, 0.459789],
    'Zn': [1.95997717991, 0.0125],
    'F': [3.11814551349, 0.061],
    'Cl': [3.47094140587, 0.265],
    'Br': [3.95559030854, 0.32],
    'I': [4.18722397526, 0.4],
    'IB': [8.9089871814, 0.1],
}

ELEM_WEIGHTS = {
    "C": 12.01,
    "N": 14.01,
    "O": 16.00,
    "P": 30.97,
    "S": 32.06,
    "H": 1.01,
}

def get_kappa(ionic_strength, sdie, temp):
    """Calculate kappa parameter"""
    N_A = 6.022045000e+23
    e_c = 4.803242384e-10
    k_B = 1.380662000e-16
    Pi = 3.14159265358979323846
    return math.sqrt(ionic_strength * 1.0e-16 * (8.0 * Pi * N_A * e_c * e_c) / (1000.0 * sdie * k_B * temp))

def get_kbt(temp):
    """Calculate kBT parameter"""
    N_A = 6.022045000e+23
    k_B = 3.297623030e-24  # cal/K
    return N_A * k_B * temp / 1000.0

# https://nvlpubs.nist.gov/nistpubs/jres/56/jresv56n1p1_a1b.pdf
def get_water_die(temp):
    """Calculate water dielectric constant,
    Malmberg, Cyrus G. and Arthur A. Maryott. “Dielectric constant of water from 0 to 100 C.” Journal of research of the National Bureau of Standards 56 (1956): 1.
    """
    c2k = 273.15
    t = temp - c2k
    sdie = 87.740 - 0.40008 * t + 9.398e-4 * t * t - 1.410e-6 * t * t * t
    return sdie

def get_xyz(line):
    """Extract xyz coordinates from PDB line"""
    return [float(line[30:38]), float(line[38:46]), float(line[46:54])]

def get_center(fname):
    """Get center of coordinates"""
    cen = [0.0, 0.0, 0.0]
    count = 0.0
    with open(fname) as fp:
        for line in fp:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
                xyz = get_xyz(line)
                cen = list(map(add, cen, xyz))
                count += 1.0
    return cen[0]/count, cen[1]/count, cen[2]/count

def get_maxrad(fname):
    """Get maximum radius assuming file is centered"""
    rad = 0.0
    with open(fname) as fp:
        for line in fp:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
                xyz = get_xyz(line)
                cur = xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]
                if cur > rad:
                    rad = cur
    return sqrt(rad)

def get_mass(fname):
    """Calculate molecular mass"""
    mass = 0.0
    with open(fname) as fp:
        for line in fp:
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                if line[12] == 'H':
                    elem = 'H'
                else:
                    elem = line[13]
                if elem in ELEM_WEIGHTS:
                    mass += ELEM_WEIGHTS[elem]
    return mass

def get_charge(fname):
    """Get total charge from PQR file"""
    charge = 0.0
    with open(fname) as fp:
        for line in fp:
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                charge += float(line[54:62])
    return charge

def get_max_atom_radius(vdw_fname):
    """Get maximum atom radius from VDW file"""
    max_rad = 0.0
    with open(vdw_fname) as fp:
        for line in fp:
            if line[:4] == "ATOM":
                rad = float(line[62:70])
                if rad > max_rad:
                    max_rad = rad
    return max_rad

def calculate_vscl(mass):
    """Calculate vscl from mass using the formula"""
    vs = 0.18 * (80 + 0.27 * mass / 1000.0) / (80 + mass / 1000.0)
    return round(vs, 3)

def vdw_assign_center(pdb_fname, vdw_fname):
    """Assign VDW parameters directly from PQR file, Center"""
    coords = []
    lines = []
    
    # Read all atom lines and coordinates
    with open(pdb_fname) as pdb_file:
        for line in pdb_file:
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                lines.append(line)
                coords.append(get_xyz(line))
    
    # Calculate center
    center = [sum(coord[i] for coord in coords) / len(coords) for i in range(3)]
    
    # Write VDW file with centered coordinates
    with open(vdw_fname, 'w') as vdw_file:
        for line in lines:
            # Fix H* and determine atom type
            if line[12] == "H":
                info=line[:13]+'H'+line[14:30]
            else:
                info=line[:30]
            key = info[13]
            
            # Get VDW parameters
            sig, eps = VDW_DATA[key]
            sig6 = math.pow(sig, 6)
            sqA = math.sqrt(4.0 * eps * sig6 * sig6)
            sqB = math.sqrt(4.0 * eps * sig6)
            
            # Get original coordinates and center them
            x, y, z = get_xyz(line)
            x_cen = x - center[0]
            y_cen = y - center[1]
            z_cen = z - center[2]
            
            # Write line with centered coordinates and VDW parameters
            vdw_file.write(info + "%8.3f%8.3f%8.3f" % (x_cen, y_cen, z_cen) + 
                          line[54:62] + "%8.4f%8s%16.11f%16.11f%16.8e%16.8e\n" % 
                          (sig/2.0, key, sig, eps, sqA, sqB))


def main():
    if len(sys.argv) < 6:
        print("Usage: %s A.pqr ion temp escl vscl [B.pqr] >parms.txt" % sys.argv[0])
        print("Extra output: subA.vdw [subB.vdw]")
        print("vscl: use '-' to auto-calculate from mass")
        sys.exit(1)
    
    pqr_a = sys.argv[1]
    ion = float(sys.argv[2])
    temp_c = float(sys.argv[3])
    escl = 1.0 if (sys.argv[4] == "-") else float(sys.argv[4])
    vscl_input = sys.argv[5]
    
    # Check if we have molecule B
    has_mol_b = len(sys.argv) > 6
    pqr_b = sys.argv[6] if has_mol_b else None
    
    # Convert temperature to Kelvin
    temp_k = temp_c + 273.15
    
    # Process molecule A
    charge_a = get_charge(pqr_a)
    mass_a = get_mass(pqr_a)
    
    # Process molecule B if present
    if has_mol_b:
        charge_b = get_charge(pqr_b)
        mass_b = get_mass(pqr_b)
        mass = math.sqrt(mass_a * mass_b)
    else:
        mass = mass_a
    
    # Calculate or use provided vscl
    if vscl_input == '-':
        vscl = calculate_vscl(mass)
    else:
        vscl = float(vscl_input)
 
    # Calculate physical parameters
    sdie = get_water_die(temp_k)
    kappa = get_kappa(ion, sdie, temp_k)
    kbt = get_kbt(temp_k)
    qscl = 1.0 + 0.025 * pow(ion, -0.4)
    
    # Output parameters
    print("%-12s : %16.8E K" % ("Temperature", temp_k))
    print("%-12s : %16.8E" % ("kbt", kbt))
    print("%-12s : %16.8E M" % ("ion", ion))
    print("%-12s : %16.8E" % ("qscl", qscl))
    print("%-12s : %16.8E" % ("sdie", sdie))
    print("%-12s : %16.8E" % ("kappa", kappa))
    print("%-12s : %16.8E angstrom" % ("lambda", 1.0/kappa))
    print("%-12s : %16.3f" % ("rscl", 1.08))
    print("%-12s : %16.3f" % ("escl", escl))
    print("%-12s : %16.3f" % ("vscl", vscl))
    
    print("%-12s : %16.0f" % ("ChargeA", charge_a))
    print("%-12s : %16.3f Da" % ("MassA", mass_a))
    
    # Create VDW file for molecule A (centered internally)
    vdw_assign_center(pqr_a, "subA.vdw")
    max_rad_a = get_maxrad("subA.vdw")
    print("%-12s : %16.3f angstrom" % ("MaxRA", max_rad_a))
    
    # Process molecule B if present
    if has_mol_b:
        print("%-12s : %16.0f" % ("ChargeB", charge_b))
        print("%-12s : %16.3f Da" % ("MassB", mass_b))
        
        # Create VDW file for molecule B (centered internally)
        vdw_assign_center(pqr_b, "subB.vdw")
        max_rad_b = get_maxrad("subB.vdw")
        print("%-12s : %16.3f angstrom" % ("MaxRB", max_rad_b))
        
        # Calculate combined parameters
        max_r = (max_rad_a + max_rad_b) / 2.0
        
        # Get maximum atom radius from both files
        atom_rad_a = get_max_atom_radius("subA.vdw")
        atom_rad_b = get_max_atom_radius("subB.vdw")
        atom_r = max(atom_rad_a, atom_rad_b)
    else:
        max_r = max_rad_a
        atom_r = get_max_atom_radius("subA.vdw")
    
    print("%-12s : %16.3f Da" % ("mass", mass))
    print("%-12s : %16.3f angstrom" % ("MaxR", max_r))
    
    # Calculate rcut
    rcut = int(1.0/kappa * 3.0)
    if rcut < 36:
        rcut = 36
    print("%-12s : %16.3f angstrom" % ("rcut", rcut))
    
    # Calculate blen
    spacing = 0.6
    blen = 2 * (int((2 * (max_r + atom_r) + rcut) / spacing) + 1)
    print("%-12s : %16.1f" % ("spacing", spacing))
    print("%-12s : %16.0f" % ("blen", blen))

if __name__ == "__main__":
    main()
