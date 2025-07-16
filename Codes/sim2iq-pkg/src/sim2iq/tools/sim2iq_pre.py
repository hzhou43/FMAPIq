#!/usr/bin/env python
#
# sim2iq_pre.py
# A simplified tool to read PDB files, load unique atomic volumes
# from a JSON file, and output a PQR file with scaled radii
#

import os
import sys
import numpy as np
import json
import argparse

# van der Waals radii (in Å)
vdW = {
    'H': 1.2, 'HE': 1.4, 'LI': 1.82, 'BE': 1.53, 'B': 1.92, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'NE': 1.54,
    'NA': 2.27, 'MG': 1.73, 'AL': 1.84, 'SI': 2.1, 'P': 1.8, 'S': 1.8, 'CL': 1.75, 'AR': 1.88, 'K': 2.75, 'CA': 2.31
}


class PDB(object):
    """Load and process PDB file."""
    
    def __init__(self, filename):
        self.filename = filename
        self.natoms = 0
        self.read_pdb(filename)
        self.unique_radius = None
        self.unique_volume = None
    
    def read_pdb(self, filename):
        # Store raw PDB lines for volume assignment and PQR output
        self.pdb_lines = []
        self.atomname = []
        self.resname = []
        self.vdW = []
        
        with open(filename) as f:
            for line in f:
                if line[0:6] == "ENDMDL":
                    break
                if line[0:4] != "ATOM" and line[0:4] != "HETA":
                    continue
                
                # Keep first 54 characters as input for PQR output
                self.pdb_lines.append(line[:54])
                
                # Extract only info needed for radius assignment
                atomname = line[12:16].strip()
                resname = line[17:20]
                
                self.atomname.append(atomname)
                self.resname.append(resname)
                
                # Determine atom type from atom name for VdW fallback
                if len(atomname) > 0 and atomname[0].isalpha():
                    atomtype = atomname[0]
                    if len(atomname) > 1 and atomname[1].isalpha() and atomname[1].islower():
                        atomtype += atomname[1]
                else:
                    atomtype = "C"
                
                atomtype = atomtype.upper()
                
                # Set VDW radius as fallback
                try:
                    dr = vdW[atomtype]
                except:
                    try:
                        dr = vdW[atomtype[0]]
                    except:
                        dr = vdW['C']
                
                self.vdW.append(dr)
        
        self.natoms = len(self.pdb_lines)
    
    def load_unique_volumes(self, volume_file):
        """Load unique volumes from a JSON file"""
        print(f"Loading unique volumes from {volume_file}")
        
        # Residue variant mapping dictionary
        residue_mapping = {
            # Histidine variants
            "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSE": "HIS", "HSP": "HIS", "HSC": "HIS", "HSD": "HIS",
            # Cysteine variants
            "CYX": "CYS", "CYM": "CYS", "CSO": "CYS", "CSS": "CYS",
            # Aspartic acid variants
            "ASH": "ASP",
            # Glutamic acid variants
            "GLH": "GLU",
            # Lysine variants
            "LYN": "LYS",
            # Arginine variants
            "ARN": "ARG",
            # Terminal variants
            "NALA": "ALA", "NGLY": "GLY", "NSER": "SER", "NTHR": "THR", "NLEU": "LEU", "NILE": "ILE",
            "NVAL": "VAL", "NASN": "ASN", "NGLN": "GLN", "NARG": "ARG", "NHID": "HIS", "NHIE": "HIS",
            "NPRO": "PRO", "NPHE": "PHE", "NTRP": "TRP", "NTYR": "TYR", "NGLU": "GLU", "NASP": "ASP",
            "NLYS": "LYS", "NMET": "MET", "NCYS": "CYS",
            "CALA": "ALA", "CGLY": "GLY", "CSER": "SER", "CTHR": "THR", "CLEU": "LEU", "CILE": "ILE",
            "CVAL": "VAL", "CASN": "ASN", "CGLN": "GLN", "CARG": "ARG", "CHID": "HIS", "CHIE": "HIS",
            "CPRO": "PRO", "CPHE": "PHE", "CTRP": "TRP", "CTYR": "TYR", "CGLU": "GLU", "CASP": "ASP",
            "CLYS": "LYS", "CMET": "MET", "CCYS": "CYS",
        }
        
        try:
            with open(volume_file, 'r') as f:
                volume_data = json.load(f)
            
            # Initialize unique_radius list
            self.unique_radius = []
            
            # Assign volumes based on residue and atom names
            for i in range(self.natoms):
                resname = self.resname[i]
                atomname = self.atomname[i]
                
                # Try exact match first
                if resname in volume_data and atomname in volume_data[resname]:
                    volume = volume_data[resname][atomname]
                    radius = sphere_radius_from_volume(volume)
                
                # Try using residue mapping for non-standard residues
                elif resname in residue_mapping:
                    mapped_resname = residue_mapping[resname]
                    if mapped_resname in volume_data and atomname in volume_data[mapped_resname]:
                        volume = volume_data[mapped_resname][atomname]
                        radius = sphere_radius_from_volume(volume)
                    else:
                        # Use VdW radius as fallback
                        radius = self.vdW[i]
                
                else:
                    # Use VdW radius as fallback
                    radius = self.vdW[i]
                
                self.unique_radius.append(radius)
            
            return True
            
        except Exception as e:
            print(f"Error loading volume data: {e}")
            return False


def sphere_volume_from_radius(R):
    """Calculate volume of a sphere given its radius"""
    V_sphere = 4 * np.pi / 3 * R**3
    return V_sphere


def sphere_radius_from_volume(V):
    """Calculate radius of a sphere given its volume"""
    R_sphere = (3 * V / (4 * np.pi)) ** (1. / 3)
    return R_sphere


def write_pqr(pdb, output_filename):
    """Write a PQR file with the calculated unique radii"""
    print(f"Writing PQR file to {output_filename}")
    
    with open(output_filename, 'w') as f:
        f.write("REMARK   Generated by sim2iq_pre.py\n")
        f.write("REMARK   Radii are from precomputed unique volumes\n")
        
        for i in range(pdb.natoms):
            # Use the stored PDB line (first 54 chars) and append charge and radius
            pdb_line = pdb.pdb_lines[i]
            charge = 0.0
            radius = pdb.unique_radius[i]
            
            f.write(f"{pdb_line}{charge:7.4f} {radius:8.6f}\n")
        
        f.write("END\n")


def main():
    parser = argparse.ArgumentParser(description="Use precomputed unique atomic volumes and output as PQR")
    parser.add_argument("-f", "--file", type=str, required=True, help="Input PDB file")
    parser.add_argument("-v", "--volumes", type=str, help="JSON file containing precomputed volumes (default: ../data/mean_volumes.json)")
    parser.add_argument("-o", "--output", type=str, help="Output PQR filename (default: input_uniqueradii.pqr)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.file):
        print(f"Error: Input file {args.file} not found")
        return 1
    
    # Set volumes file
    if args.volumes:
        volumes_file = args.volumes
    else:
        # Look for mean_volumes.json in ../data/ relative to script location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        volumes_file = os.path.join(script_dir, "..", "data", "mean_volumes.json")
    
    if not os.path.exists(volumes_file):
        print(f"Error: Volume file {volumes_file} not found")
        return 1
    
    # Set output filename
    if args.output:
        output_filename = args.output
    else:
        base_name = os.path.splitext(args.file)[0]
        output_filename = f"{base_name}_uniqueradii.pqr"
    
    # Load PDB
    pdb = PDB(args.file)
    print(f"Loaded PDB with {pdb.natoms} atoms")
    
    # Load precomputed volumes
    if not pdb.load_unique_volumes(volumes_file):
        print("Failed to load unique volumes")
        return 1
    
    # Write output with unique radii
    write_pqr(pdb, output_filename)
    
    print(f"Output written to {output_filename}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
