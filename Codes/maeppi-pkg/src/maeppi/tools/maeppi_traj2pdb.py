#!/usr/bin/env python3
"""
euler2pdb.py - Convert trajectory files to PDB structures

Usage:
  euler2pdb.py template.pdb trj [--onlyCA] [--exyz] > output
"""

import sys
import argparse
import numpy as np
from typing import List, Tuple, Dict, Optional

# Constants
MAXLENLINE = 1000
CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

class Atom:
    def __init__(self, line: str):
        """Parse an ATOM/HETATM line from PDB format"""
        self.line = line.rstrip('\n')
        
        # Store the original line parts
        self.prefix = line[:30]  # First 30 characters
        self.suffix = line[54:] if len(line) > 54 else ""  # Everything after coordinates
        
        # Extract only what we need
        try:
            self.x = float(line[30:38])
            self.y = float(line[38:46])
            self.z = float(line[46:54])
        except ValueError:
            self.x = self.y = self.z = 0.0
        
        # Extract chain ID at position 21
        self.chain_id = line[21:22] if len(line) > 53 else " "
        
        # Extract element at position 13 for exyz format
        self.element = line[13:14] if len(line) > 53 else " "
        
        # Extract atom name to identify CA atoms
        self.atom_name = line[12:16].strip() if len(line) > 53 else ""
    
    def to_pdb_line(self, new_chain_id: str = None) -> str:
        """Return the atom in PDB format with optional modifications"""
        # Start with the original prefix
        prefix = list(self.prefix)
        
        # Change chain ID if needed
        if new_chain_id is not None and len(prefix) > 21:
            prefix[21] = new_chain_id
        
        # Combine with the updated coordinates and original suffix
        return ''.join(prefix) + f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}" + self.suffix

def euler_to_rot(psi: float, theta: float, phi: float) -> np.ndarray:
    """Convert Euler angles to rotation matrix"""
    r11 = np.cos(psi) * np.cos(phi) - np.sin(psi) * np.cos(theta) * np.sin(phi)
    r21 = np.sin(psi) * np.cos(phi) + np.cos(psi) * np.cos(theta) * np.sin(phi)
    r31 = np.sin(theta) * np.sin(phi)

    r12 = -np.cos(psi) * np.sin(phi) - np.sin(psi) * np.cos(theta) * np.cos(phi)
    r22 = -np.sin(psi) * np.sin(phi) + np.cos(psi) * np.cos(theta) * np.cos(phi)
    r32 = np.sin(theta) * np.cos(phi)

    r13 = np.sin(psi) * np.sin(theta)
    r23 = -np.cos(psi) * np.sin(theta)
    r33 = np.cos(theta)

    rot = np.array([
        [r11, r12, r13],
        [r21, r22, r23],
        [r31, r32, r33]
    ])
    
    return rot

def read_pdb_file(filename: str, only_ca: bool = False) -> List[Atom]:
    """Read a PDB file and return a list of Atom objects"""
    atoms = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = Atom(line)
                # Skip if only_ca is True and this is not a CA atom
                if only_ca and atom.atom_name != "CA":
                    continue
                atoms.append(atom)
    
    return atoms

def count_atoms_in_traj(traj_file: str, is_txt: bool) -> int:
    """Count the number of atoms/transformations in the trajectory file"""
    count = 0
    prev_id = 0
    
    with open(traj_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if is_txt:
                if line.startswith("MODEL") or line.startswith("ENDMDL"):
                    continue
                
                parts = line.split()
                if len(parts) >= 7:  # Expect at least 7 columns in text format
                    try:
                        atom_id = int(parts[0])
                        if atom_id > prev_id:
                            prev_id = atom_id
                            count += 1
                        else:
                            break
                    except ValueError:
                        continue
            else:  # PDB format
                if line.startswith("ATOM"):
                    try:
                        # Extract atom ID from positions 6-11
                        atom_id = int(line[6:11])
                        if atom_id > prev_id:
                            prev_id = atom_id
                            count += 1
                        else:
                            break
                    except ValueError:
                        continue
    
    return count

def read_transformations(traj_file: str, num_atoms: int, is_txt: bool) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Read transformation and rotation data from trajectory file"""
    transformations = []
    
    with open(traj_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if is_txt:
                if line.startswith("MODEL") or line.startswith("ENDMDL"):
                    continue
                
                parts = line.split()
                if len(parts) >= 7:  # Expect at least 7 columns in text format
                    try:
                        # Format: ID TR1 TR2 TR3 XR1 XR2 XR3
                        tr = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
                        xr = np.array([float(parts[4]), float(parts[5]), float(parts[6])])
                        transformations.append((tr, xr))
                    except (ValueError, IndexError):
                        continue
            else:  # PDB format
                if line.startswith("ATOM"):
                    try:
                        # Format based on PD6DReadline function
                        # sscanf(line+30, "%8lf%8lf%8lf%8lf%8lf%8lf", &tr[0], &tr[1], &tr[2], &xr[0], &xr[1], &xr[2]);
                        values = [line[30:38],line[38:46],line[46:54],line[54:62],line[62:70],line[70:78]]
                        if len(values) >= 6:
                            tr = np.array([float(values[0]), float(values[1]), float(values[2])])
                            xr = np.array([float(values[3]), float(values[4]), float(values[5])])
                            transformations.append((tr, xr))
                    except (ValueError, IndexError):
                        continue
            
            if len(transformations) >= num_atoms:
                break
    
    return transformations

def apply_transformation(atoms: List[Atom], tr: np.ndarray, xr: np.ndarray) -> List[Atom]:
    """Apply transformation and rotation to atoms"""
    # Convert Euler angles to rotation matrix
    rot_matrix = euler_to_rot(xr[0], xr[1], xr[2])
    
    transformed_atoms = []
    for atom in atoms:
        # Create a copy of the atom
        new_atom = Atom(atom.to_pdb_line())
        
        # Get original coordinates
        orig_coords = np.array([atom.x, atom.y, atom.z])
        
        # Apply rotation and translation
        new_coords = tr + np.dot(orig_coords, rot_matrix.T)
        
        # Update coordinates
        new_atom.x = new_coords[0]
        new_atom.y = new_coords[1]
        new_atom.z = new_coords[2]
        
        transformed_atoms.append(new_atom)
    
    return transformed_atoms

def write_exyz(atoms: List[Atom]) -> str:
    """Write atoms in exyz format using the element from position 13 of the PDB line"""
    output = []
    
    for atom in atoms:
        # Use element (char at position 13) as specified
        output.append(f"{atom.element} {atom.x:16.6f}{atom.y:16.6f}{atom.z:16.6f}")
    
    return "\n".join(output)

def main():
    parser = argparse.ArgumentParser(description='Convert trajectory files to PDB structures')
    parser.add_argument('template', help='Template PDB file')
    parser.add_argument('traj', help='Trajectory file')
    parser.add_argument('--onlyCA', action='store_true', help='Output only CA atoms')
    parser.add_argument('--exyz', action='store_true', help='Output in exyz format')
    parser.add_argument('--txt', action='store_true', help='Trajectory is in text format')
    
    args = parser.parse_args()
    
    # Read template PDB file
    template_atoms = read_pdb_file(args.template, args.onlyCA)
    
    # Count atoms in trajectory file
    num_transformations = count_atoms_in_traj(args.traj, args.txt)
    #print(f"Found {num_transformations} transformations in trajectory file", file=sys.stderr)
    
    if num_transformations == 0:
        print("Error: No transformations found in trajectory file", file=sys.stderr)
        sys.exit(1)
    
    # Read transformations
    transformations = read_transformations(args.traj, num_transformations, args.txt)
    #print(f"Read {len(transformations)} transformations", file=sys.stderr)
    
    # Generate all structures
    all_atoms = []
    
    for i, (tr, xr) in enumerate(transformations):
        # Get chain ID from cycle through CHAIN_IDS
        chain_id = CHAIN_IDS[i % len(CHAIN_IDS)]
        
        # Apply transformation
        transformed_atoms = apply_transformation(template_atoms, tr, xr)
        
        # Set chain ID for all atoms in this transformation
        for atom in transformed_atoms:
            atom.chain_id = chain_id
        
        all_atoms.extend(transformed_atoms)
    
    # Output
    if args.exyz:
        print(write_exyz(all_atoms))
    else:
        # Write PDB header
        #print("HEADER    GENERATED BY EULER2PDB.PY")
        #print("REMARK    TRANSFORMATIONS APPLIED: ", num_transformations)
        
        # Write all atoms
        atom_num = 1
        for atom in all_atoms:
            # Set atom number and chain ID without processing the rest of the line
            line = atom.to_pdb_line(new_chain_id=atom.chain_id)
            print(line,end="")
             
        # Write PDB footer
        print("END")

if __name__ == "__main__":
    main()
