#!/usr/bin/env python3
"""
Simple script to read final_config.dat and recalculate total energy.
"""

import sys

# Import from the package
try:
    import maeppi
except ModuleNotFoundError:
    # Fall back to local version, add the directory containing fmapb2
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    #print(sys.path)
    import maeppi

class Particle:
    """Simple particle class with position and orientation"""
    def __init__(self, pos=None, orient=None):
        self.coord = [0.0] * 6  # [x, y, z, theta1, theta2, theta3]
        if pos is not None:
            self.coord[:3] = pos
        if orient is not None:
            self.coord[3:] = orient

def read_configuration(filename):
    """Read particle configuration from file"""
    particles = []
    box_energy = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                # Try to extract energy from comment
                if 'Energy' in line:
                    try:
                        # Look for "Energy X.XXXXXX" pattern
                        parts = line.split()
                        for i, part in enumerate(parts):
                            if part == 'Energy' and i + 1 < len(parts):
                                box_energy = float(parts[i + 1].rstrip(','))
                                break
                    except:
                        pass
                continue
            
            if line:  # Non-empty line
                parts = line.split()
                if len(parts) >= 7:  # id x y z euler_x euler_y euler_z
                    coords = [float(x) for x in parts[1:7]]  # Skip particle ID
                    particle = Particle()
                    particle.coord = coords
                    particles.append(particle)
    
    return particles, box_energy

def calculate_total_energy(particles, box_len, ern_calc):
    """Calculate total energy of the configuration"""
    total_energy = 0.0
    n_particles = len(particles)
    
    for i in range(n_particles):
        for j in range(i + 1, n_particles):
            energy = ern_calc.energy(
                particles[i].coord[:3], particles[i].coord[3:],
                particles[j].coord[:3], particles[j].coord[3:],
                box_len, 'euler'
            )
            total_energy += energy
    
    return total_energy

def main():
    """Main function to recalculate energy"""
    if len(sys.argv) < 4:
        print("Usage: python maeppi_ern.py <config_file> <box_length> <key_string>")
        print("Example: python maeppi_ern.py final_config.dat 10.0 some_key")
        sys.exit(1)
    
    config_file = sys.argv[1]
    box_length = float(sys.argv[2])
    key_string = sys.argv[3]
    
    # Box dimensions
    box_len = [box_length, box_length, box_length]
    
    print(f"Reading configuration from: {config_file}")
    print(f"Box length: {box_length}")
    
    # Read configuration
    try:
        particles, original_energy = read_configuration(config_file)
        print(f"Read {len(particles)} particles")
        
        if original_energy is not None:
            print(f"Original energy from file: {original_energy:.6f}")
        
    except FileNotFoundError:
        print(f"Error: Configuration file '{config_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading configuration: {e}")
        sys.exit(1)
    
    if not particles:
        print("No particles found in configuration file")
        sys.exit(1)
    
    # Initialize energy calculator
    ern_calc = maeppi.Ern(key_string)
    
    # Recalculate total energy
    print("Recalculating total energy...")
    recalc_energy = calculate_total_energy(particles, box_len, ern_calc)
    
    print("-" * 50)
    print(f"Recalculated total energy: {recalc_energy:.6f}")
    
    if original_energy is not None:
        diff = abs(recalc_energy - original_energy)
        print(f"Difference from original: {diff:.6f}")
        if diff < 1e-6:
            print("✓ Energies match (within numerical precision)")
        else:
            print("⚠ Energies differ significantly")
    
    print(f"Number of particle pairs: {len(particles) * (len(particles) - 1) // 2}")

if __name__ == "__main__":
    main()
