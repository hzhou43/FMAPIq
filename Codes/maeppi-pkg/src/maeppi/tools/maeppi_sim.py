#!/usr/bin/env python3
"""
Python wrapper for Monte Carlo simulation using CFFI.
This provides a high-level interface to the C simulation engine.
"""

import numpy as np
import os
from typing import Optional, Tuple, List
import time

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


class MCSimulation:
    """Monte Carlo simulation wrapper using maeppi interface."""
    
    def __init__(self):
        self._sim = None
    
    def create_simulation(self, n_particles: int, box_length: float, temperature: float,
                         max_trans_param: float = 1.0, max_rot_param: float = 1.0,
                         seed: int = 12345, fmap_key: str = "abc123", zratio: float = 1.0):
        """Create a new simulation box."""
        if self._sim is not None:
            del self._sim
        
        self._sim = maeppi.Sim(n_particles, box_length, temperature,
                              max_trans_param, max_rot_param, seed, fmap_key, zratio)
        
        print(f"Created simulation box with {n_particles} particles")
        return self
    
    def initialize_configuration(self, start_file: str = "-") -> bool:
        """Initialize particle configuration from file or randomly."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        if start_file == "-":
            print("Using random initialization...")
            return bool(maeppi.sim.sim_initialize_random(self._sim._box))
        else:
            print(f"Attempting to load configuration from {start_file}...")
            success = bool(maeppi.sim.sim_load_configuration(self._sim._box, start_file.encode('utf-8')))
            if success:
                print(f"Successfully loaded configuration from {start_file}")
                return True
            else:
                print(f"Failed to load from {start_file}, falling back to random initialization...")
                return bool(maeppi.sim.sim_initialize_random(self._sim._box))
    
    def save_configuration(self, filename: str) -> bool:
        """Save current particle configuration to file."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        return bool(maeppi.sim.sim_save_configuration(self._sim._box, filename.encode('utf-8')))
    
    def run_steps(self, n_steps: int):
        """Run simulation for specified number of Monte Carlo steps."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        maeppi.sim.sim_run_steps(self._sim._box, n_steps)
    
    def get_stats(self) -> dict:
        """Get current simulation statistics."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        energy, moves, accepted, acceptance = self._sim.get_stats()
        
        return {
            'total_energy': energy,
            'n_moves': moves,
            'n_accepted': accepted,
            'acceptance_rate': acceptance
        }
    
    def get_particle_coords(self, particle_id: int) -> np.ndarray:
        """Get coordinates of a specific particle."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        coords = maeppi.sim_ffi.new("double[6]")
        result = maeppi.sim.sim_get_particle_coords(self._sim._box, particle_id, coords)
        
        if not result:
            raise ValueError(f"Invalid particle ID: {particle_id}")
        
        return np.array([coords[i] for i in range(6)])
    
    def set_particle_coords(self, particle_id: int, coords: np.ndarray) -> bool:
        """Set coordinates of a specific particle."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        if len(coords) != 6:
            raise ValueError("coords must have length 6")
        
        coords_c = maeppi.sim_ffi.new("double[6]", coords.tolist())
        result = maeppi.sim.sim_set_particle_coords(self._sim._box, particle_id, coords_c)
        return bool(result)
    
    def get_all_coords(self) -> np.ndarray:
        """Get coordinates of all particles."""
        params = self.get_parameters()
        n_particles = params['n_particles']
        
        coords = np.zeros((n_particles, 6))
        for i in range(n_particles):
            coords[i] = self.get_particle_coords(i)
        
        return coords
    
    def set_all_coords(self, coords: np.ndarray):
        """Set coordinates of all particles."""
        params = self.get_parameters()
        n_particles = params['n_particles']
        
        if coords.shape != (n_particles, 6):
            raise ValueError(f"coords must have shape ({n_particles}, 6)")
        
        for i in range(n_particles):
            self.set_particle_coords(i, coords[i])
    
    def get_box_dimensions(self) -> np.ndarray:
        """Get simulation box dimensions."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        box_len = maeppi.sim_ffi.new("double[3]")
        maeppi.sim.sim_get_box_dimensions(self._sim._box, box_len)
        
        return np.array([box_len[i] for i in range(3)])
    
    def get_parameters(self) -> dict:
        """Get simulation parameters."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        n_particles = maeppi.sim_ffi.new("int *")
        temperature = maeppi.sim_ffi.new("double *")
        max_trans_param = maeppi.sim_ffi.new("double *")
        max_rot_param = maeppi.sim_ffi.new("double *")
        
        maeppi.sim.sim_get_parameters(self._sim._box, n_particles, temperature, 
                                     max_trans_param, max_rot_param)
        
        return {
            'n_particles': n_particles[0],
            'temperature': temperature[0],
            'max_trans_param': max_trans_param[0],
            'max_rot_param': max_rot_param[0]
        }
    
    def run_with_trajectory(self, n_steps: int, sample_freq: int = 1000) -> Tuple[List[dict], np.ndarray]:
        """Run simulation and collect trajectory data."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")
        
        params = self.get_parameters()
        n_particles = params['n_particles']
        
        # Calculate number of samples
        n_samples = (n_steps + sample_freq - 1) // sample_freq + 1
        
        # Pre-allocate trajectory array
        trajectory = np.zeros((n_samples, n_particles, 6))
        stats_list = []
        
        print(f"Running {n_steps} steps with trajectory collection...")
        print(f"Sampling every {sample_freq} steps ({n_samples} total samples)")
        
        start_time = time.time()
        sample_idx = 0

        # Collect initial configuration
        initial_stats = self.get_stats()
        initial_stats['step'] = 0
        stats_list.append(initial_stats)
        trajectory[sample_idx] = self.get_all_coords()
        sample_idx += 1
        
        for step in range(0, n_steps, sample_freq):
            # Determine how many steps to run
            steps_to_run = min(sample_freq, n_steps - step)
            
            # Run simulation steps
            self.run_steps(steps_to_run)
            
            # Collect statistics
            stats = self.get_stats()
            stats['step'] = step + steps_to_run
            stats_list.append(stats)
            
            # Collect trajectory
            trajectory[sample_idx] = self.get_all_coords()
            sample_idx += 1
            
            # Progress update
            if (step + steps_to_run) % (n_steps // 10) == 0 or step + steps_to_run == n_steps:
                elapsed = time.time() - start_time
                progress = (step + steps_to_run) / n_steps
                eta = elapsed / progress - elapsed if progress > 0 else 0
                
                print(f"Step {step + steps_to_run:8d}/{n_steps}: "
                      f"Energy = {stats['total_energy']:12.6f}, "
                      f"Accept = {stats['acceptance_rate']:.3f}, "
                      f"ETA = {eta:.1f}s")
        
        total_time = time.time() - start_time
        print(f"Simulation completed in {total_time:.2f} seconds")
        
        return stats_list, trajectory[:sample_idx]

    def save_trajectory(self, filename: str, stats_list: List[dict], trajectory: np.ndarray, sel: List = None):
        """Save trajectory data to file with C-compatible format plus stats in comments."""
        if self._sim is None:
            raise RuntimeError("No simulation box created")

        if len(stats_list) != len(trajectory):
            raise ValueError("stats_list and trajectory must have same length")

        if sel is None:
            sel = list(range(len(trajectory)))  
        print(f"Saving trajectory with {len(sel)} frames to {filename}")

        with open(filename, 'w') as f:
            # Write file header
            params = self.get_parameters()
            box_dims = self.get_box_dimensions()
            f.write(f"# Particles: {params['n_particles']} ")
            f.write(f"Temperature: {params['temperature']:.6f} ")
            f.write(f"Box dimensions: {box_dims[0]:.3f} x {box_dims[1]:.3f} x {box_dims[2]:.3f}\n")

            for frame_idx, (stats, coords) in enumerate(zip(stats_list, trajectory)):
                if frame_idx not in sel:
                    continue
                # Write frame header with statistics
                f.write(f"# Frame {frame_idx}: Step {stats['step']}, Energy {stats['total_energy']:.6f}, "
                       f"Accept {stats['acceptance_rate']:.6f}, Moves {stats['n_moves']}\n")

                # Write particle coordinates in C-compatible format
                for particle_idx, particle_coords in enumerate(coords):
                    # Format matches C PARTICLE_SAVE: x y z euler_x euler_y euler_z
                    f.write(f"{particle_idx+1:16d}"
                            f"{particle_coords[0]:16.10f}{particle_coords[1]:16.10f}{particle_coords[2]:16.10f}"
                           f"{particle_coords[3]:16.10f}{particle_coords[4]:16.10f}{particle_coords[5]:16.10f}\n")

        print(f"Trajectory successfully saved to {filename}")
        print(f"File contains {len(sel)} frames with {trajectory.shape[1]} particles each")

def main():
    """Main function matching simulation.c argument structure."""
    import sys
    
    if len(sys.argv) < 10:
        print(f"Usage: {sys.argv[0]} <n_particles> <box_length> <temperature> <n_steps> <seed> <start.dat> <stop.dat> <sample_freq> <keyp> [tx] [rx] [zx]")
        print("  start.dat: '-' for random initialization, or filename to read initial configuration")
        print("  stop.dat:  filename to save final configuration")
        print("  keyp:      shared memory key for FMAP data (hex string, e.g., 'abc123')")
        return 1
    
    n_particles = int(sys.argv[1])
    box_length = float(sys.argv[2])
    temperature = float(sys.argv[3])
    n_steps = int(sys.argv[4])
    seed = int(sys.argv[5])
    start_file = sys.argv[6]
    stop_file = sys.argv[7]
    sample_freq = int(sys.argv[8])
    key_string = sys.argv[9]
    
    max_trans_param = float(sys.argv[10]) if len(sys.argv) > 10 else 1.0
    max_rot_param = float(sys.argv[11]) if len(sys.argv) > 11 else 1.0
    zratio = float(sys.argv[12]) if len(sys.argv) > 12 else 1.0
    
    print("Starting NVT Monte Carlo simulation")
    print(f"Parameters: {n_particles} particles, box={box_length}, T={temperature}")
    print(f"Start file: {start_file}")
    print(f"Stop file: {stop_file}")
    print("-" * 40)
    
    # Create simulation instance
    sim = MCSimulation()
    
    # Create simulation box
    sim.create_simulation(
        n_particles=n_particles,
        box_length=box_length,
        temperature=temperature,
        max_trans_param=max_trans_param,
        max_rot_param=max_rot_param,
        seed=seed,
        fmap_key=key_string,
        zratio=zratio
    )
    
    # Initialize configuration
    if not sim.initialize_configuration(start_file):
        print("Failed to initialize configuration")
        return 1

    # Get and display initial energy
    initial_stats = sim.get_stats()
    print(f"Initial energy: {initial_stats['total_energy']:.6f}")
    print(f"Temperature: {temperature:.3f}")
    
    params = sim.get_parameters()
    box_dims = sim.get_box_dimensions()
    print(f"Number of particles: {params['n_particles']}")
    print(f"Box dimensions: {box_dims[0]:.3f} x {box_dims[1]:.3f} x {box_dims[2]:.3f}")
    print("----------------------------------------")

    # Run simulation with trajectory collection
    start_time = time.time()
    
    for step in range(0, n_steps, sample_freq):
        steps_to_run = min(sample_freq, n_steps - step)
        sim.run_steps(steps_to_run)
        
        stats = sim.get_stats()
        elapsed = time.time() - start_time
        steps_per_sec = (step + steps_to_run) / elapsed if elapsed > 0 else 0
        
        print(f"Step {step + steps_to_run:8d}: Energy = {stats['total_energy']:12.6f}, "
              f"Accept = {stats['acceptance_rate']:.3f}, Speed = {steps_per_sec:.0f} steps/s")
    
    # Print final statistics
    final_stats = sim.get_stats()
    total_time = time.time() - start_time
    
    print("----------------------------------------")
    print("Simulation completed!")
    print(f"Final energy: {final_stats['total_energy']:.6f}")
    print(f"Overall acceptance rate: {final_stats['acceptance_rate']:.3f}")
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Average speed: {n_steps / total_time:.0f} steps/second")
    
    # Save configuration
    sim.save_configuration(stop_file)
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
