#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FMAPIq Pipeline - Core Module
"""

import os
import sys
import shutil
import tempfile
import subprocess
import numpy as np
import glob
import tarfile
import json
from concurrent.futures import ThreadPoolExecutor
from io import StringIO

# Optional plotting
MATPLOTLIB_AVAILABLE = False
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    pass

# Tool imports
try:
    from fmapb2.tools.fmappre import main as prep_main
    from fmapb2.tools.fmap import main as fmap_main
    from fmapb2.tools.fmapstat import main as stat_main
    from sim2iq.tools.sim2iq_pre import main as iq_prep_main
    from sim2iq.tools.sim2iq_calc import main as iq_calc_main
    from maeppi.tools.maeppi_sim import main as sim_main
    from maeppi.tools.maeppi_shm_server import main as server_main
    import maeppi
    ANGLE_DATA = os.path.join(maeppi.__path__[0], "data/sspNF%05d/rot.eul" % 72)
except ImportError as e:
    print(f"Error importing tools: {e}")
    print('Go to parent directory and make; then go back use "uv sync --reinstall" to update.')
    sys.exit(1)

# Default parameters
DEFAULTS = {
    'n_proteins': 200,
    'box_size': 619.65,
    'ion': 0.016,
    'temp': 25,
    'nrot': 72,
    'escl': 1.0,
    'vscl': -1.0,
    'n_reps': 3,
    'n_cycles': 20000,
    'spacing': 1.0,
    # Plot defaults
    'show_individual': True,
    'show_pq': True,
    'q_min': 0.01,
    'q_max': 0.4,
    'y_min': 0.01,
    'y_max': 1.5,
    'log_scale_x': True,
    'log_scale_y': False,
}

class Pipeline:
    """Main pipeline class for FMAPIq workflow"""

    def __init__(self, work_dir=None):
        self.work_dir = work_dir
        if self.work_dir:
            self._create_work_dir()

        # Pipeline state
        self.pqr_file = None
        self.temp = None
        self.mass = None
        self.params = None
        self.configs = []
        self.iq_data = None
        self.config_data = DEFAULTS.copy()

        # Thread settings
        self.threads = os.cpu_count() // 2
        self.sim_threads = 1

    def _get_config_file_path(self):
        """Get path to configuration JSON file"""
        if self.work_dir:
            self._create_work_dir()
            return os.path.join(self.work_dir, "config.json")
        return None

    def save_config(self):
        """Save current configuration to JSON file"""
        config_file = self._get_config_file_path()
        if config_file:
            try:
                with open(config_file, 'w') as f:
                    json.dump(self.config_data, f, indent=2)
            except Exception as e:
                print(f"Warning: Could not save config: {e}")

    def load_config(self):
        """Load configuration from JSON file"""
        config_file = self._get_config_file_path()
        if config_file and os.path.exists(config_file):
            try:
                with open(config_file, 'r') as f:
                    loaded_config = json.load(f)
                    self.config_data.update(loaded_config)
                    return True
            except Exception as e:
                print(f"Warning: Could not load config: {e}")
        return False

    def update_config(self, **kwargs):
        """Update configuration with new values and save"""
        self.config_data.update(kwargs)
        self.save_config()

    def _set_openmp_threads(self, threads):
        """Configure OpenMP threading"""
        os.environ['OMP_NUM_THREADS'] = str(threads)

        # Attempt to dynamically set OpenMP threads
        try:
            import ctypes
            openmp_libs = ['libgomp.so.1', 'libiomp5.so', 'libomp.so']
            for lib_name in openmp_libs:
                try:
                    libgomp = ctypes.CDLL(lib_name, mode=ctypes.RTLD_GLOBAL)
                    if hasattr(libgomp, 'omp_set_num_threads'):
                        libgomp.omp_set_num_threads(int(threads))
                        print(f"Set OpenMP threads to {threads} via {lib_name}")
                        break
                except OSError:
                    continue
        except Exception as e:
            print(f"Could not set OpenMP threads via ctypes: {e}")

    def _create_work_dir(self):
        """Create work directory if it doesn't exist"""
        if not self.work_dir:
            tmp_dir = "tmp"
            if not os.path.exists(tmp_dir):
                os.mkdir(tmp_dir)
            self.work_dir = tempfile.mkdtemp(prefix="maeppi_", dir=tmp_dir)
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)
        return self.work_dir

    def _unpack_tgz(self, tgz_path):
        """Unpack .tgz archive to work directory"""
        if not os.path.exists(tgz_path):
            return None, "Archive file not found"

        try:
            self._create_work_dir()

            with tarfile.open(tgz_path, "r:gz") as tar:
                tar.extractall(self.work_dir)

            return self.work_dir, f"Unpacked archive to: {self.work_dir}"

        except Exception as e:
            return None, f"Failed to unpack archive: {str(e)}"

    def _parse_params_file(self):
        """Parse parameters from parms.txt file"""
        if not self.params:
            return

        try:
            lines = self.params.split('\n')
            for line in lines:
                if 'kbt' in line:
                    self.temp = float(line.split()[2])
                elif 'mass' in line:
                    self.mass = float(line.split()[2])
        except (IndexError, ValueError) as e:
            print(f"Warning: Could not parse parameters: {e}")

    def load_data(self, stage=None):
        """Load existing data from work directory"""
        if not self.work_dir or not os.path.exists(self.work_dir):
            return False

        # Load configuration first
        self.load_config()

        try:
            # Stage 1: Parameters and energy maps
            if stage in [None, 1]:
                params_file = os.path.join(self.work_dir, "parms.txt")
                if os.path.exists(params_file):
                    self.pqr_file = os.path.join(self.work_dir, "subA.pqr")
                    with open(params_file, "r") as f:
                        self.params = f.read()
                    self._parse_params_file()

            # Stage 2: Simulation configurations
            if stage in [None, 2]:
                config_files = glob.glob(os.path.join(self.work_dir, "config_*.dat"))
                self.configs = [os.path.basename(c) for c in config_files]

            # Stage 3: I(q) data
            if stage in [None, 3]:
                iqs_file = os.path.join(self.work_dir, "Iqs.txt")
                if os.path.exists(iqs_file):
                    data = np.loadtxt(iqs_file)
                    if len(data) > 1:
                        self.iq_data = {
                            'q': data[:, 0],
                            'pq': data[:, 1],
                            'curves': data[:, 2:],
                            'n_curves': data.shape[1] - 2
                        }
            return True
        except Exception as e:
            print(f"Load error: {e}")
            return False

    def stage1(self, pqr_file, ion, temp, escl, vscl, threads, nrot=72):
        """Stage 1: Prepare energy maps"""
        if not (pqr_file or self.pqr_file):
            return "Error: No PQR file provided"

        # Update configuration
        self.update_config(ion=ion, temp=temp, escl=escl, vscl=vscl, nrot=nrot)

        try:
            self._create_work_dir()

            if not self.pqr_file:
                self.pqr_file = os.path.join(self.work_dir, "subA.pqr")
                shutil.copy2(pqr_file, self.pqr_file)

            original_dir = os.getcwd()
            os.chdir(self.work_dir)

            try:
                # Run fmappre
                vscl_arg = "-" if vscl < 0.0 else str(vscl)
                result = prep_main(['-o', "parms.txt", "subA.pqr", str(ion), str(temp), str(escl), vscl_arg])
                if result is not None and result != 0:
                    return "Error: fmappre failed"

                # Append nrot to parms.txt if nrot != 72
                if nrot != 72:
                    with open("parms.txt", "a") as f:
                        f.write("%-12s : %16.0f\n" % ("nrot", nrot))
                    ANGLE_DATA_N = os.path.join(maeppi.__path__[0], "data/sspNF%05d/rot.eul" % nrot)
                else:
                    ANGLE_DATA_N = ANGLE_DATA

                # Run fmap
                self._set_openmp_threads(threads)


                # Redirect stdout to capture fmap output
                original_stdout = os.dup(1)
                with open('fmap.ern', 'w') as temp_file:
                    os.dup2(temp_file.fileno(), 1)
                    try:
                        result = fmap_main(['subA.vdw', 'parms.txt', ANGLE_DATA_N])
                        if result is not None and result != 0:
                            return "Error: fmap failed"
                    finally:
                        os.dup2(original_stdout, 1)
                        os.close(original_stdout)

                # Run fmapstat
                result = stat_main(['stats'])
                if result is not None and result != 0:
                    return "Error: fmapstat failed"

            finally:
                os.chdir(original_dir)

            self.load_data(stage=1)
            return f"Stage 1 Complete!\n{self.params}"

        except Exception as e:
            return f"Stage 1 Error: {str(e)}"

    def _run_single_simulation(self, args):
        """Execute single simulation job"""
        sim_id, n_proteins, box_size, temp, n_cycles, server_key = args

        try:
            config_file = f"config_{sim_id:03d}.dat"
            n_steps = int(n_cycles) * int(n_proteins)

            result = sim_main([
                str(n_proteins), str(box_size), str(temp), str(n_steps),
                str(sim_id), '-', config_file, str(n_steps), server_key
            ])
            return config_file if result == 0 else None
        except Exception as e:
            print(f"Simulation {sim_id} failed: {e}")
            return None

    def stage2(self, n_reps, n_cycles, n_proteins, box_size, threads, sim_threads):
        """Stage 2: Run molecular dynamics simulations"""
        if not self.work_dir or not os.path.exists(self.work_dir):
            return "Error: Run Stage 1 first"

        # Update configuration
        self.update_config(n_reps=n_reps, n_cycles=n_cycles, n_proteins=n_proteins, box_size=box_size)

        try:
            original_dir = os.getcwd()
            os.chdir(self.work_dir)

            try:
                # Start simulation server
                if server_main(['start', 'parms.txt', '--load-only']) != 0:
                    return "Error: Failed to start simulation server"

                # Get server key
                original_stdout = sys.stdout
                output_buffer = StringIO()
                sys.stdout = output_buffer
                server_main(['status', 'parms.txt', '--key'])
                sys.stdout = original_stdout
                server_key = output_buffer.getvalue().strip().split('\n')[-1]

                # Prepare simulation arguments
                sim_args = [(i, n_proteins, box_size, self.temp, n_cycles, server_key)
                           for i in range(n_reps)]

                # Run simulations in parallel
                self._set_openmp_threads(sim_threads)
                max_jobs = int(threads) // int(sim_threads)
                with ThreadPoolExecutor(max_workers=max_jobs) as executor:
                    results = list(executor.map(self._run_single_simulation, sim_args))

                # Cleanup server
                server_main(['cleanup', 'parms.txt'])

            finally:
                os.chdir(original_dir)

            self.load_data(stage=2)
            successful_configs = [r for r in results if r is not None]
            return f"Stage 2 Complete!\nGenerated {len(successful_configs)} configurations"

        except Exception as e:
            return f"Stage 2 Error: {str(e)}"

    def stage3(self, box_size, spacing, threads):
        """Stage 3: Calculate I(q) scattering intensities"""
        if not self.work_dir or not self.configs:
            return "Error: Run Stage 2 first"

        # Update configuration
        self.update_config(spacing=spacing)

        try:
            n_proteins = int(self.config_data["n_proteins"])
            n_grid_points = int(box_size / spacing)
            original_dir = os.getcwd()
            os.chdir(self.work_dir)

            try:
                self._set_openmp_threads(threads)

                # Prepare input for I(q) calculation
                if iq_prep_main(['-f', 'subA.vdw', '-o', 'subA.vol']) != 0:
                    return "Error: sim2iq_pre failed"

                # Calculate P(q) reference
                pq_file = "Pq.dat"
                result = iq_calc_main(['subA.vol', str(n_grid_points), str(box_size), pq_file, '--steps', '4'])
                if not os.path.exists(pq_file):
                    return "Error: P(q) calculation failed"

                # Calculate I(q) for each configuration
                iq_results = []
                for i, config_file in enumerate(self.configs):
                    iq_output = f"iq_{i:03d}.dat"
                    iq_calc_main(['subA.vol', str(n_grid_points), str(box_size), iq_output,
                                   '--steps', '4', '--traj', config_file])
                    if os.path.exists(iq_output):
                        try:
                            data = np.loadtxt(iq_output)
                            iq_results.append({'q': data[:, 0], 'iq': data[:, 1]})
                        except Exception:
                            continue

                if not iq_results:
                    return "Error: No I(q) data generated"

                # Combine all I(q) results
                pq_data = np.loadtxt(pq_file)
                q_values = iq_results[0]['q']

                # Normalize I(q) curves by number of proteins
                iq_curves = []
                for result in iq_results:
                    if len(result['iq']) == len(q_values):
                        iq_curves.append(result['iq'] / n_proteins)

                if not iq_curves:
                    return "Error: No valid I(q) curves"

                # Save combined data: q, P(q), I(q) curves
                combined_data = np.column_stack([q_values, pq_data[:, 1], np.column_stack(iq_curves)])
                np.savetxt("Iqs.txt", combined_data, fmt='%.6e')

            finally:
                os.chdir(original_dir)

            self.load_data(stage=3)
            return f"Stage 3 Complete!\nCalculated I(q) for {len(iq_results)} configurations"

        except Exception as e:
            return f"Stage 3 Error: {str(e)}"

    def create_plot(self, show_individual=True, show_pq=True, q_min=0.01, q_max=0.4,
                   y_min=0, y_max=1.2, log_scale_x=True, log_scale_y=False):
        """Create I(q) plot with statistics"""
        if not self.iq_data:
            return None, "No I(q) data available"

        # Update configuration
        self.update_config(show_individual=show_individual, show_pq=show_pq,
                          q_min=q_min, q_max=q_max, y_min=y_min, y_max=y_max,
                          log_scale_x=log_scale_x, log_scale_y=log_scale_y)

        try:
            q = self.iq_data['q']
            pq = self.iq_data['pq']
            curves = self.iq_data['curves']

            # Normalize data
            pq_norm_factor = pq[0]
            pq_normalized = pq / pq_norm_factor
            curves_normalized = curves / pq_norm_factor

            # Calculate statistics
            mean_intensity = np.mean(curves_normalized, axis=1)
            std_error = np.std(curves_normalized, axis=1) / np.sqrt(curves_normalized.shape[1])

            # Save statistics file
            stats_data = np.column_stack([q, pq_normalized, mean_intensity, std_error])
            stats_file = os.path.join(self.work_dir, "Iqs.stat.txt")
            np.savetxt(stats_file, stats_data, fmt='%.6e')

            if not MATPLOTLIB_AVAILABLE:
                return None, "Matplotlib not available for plotting"

            # Create the plot
            plt.figure(figsize=(10, 6))

            if show_pq:
                plt.plot(q, pq_normalized, 'lightgray', alpha=0.7, linewidth=1, label='P(q)')

            if show_individual:
                for i in range(curves_normalized.shape[1]):
                    plt.plot(q[1:], curves_normalized[1:, i], 'red', alpha=0.3, linewidth=0.5)

            plt.errorbar(q[1:], mean_intensity[1:], yerr=std_error[1:],
                        label=f'Mean ± SEM (n={self.iq_data["n_curves"]})',
                        capsize=2, linewidth=2)

            plt.xlabel('q (Å⁻¹)')
            plt.ylabel('I(q) / n_proteins')
            plt.title('Small Angle Scattering Intensity I(q)')
            plt.legend()
            plt.grid(True, alpha=0.3)

            plt.xlim(q_min, q_max)
            plt.ylim(y_min, y_max)

            if log_scale_x:
                plt.xscale('log')
            if log_scale_y:
                plt.yscale('log')

            plt.tight_layout()

            plot_path = os.path.join(self.work_dir, "iq_plot.png")
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()

            return plot_path, "Plot created successfully"

        except Exception as e:
            return None, f"Plot error: {str(e)}"

    def pack_results(self, skip_energy_maps=True, skip_raw_data=True):
        """Pack work directory into compressed archive"""
        if not self.work_dir or not os.path.exists(self.work_dir):
            return None, "No work directory to pack"

        try:
            work_dir_name = os.path.basename(self.work_dir)
            archive_path = os.path.join(self.work_dir, f"{work_dir_name}.tar.gz")

            with tarfile.open(archive_path, "w:gz") as tar:
                for root, dirs, files in os.walk(self.work_dir):
                    # Skip energy maps directory if requested
                    if skip_energy_maps and 've' in dirs:
                        dirs.remove('ve')

                    for filename in files:
                        # Skip the archive itself
                        if filename.endswith('.tar.gz'):
                            continue
                        # Skip raw data files if requested
                        if skip_raw_data and filename.endswith('.dat'):
                            continue

                        file_path = os.path.join(root, filename)
                        archive_name = os.path.relpath(file_path, self.work_dir)
                        tar.add(file_path, arcname=archive_name)

            size_mb = os.path.getsize(archive_path) / (1024 * 1024)
            return archive_path, f"Archive created: {os.path.basename(archive_path)} ({size_mb:.1f} MB)"

        except Exception as e:
            return None, f"Pack error: {str(e)}"

def main():
    """Command line interface for the pipeline"""
    import argparse

    parser = argparse.ArgumentParser(description="FMAPIq Pipeline - Small Angle Scattering Calculation")
    parser.add_argument('pqr_file', nargs='?', help='Input PQR structure file')
    parser.add_argument('--work-dir', help='Work directory (new or existing)')
    parser.add_argument('--load', help='Load existing work directory')
    parser.add_argument('--pack', action='store_true', help='Pack results into archive')
    parser.add_argument('--skip-energy-maps', action='store_true', help='Skip energy maps when packing')
    parser.add_argument('--skip-raw-data', action='store_true', help='Skip raw data when packing')

    # Threading parameters
    parser.add_argument('--threads', type=int, help='Total threads (default: CPU cores)')
    parser.add_argument('--sim-threads', type=int, default=2, help='Threads per simulation job (default: 2)')

    # Pipeline parameters - only add non-plotting parameters
    excluded_params = ['show_individual', 'show_pq', 'q_min', 'q_max', 'y_min', 'y_max', 'log_scale_x', 'log_scale_y']
    for key, default in DEFAULTS.items():
        if key not in excluded_params:
            parser.add_argument(f'--{key.replace("_", "-")}',
                              type=type(default), default=default,
                              help=f'{key} (default: {default})')

    args = parser.parse_args()

    # Initialize pipeline
    pipeline = Pipeline(args.work_dir or args.load)

    if args.threads:
        pipeline.threads = args.threads

    # Load existing data if work directory exists
    if pipeline.work_dir and os.path.exists(pipeline.work_dir):
        if pipeline.load_data():
            print(f"Loaded existing work directory: {pipeline.work_dir}")

            # Check stage completion status
            stage1_complete = os.path.exists(os.path.join(pipeline.work_dir, "parms.txt"))
            stage2_complete = len(glob.glob(os.path.join(pipeline.work_dir, "config_*.dat"))) > 0
            stage3_complete = os.path.exists(os.path.join(pipeline.work_dir, "Iqs.txt"))

            print(f"Stage 1: {'OK' if stage1_complete else 'Not completed'}")
            print(f"Stage 2: {'OK' if stage2_complete else 'Not completed'}")
            print(f"Stage 3: {'OK' if stage3_complete else 'Not completed'}")

            if stage3_complete and pipeline.iq_data:
                print(f"I(q) datasets available: {pipeline.iq_data['n_curves']}")

            # Handle load-only or pack-only operations
            if args.load and not args.pqr_file:
                if args.pack:
                    file_path, message = pipeline.pack_results(args.skip_energy_maps, args.skip_raw_data)
                    if file_path:
                        print(f"Results packed: {file_path}")
                    else:
                        print(f"Packing failed: {message}")
                return True
        else:
            print(f"Failed to load work directory: {pipeline.work_dir}")
            return False

    # Validate PQR file for new runs
    pqr_file = args.pqr_file
    if not pqr_file and pipeline.work_dir:
        potential_pqr = os.path.join(pipeline.work_dir, "subA.pqr")
        if os.path.exists(potential_pqr):
            pqr_file = potential_pqr
            print(f"Using existing PQR file: {pqr_file}")

    if not pqr_file:
        print("Error: PQR file required for new pipeline run")
        return False

    if not os.path.exists(pqr_file):
        print(f"Error: PQR file not found: {pqr_file}")
        return False

    # Determine which stages need to be executed
    stage1_needed = not (pipeline.work_dir and os.path.exists(os.path.join(pipeline.work_dir, "parms.txt")))
    stage2_needed = not (pipeline.work_dir and glob.glob(os.path.join(pipeline.work_dir, "config_*.dat")))
    stage3_needed = not (pipeline.work_dir and os.path.exists(os.path.join(pipeline.work_dir, "Iqs.txt")))

    # Execute pipeline stages
    if stage1_needed:
        print("Running Stage 1: Energy map calculation...")
        result = pipeline.stage1(pqr_file, args.ion, args.temp, args.escl, args.vscl, pipeline.threads, args.nrot)
        print(result)
        if "Error" in result:
            return False
    else:
        print("Stage 1: Already completed, skipping")

    if stage2_needed:
        print("Running Stage 2: Molecular dynamics simulations...")
        result = pipeline.stage2(args.n_reps, args.n_cycles, args.n_proteins, args.box_size,
                                 pipeline.threads, args.sim_threads)
        print(result)
        if "Error" in result:
            return False
    else:
        print("Stage 2: Already completed, skipping")

    if stage3_needed:
        print("Running Stage 3: I(q) intensity calculation...")
        result = pipeline.stage3(args.box_size, args.spacing, pipeline.threads)
        print(result)
        if "Error" in result:
            return False
    else:
        print("Stage 3: Already completed, skipping")

    print(f"Pipeline completed successfully! Results in: {pipeline.work_dir}")

    # Generate plot if matplotlib is available
    if MATPLOTLIB_AVAILABLE:
        plot_path, message = pipeline.create_plot()
        if plot_path:
            print(f"Plot saved: {plot_path}")
        else:
            print(f"Plot generation failed: {message}")

    # Pack results if requested
    if args.pack:
        file_path, message = pipeline.pack_results(args.skip_energy_maps, args.skip_raw_data)
        if file_path:
            print(f"Results packed: {file_path}")
        else:
            print(f"Packing failed: {message}")

    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
