#!/usr/bin/env python3
"""
Simple CFFI wrapper for SAXS calculator
"""

import argparse
import numpy as np
import sys

try:
    import sim2iq
except ModuleNotFoundError:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    import sim2iq

def main(argv=None):
    """Simple command line interface"""
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser(description='Calculate SAXS profile from PQR file')
    parser.add_argument('input_pqr', help='Input PQR file')
    parser.add_argument('ngrid', type=int, help='Number of grid points')
    parser.add_argument('box_length', type=float, help='Box length (Å)')
    parser.add_argument('output_dat', help='Output data file')
    
    parser.add_argument('-b', '--b-factor', type=float, default=0.0, 
                       help='B-factor (default: auto)')
    parser.add_argument('--offset', type=float, default=2.8,
                       help='Cutoff offset (default: 2.8)')
    parser.add_argument('--steps', type=int, default=1, choices=[1,2,3,4],
                       help='Number of calculation steps (default: 1)')
    parser.add_argument('--se', type=float, default=-1.0,
                       help='Solvent excluded weight (default: -1.0)')
    parser.add_argument('--shell', type=float, default=0.011,
                       help='Shell weight (default: 0.011)')
    parser.add_argument('--traj', help='Trajectory file')
    
    args = parser.parse_args(argv)
    
    # Create calculator and run
    calc = sim2iq.Calc()
    spacing = args.box_length / args.ngrid
    q, intensity = calc(
        pqr_file=args.input_pqr,
        grid_size=args.ngrid,
        grid_spacing=spacing,
        steps=args.steps,
        cutoff_offset=args.offset,
        se_weight=args.se,
        ou_weight=args.shell,
        in_weight=-args.shell,
        b_factor=args.b_factor,
        traj_file=args.traj,
        traj_is_txt=True
    )
    
    # Write output
    np.savetxt(args.output_dat, np.column_stack([q, intensity]),  # Fixed: now uses output_dat
               fmt='%.12e', delimiter='\t')
    
    print(f"Wrote {len(q)} data points to {args.output_dat}")

if __name__ == '__main__':
    main()
