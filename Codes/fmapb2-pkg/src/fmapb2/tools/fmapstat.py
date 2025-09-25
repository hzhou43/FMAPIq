#!/usr/bin/env python

# Copyright 2014, Jerome Fung, Rebecca W. Perry, Thomas G. Dimiduk
#
# flyvbjerg_petersen_std_err is free software: you can redistribute it 
# and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the 
# License, or (at your option) any later version.
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with flyvbjerg_petersen_std_err.  If not, see 
# <http://www.gnu.org/licenses/>.

'''
Apply Flyvbjerg-Petersen block decorrelation method for estimating the 
standard error on the mean of a (possibly correlated) time series of 
values.

.. moduleauthor:: Jerome Fung <jerome.fung@gmail.com>
.. moduleauthor:: Rebecca W. Perry <rperry@seas.harvard.edu>
.. moduleauthor:: Tom Dimiduk <tom@dimiduk.net>

Reference: H. Flyvbjerg and H. G. Petersen, "Error estimates on correlated
data", J. Chem. Phys. 91, 461--466 (1989).
'''

import numpy as np
import warnings
import subprocess
import os
import sys
import gzip
import random
import math

def block_transformation(series):
    """
    Do a single step of fp block averaging.

    Parameters
    ----------
    series : ndarray
        Things we want to average: e.g. squared displacements to calculate the
        mean squared displacement of a Brownian particle.

    Returns
    -------
    blocked_series : ndarray
        an array of half the length of series with adjacent terms averaged

    Notes
    -----
    Flyvbjerg & Peterson 1989, equation 20

    """
    n_steps = series.size
    n_steps_p = int(np.floor(n_steps/2.))
    output = 0.5 * (series[::2][:n_steps_p] + series[1::2][:n_steps_p])
    return output

def calculate_blocked_variances(series, npmin = 15):
    """
    Compute a series of blocks and variances.

    Parameters
    ----------
    series : ndarray
        the thing we want to average: e.g. squared
        displacements for a Brownian random walk.
    npmin : int
        cutoff number of points to stop blocking

    Returns
    -------
    output_var, var_stderr : ndarray
        The variance and stderr of the variance at each blocking level

    Notes
    -----
    Flyvbjerg & Peterson suggest continuing blocking down to 2 points, but the
    last few blocks are very noisy, so we default to cutting off before that.

    """
    n_steps = series.size

    def var(d, n):
        # see eq. 27 of FP paper
        return d.var()/(n-1)
    def stderr_var(n):
        # see eq. 27 of FP paper
        return np.sqrt(2./(n-1))

    output_var = np.array([var(series, n_steps)]) # initialize
    var_stderr = np.array([stderr_var(n_steps)])

    while n_steps > npmin:
        series = block_transformation(series)
        n_steps = series.size
        # TODO: precompute size of output_var and var_stderr from n_steps
        # rather than appending
        output_var = np.append(output_var, var(series, n_steps))
        var_stderr = np.append(var_stderr, stderr_var(n_steps))

    return output_var, var_stderr

def detect_fixed_point(fp_var, fp_sev, full_output = False):
    """
    Find whether the block averages decorrelate the data series to a fixed
    point.

    Parameters
    ----------
    fp_var: ndarray
        FP blocked variance
    fp_sev: ndarray
        FP standard error of the variance.

    Returns
    -------
    best_var : float
        best estimate of the variance
    converged : bool
        did the series converge to a fixed point?
    bounds : (int, int) only if full_output is True
        range of fp_var averaged to compute best_var

    Notes
    -----
    Expects both fp_var and fp_sev will have been
    truncated to cut off points with an overly small n_p and
    correspondingly large standard error of the variance.

    """
    n_trans = fp_var.size # number of block transformations and index

    left_index = 0
    right_index = 0

    # Detect left edge
    for i in np.arange(n_trans-1):
        # ith point inside error bars of next point
        if np.abs(fp_var[i + 1] - fp_var[i]) < fp_var[i + 1] * fp_sev[i + 1]:
            left_index = i
            break

    # Check right edge
    for i in np.arange(n_trans)[::-1]:
        if np.abs(fp_var[i] - fp_var[i - 1]) < fp_var[i - 1] * fp_sev[i - 1]:
            right_index = i
            break

    # if search succeeds
    if (left_index >= 0) and (right_index >= 0) and \
            (right_index >= left_index):
        best_var = np.average(fp_var[left_index:right_index + 1],
                              weights = 1./fp_sev[left_index:right_index + 1])
        converged = True
    else:
        best_var = fp_var.max()
        converged = False

    if full_output is True:
        return best_var, converged, (left_index, right_index)
    else:
        return best_var, converged


def fp_stderr(data):
    '''
    Compute standard error using Flyvbjerg-Petersen blocking.

    Computes the standard error on the mean of a possibly correlated timeseries
    of measurements.

    Parameters
    ----------
    data: ndarray
        data whose mean is to be calculated, and for which we need
        a standard error on the mean

    Returns
    -------
    stderr : float
        Standard error on the mean of data

    Notes
    -----

    Uses the technique described in H. Flyvbjerg and H. G. Petersen,
    "Error estimates on correlated data", J. Chem. Phys. 91, 461--466 (1989).
    section 3.

    '''
    block_trans_var, block_trans_sev = calculate_blocked_variances(data)
    var_mean, conv, bounds = detect_fixed_point(block_trans_var,
                                                block_trans_sev, True)
    #print(bounds)

    if not conv:
        warnings.warn("Fixed point not found, returned value is a lower bound on the standard error")
    return np.sqrt(var_mean)

#  LocalWords:  Flyvbjerg

def getdata(fname,pos):
    data=[]
    if fname=="-":
        fp=sys.stdin
    else:
        fp=open(fname)
    for line in fp:
        if line[0]!="#":
            items=line.split()
            v=float(items[pos])
            data.append(v)
    return data

# New functions for stats.sh functionality
def read_parms(filename='parms.txt'):
    """Read parameters from parms.txt file"""
    parms = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                items=line.split(":")
                key=items[0].strip()
                if key == "mass":
                    parms['mass'] = float(items[1].split()[0])
                elif key == "blen":
                    parms['blen'] = float(items[1].split()[0])
                elif key == "spacing":
                    parms['spacing'] = float(items[1].split()[0])
    except FileNotFoundError:
        print(f"Warning: {filename} not found")
    return parms

def process_ern_data(ern_type):
    """Process fmap.ern data in memory for given ern_type"""
    data = []
    
    # Read from fmap.ern or fmap.ern.gz
    if os.path.exists('fmap.ern'):
        with open('fmap.ern', 'r') as f:
            for line in f:
                if ern_type in line:
                    parts = line.split()
                    if len(parts) >= 3:
                        data.append(float(parts[2]))  # column 3 (0-indexed as 2)
    elif os.path.exists('fmap.ern.gz'):
        with gzip.open('fmap.ern.gz', 'rt') as f:
            for line in f:
                if ern_type in line:
                    parts = line.split()
                    if len(parts) >= 3:
                        data.append(float(parts[2]))
    
    if not data:
        return None
        
    # Shuffle with seed 2019
    random.seed(2019)
    random.shuffle(data)
    
    return data

def compute_msf_stats(data):
    """Compute MSF statistics using FP method"""
    results = []
    if len(data)<100:
        block=2
    else:
        block=100
    for i in range(block, len(data)+block, block):
        if i <= len(data):
            subset = np.array(data[:i])
        else:
            subset = np.array(data)
        std = fp_stderr(subset)
        mean_val = np.mean(subset)
        std_val = np.std(subset)
        results.append((len(subset), mean_val, std_val, std))
    
    return results

def run_stats():
    """Main function to replicate stats.sh functionality"""
    # Read parameters
    parms = read_parms()
    if 'mass' not in parms or 'blen' not in parms:
        print("Error: mass and blen must be defined in parms.txt")
        return
    
    mass = parms['mass']
    blen = parms['blen']
    conv = 6.022140857e+23 / 1e+27 / (mass**2) * 1000 * 10000
    try:
        dx = parms['spacing']
    except:
        dx = 0.6
    
    results = {}
    
    # Process both ern types
    for ern_type in ['v+e', 'vol']:
        data = process_ern_data(ern_type)
        if data:
            msf_results = compute_msf_stats(data)
            if msf_results:
                # Write MSF file
                #msf_file = f"fmap.ern.{ern_type}.msf"
                #with open(msf_file, 'w') as f:
                #    for n, mean_val, std_val, stderr in msf_results:
                #        f.write(f"{n:16d}{mean_val:16.5f} {std_val:16.5f} {stderr:16.5f}\n")
                
                # Write B22 file
                b22_file = f"fmap.ern.{ern_type}.msf.b22"
                with open(b22_file, 'w') as f:
                    for n, mean_val, std_val, stderr in msf_results:
                        a = 0.5 * ((blen**3 - mean_val) * conv * dx**3)
                        b = 0.5 * (std_val * conv * dx**3)
                        c = 0.5 * (stderr * conv * dx**3)
                        f.write(f"{n} {a} {b} {c}\n")
                
                # Get last result for final calculations
                n, mean_val, std_val, stderr = msf_results[-1]
                
                # Calculate volumes and B22 values
                if ern_type == 'vol':
                    vco = (blen**3 - mean_val) * dx**3
                    results['FMAPvco'] = vco
                    results['FMAPexb22'] = 0.5 * ((blen**3 - mean_val) * conv * dx**3)
                    results['FMAPexb22std'] = 0.5 * (stderr * conv * dx**3)
                elif ern_type == 'v+e':
                    vs = (blen**3 - mean_val) * dx**3
                    results['FMAPvs'] = vs
                    results['FMAPb22'] = 0.5 * ((blen**3 - mean_val) * conv * dx**3)
                    results['FMAPb22std'] = 0.5 * (stderr * conv * dx**3)
    
    # Calculate vr
    if 'FMAPvco' in results:
        vr = (results['FMAPvco'] / (32.0/3 * math.pi))**(1.0/3.0) * 2.0
        results['FMAPvr'] = vr
    
    # Write results to parms.txt
    with open('parms.txt', 'a') as f:
        for key, value in results.items():
            if 'vr' in key:
                f.write(f"{key:<12} :{value:16.3f} angstrom\n")
            elif 'vco' in key or 'vs' in key:
                f.write(f"{key:<12} :{value:16.3f} angstrom^3\n")
            elif 'std' in key:
                f.write(f"{key:<12} :{value:16.4f} x 10^-4 mL mol / g^2\n")
            else:
                f.write(f"{key:<12} :{value:16.3f} x 10^-4 mL mol / g^2\n")
    
    #print("Stats processing complete!")

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    if len(argv) > 0 and argv[0] == "stats":
        run_stats()
    elif len(argv) >= 2:
        fname = argv[0]
        pos = int(argv[1])
        dat = getdata(fname, pos)
        npdat = np.array(dat)
        std = fp_stderr(npdat)
        print(npdat.size, np.mean(npdat), np.std(npdat), std)
    else:
        print("Usage:")
        print("  python script.py <filename> <column>  # Original FP functionality")
        print("  python script.py stats                # Run stats.sh functionality")

if __name__ == "__main__":
    main()
