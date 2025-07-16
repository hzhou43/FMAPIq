# FMAPB2

FMAPB2 implements FMAP, a method based on fast Fourier transform (FFT), to calculate second virial coefficients (B2) for proteins represented at the all-atom level in implicit solvent. FMAP stands for FFT-based Modeling of Atomistic Protein-protein interactions. In FMAPB2, we express terms of the protein-protein inteaction energy as correlation functions, and evaluate them by FFT. These terms include steric repulsion, nonpolar attraction (in the form of a Lennard-Jones potential), and electrostatic interactions (in the form of a Debye-Hückel potential).
The input to FMAPB2 is the structure of the protein, in PQR format. The user also sets the solvent condition, including ionic strength and temperature.

## Requirement
fftw3, openmp, cffi

## Install
from pip

    pip install fmapb2

from source

    make build
    make install

## Usage
### command line

    fmapb2pre A.pqr ion tempC escl vscl [B.pqr]
    fmapb2    subA.vdw parms.txt ang.dat [subB.vdw]

See tests/test.sh for detials

### library via cffi
See fmapb2.tools.fmap for details