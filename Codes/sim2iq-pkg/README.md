#Sim2Iq
Sim2Iq calculates the small-angle X-ray scattering intensity profile from atomic coordinates of protein molecules in a periodic box. It reimplements the core functionality of DENSS to get intensity profile from pdb, with major changes. It handles periodic boundary conditions and uses OpenMP for parallelization.

## Requirement
fftw3, openmp, cffi

## Install
from pip

    pip install sim2iq

from source

    make build
    make install

## Usage
### command line

    sim2iqpre -f input.pdb -o input.pqr
    sim2iq input_pqr ngrid box_length output_dat --steps 4 --traj traj.txt   

See tests/test.sh for detials

### library via cffi
See sim2iq.tools.sim2iq_calc for details