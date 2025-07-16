# MAEppi
MAEPPI (many-protein atomistic energy from pre-computed pair interaction) enables simulations of atomistic proteins at the speed of Lennard-Jones particles.

## Requirement
openmp, cffi

## Install
from pip

    pip install maeppi

from source

    make build
    make install

## Usage
### command line

    maeppi_shm_server start parms.txt --load-only
    maeppi_sim <n_particles> <box_length> <kbt> <n_steps> <seed> <start.dat> <stop.dat> <sample_freq> <keyp> [tx] [rx] [zx]
    maeppi_shm_server cleanup parms.txt

See tests/test.sh for details

### library via cffi
See maeppi.tools.maeppi_sim, maeppi.tools.maeppi_ern for details