#!/bin/bash

set -e

ELIBS=${EPIPHANY_HOME}/tools/host/lib:${LD_LIBRARY_PATH}
EHDF=${EPIPHANY_HDF}

#sudo -E LD_LIBRARY_PATH=${ELIBS} EPIPHANY_HDF=${EHDF} ./bin/main

sudo -E LD_LIBRARY_PATH=${ELIBS} EPIPHANY_HDF=${EHDF} mpirun -np 4 --hostfile ~/.mpi_hostfile ./bin/main


