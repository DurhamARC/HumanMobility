#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-3
RIVER=river-3
HUMANMOBILITY=humanMobility
DATA=data-3

# Remove previously generated data and images to regenerate the new ones
# rm ${HUMANS}.hdf5
# rm ${RIVER}.hdf5
# rm data/*
# rm images/*

# Generate the initial conditions if they are not present.
if [ ! -e ${HUMANS}.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    # python3 makeIC.py -t gas -f ${HUMANS}.hdf5  # -t particles
#     python3 makeIC.py -t gas -n 10000 -b 1000000 -f ${HUMANS}.hdf5  # -t particles
    python3 makeIC.py -t gas -n 100 -b 10000 -f ${HUMANS}.hdf5  # -t particles
fi

# Generate acceleration field for river
if [ ! -e ${RIVER}.hdf5 ]
then
    echo "Generating acceleration field for the river..."
    # python3 makeRandomRiver.py -b 1000000 -g 100000 -f ${RIVER}.hdf5
    python3 makeRandomRiver.py -b 10000 -g 1000 -f ${RIVER}.hdf5
fi

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVER
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

# Run SWIFT with SLURM environment variables
# mpirun -np ${SLURM_NTASKS:-8} swift_mpi --threads=${SLURM_CPUS_PER_TASK:-16} \
# likwid-mpirun -np ${SLURM_NTASKS:-8} -t ${SLURM_CPUS_PER_TASK:-16} -omp intel -g MEM -- \
likwid-mpirun -np ${SLURM_NTASKS:-8} -t ${SLURM_CPUS_PER_TASK:-16} -omp intel -g FLOPS_DP -- \
    swift_mpi --threads=${SLURM_CPUS_PER_TASK:-16} \
    -A -s -g -G \
    --hm-river \
    --hm-randomwalk \
    -n 1000 \
    ${HUMANMOBILITY}.yml

# module list

# gdb --args swift -g --threads=4 -n 10000 ${HUMANMOBILITY}.yml # -A -s
# swift --hm-river --hm-randomwalk --threads=8 -n 1000 ${HUMANMOBILITY}.yml # -A -s -g -G 
# likwid-perfctr -f -C 0 -g MEMREAD swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 ${HUMANMOBILITY}.yml # -A -s -g -G
# likwid-perfctr -f -C 0 -g MEM swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 10 ${HUMANMOBILITY}.yml # -A -s -g -G
# likwid-perfctr -f -C 0 -g FLOPS_DP swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 1000 ${HUMANMOBILITY}.yml # -A -s -g -G
# mpirun -np 8 swift_mpi -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 1000 ${HUMANMOBILITY}.yml
# swift_intel2025 -h | grep version
# likwid-mpirun -np 8 -t 16 -omp intel -g MEM -- swift_mpi --threads=16 -A -s -g -G --hm-river --hm-randomwalk -n 10000 ${HUMANMOBILITY}.yml
# likwid-perfctr -a
