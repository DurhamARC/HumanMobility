#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-rivers-3
RIVERS=river-rivers-3
HUMANMOBILITY=humanMobility
DATA=data-rivers-3
IMAGES=images-rivers-3

# Remove previously generated data and images to regenerate the new ones
rm ${HUMANS}.hdf5
rm ${RIVERS}.hdf5
rm ${DATA}/*
rm ${IMAGES}/*

# Generate acceleration field for river
if [ ! -e ${RIVERS}.hdf5 ]
then
    echo "Generating acceleration field for the river..."
    # python3 makeRandomRiver.py -b 1000000 -g 100000 -f ${RIVERS}.hdf5
    python3 makeRivers.py -b 100000 -g 10000 -f ${RIVERS}.hdf5
fi

# Generate the initial conditions if they are not present.
if [ ! -e ${HUMANS}.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    # python3 makeIC.py -t gas -f ${HUMANS}.hdf5  # -t particles
#     python3 makeIC.py -t gas -n 10000 -b 1000000 -f ${HUMANS}.hdf5  # -t particles
    python3 makeIC.py -t gas -n 1000 -b 100000 -f ${HUMANS}.hdf5  # -t particles
fi

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

# Run SWIFT with SLURM environment variables
# mpirun -np ${SLURM_NTASKS:-8} swift_mpi --threads=${SLURM_CPUS_PER_TASK:-16} \
# likwid-mpirun -np ${SLURM_NTASKS:-8} -t ${SLURM_CPUS_PER_TASK:-16} -omp intel -g MEM -- \

# likwid-mpirun -np ${SLURM_NTASKS:-8} -t ${SLURM_CPUS_PER_TASK:-16} -omp intel -g FLOPS_DP -- \
#     swift_mpi --threads=${SLURM_CPUS_PER_TASK:-16} \
#     -A -s -g -G \
#     --hm-river \
#     --hm-randomwalk \
#     -n 10000 \
#     ${HUMANMOBILITY}.yml

#OMP_NUM_THREADS=4 mpirun -n 4 bin/bt-mz.B.x
    # --number-processes=${SLURM_NTASKS:-8} \
    # --envv_OMP_NUM_THREADS=4 -- bin/bt-mz.B.x
#maqao oneview -R1 --mpi-command="mpirun -np ${SLURM_NTASKS:-8}" \

# module list

# gdb --args swift -g --threads=4 -n 10000 ${HUMANMOBILITY}.yml # -A -s
swift_intel2025 --threads=${SLURM_CPUS_PER_TASK:-16} \
    -A -s -g -G --hm-river --hm-randomwalk -n 10000 ${HUMANMOBILITY}.yml
# mpirun -np ${SLURM_NTASKS:-8} \
#     swift_mpi_intel2025 --threads=${SLURM_CPUS_PER_TASK:-16} \
#     -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
# likwid-perfctr -f -C 0 -g MEMREAD swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 ${HUMANMOBILITY}.yml # -A -s -g -G
# likwid-perfctr -f -C 0 -g MEM swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 10 ${HUMANMOBILITY}.yml # -A -s -g -G
# likwid-perfctr -f -C 0 -g FLOPS_DP swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 1000 ${HUMANMOBILITY}.yml # -A -s -g -G
# swift_intel2025 -h | grep version
# likwid-mpirun -np 8 -t 16 -omp intel -g MEM -- swift_mpi --threads=16 -A -s -g -G --hm-river --hm-randomwalk -n 10000 ${HUMANMOBILITY}.yml
# likwid-perfctr -a
