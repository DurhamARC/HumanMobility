#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-rivers-3
RIVERS=river-rivers-3
HUMANMOBILITY=humanMobility
DATA=data-rivers-3
IMAGES=images-rivers-3

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

# Run SWIFT with double-precision FLOPS measurement
likwid-perfctr -f -C 0-15 -g FLOPS_DP -- \
    swift_intel2025 --threads=${SLURM_CPUS_PER_TASK:-16} \
    -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
# likwid-mpirun -np ${SLURM_NTASKS:-1} -t ${SLURM_CPUS_PER_TASK:-16} -omp intel -g FLOPS_DP -- \
#     swift_mpi_intel2025 --threads=${SLURM_CPUS_PER_TASK:-16} \
#     -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
