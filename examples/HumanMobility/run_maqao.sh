#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-maqao
RIVERS=river-maqao
HUMANMOBILITY=humanMobility
DATA=data-maqao
IMAGES=images-maqao

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

#OMP_NUM_THREADS=4 mpirun -n 4 bin/bt-mz.B.x
    # --number-processes=${SLURM_NTASKS:-8} \
    # --envv_OMP_NUM_THREADS=4 -- bin/bt-mz.B.x
#maqao oneview -R1 --mpi-command="mpirun -np ${SLURM_NTASKS:-8}" \

# Run Maqao oneview analysis on SWIFT (single node, threaded)
# maqao oneview -R1 --output-format=all -- \
#     ${SWIFT} -A -s -g -G --hm-river --hm-randomwalk --threads=${SLURM_CPUS_PER_TASK:-16} -n 1000 ${HUMANMOBILITY}.yml

# Run Maqao oneview analysis on SWIFT MPI version with correct options
maqao oneview -R1 --output-format=all \
    --mpi-command="mpirun -np ${SLURM_NTASKS:-4}" \
    --envv_OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" -- \
    ${SWIFT_MPI} -A -s -g -G --hm-river --hm-randomwalk --threads=${SLURM_CPUS_PER_TASK:-4} -n 1000 ${HUMANMOBILITY}.yml