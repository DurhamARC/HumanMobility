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
# swift_intel2025 --threads=${SLURM_CPUS_PER_TASK:-16} \
#     -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
# mpirun -np ${SLURM_NTASKS:-8} \
#     swift_mpi_intel2025 --threads=${SLURM_CPUS_PER_TASK:-16} \
#     -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
# likwid-perfctr -f -C 0 -g MEMREAD swift -A -s -g -G --hm-river --hm-randomwalk --threads=${SLURM_CPUS_PER_TASK:-16} -n 100 ${HUMANMOBILITY}.yml
likwid-perfctr -f -C 0 -g MEM swift_intel2025 -A -s -g -G --hm-river --hm-randomwalk --threads=${SLURM_CPUS_PER_TASK:-16} -n 1000 ${HUMANMOBILITY}.yml
# likwid-perfctr -f -C 0 -g FLOPS_DP swift -A -s -g -G --hm-river --hm-randomwalk --threads=${SLURM_CPUS_PER_TASK:-16} -n 1000 ${HUMANMOBILITY}.yml
# swift_intel2025 -h | grep version
# likwid-mpirun -np ${SLURM_NTASKS:-1} -t ${SLURM_CPUS_PER_TASK:-16} -omp intel -g MEM -- \ # FLOPS_DP
    # swift_mpi_intel2025 --threads=${SLURM_CPUS_PER_TASK:-16} \
    # -A -s -g -G --hm-river --hm-randomwalk -n 10000 ${HUMANMOBILITY}.yml
# likwid-perfctr -a
