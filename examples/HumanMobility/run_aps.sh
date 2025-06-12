#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-aps
RIVERS=river-aps
HUMANMOBILITY=humanMobility
DATA=data-aps
IMAGES=images-aps

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

# Use SLURM environment variables for configuration
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

echo "Intel APS Analysis Configuration:"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"

# Enable SWIFT's built-in logging
# export SWIFT_TASK_DUMPS=1
# export SWIFT_MPIUSE_REPORTS=1
# export SWIFT_MEMUSE_REPORTS=1

# Run APS with automatic serial/parallel detection
if [ $ntasks -eq 1 ]; then
    echo "Running Intel APS on SERIAL SWIFT"
    
    echo "Running Intel APS performance snapshot analysis..."
    aps -r aps_serial_intranode \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml #--task-dumps=1

    echo "Generating APS summary report..."
    aps --report=summary --result-dir=aps_serial_intranode > aps_serial_summary.txt

    echo "Generating APS detailed report..."
    aps --report=detailed --result-dir=aps_serial_intranode > aps_serial_detailed.txt
else
    echo "Running Intel APS on PARALLEL SWIFT"
    
    echo "Running Intel APS MPI performance snapshot analysis..."
    srun -n $ntasks -c $cpus_per_task \
        aps -r aps_mpi_intranode \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml #--task-dumps=1

    echo "Generating APS MPI summary report..."
    aps --report=summary --result-dir=aps_mpi_intranode > aps_mpi_summary.txt

    echo "Generating APS MPI detailed report..."
    aps --report=detailed --result-dir=aps_mpi_intranode > aps_mpi_detailed.txt
fi

echo "Intel APS analysis complete. Results in aps_*_summary.txt and aps_*_detailed.txt files"
echo "Result directories: aps_*_intranode/"
echo ""
echo "Key APS metrics to check:"
echo "  - CPU utilization and efficiency"
echo "  - Memory bandwidth utilization"
echo "  - Elapsed time and performance characteristics"
echo "  - MPI communication overhead (for parallel runs)"
echo ""
echo "To view results:"
echo "  cat aps_*_summary.txt     # Quick overview"
echo "  cat aps_*_detailed.txt    # Detailed analysis"