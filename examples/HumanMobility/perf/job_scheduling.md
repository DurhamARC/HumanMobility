## A SLURM job scheduling suite

### Scripts

Whenever you want to submit a job for a simulation (without or with profiling measurement) or for other applications (like visualisation), use the provided job submission script.

* `./submit.sh --${type} -N ${nodes} -n ${ntasks} -c ${threads} [-p ${partition}]` .. a unified job submission script to submit a job of choice with possible types to select from: `gen` .. to generate data of humans and a map of rivers; `vis` .. to plot simulation results; `map` .. to plot only geography without humans; a set of other options to run a simulation (`run` .. without profiling; `likwid` .. profiled with _LIKWID_; `maqao` .. profiled with _MAQAO_; `aps` .. profiled with _Intel APS_; `vtune` .. profiled with _Intel VTune_; `advisor` .. profiled with _Intel Advisor_; `inspector` .. profiled with _Intel Inspector_; `scorep` .. profiled with _Score-P_)

* `job.sh` .. a unified selector script that automatically detects partition and passes SLURM parameters to appropriate simulation scripts. Works on both COSMA5 and legacy COSMA partitions.

* `gen.sh` .. a script to generate data of humans and a map of rivers (no parameters are passed)

* `visualise.sh` .. a script to plot simulation results (no parameters are passed)

* `map.sh` .. a script to plot only geography without humans (no parameters are passed)

* `run.sh` .. a unified simulation script without profiling that automatically detects partition (COSMA5 vs COSMA) and execution mode (serial vs parallel) based on SLURM parameters

* `run_perf.sh` .. a unified performance analysis script that handles all profiling tools (`--likwid`, `--maqao`, `--aps`, `--vtune`, `--advisor`, `--inspector`, `--scorep`) with automatic partition detection and resource-based directory naming

#### Output Directory Structure

The unified scripts automatically generate descriptive directory names based on:
- Tool used (for performance analysis)
- Number of nodes, MPI tasks, and threads
- Partition (COSMA5 vs COSMA)

Examples:
- `data-rivers-1-1-64` (basic run on COSMA5: 1 node, 1 task, 64 threads)
- `data-rivers-1-1-16-cosma` (basic run on COSMA: 1 node, 1 task, 16 threads)
- `data-likwid-2-4-32` (LIKWID analysis on COSMA5: 2 nodes, 4 tasks, 32 threads each)
- `data-vtune-1-2-8-cosma` (VTune analysis on COSMA: 1 node, 2 tasks, 8 threads each)

#### Usage examples

Basic usage with defaults (1 node, 1 MPI rank, 64 threads on COSMA5):
```bash
./submit.sh --run                     # Basic simulation
./submit.sh --maqao                   # MAQAO performance analysis
./submit.sh --vtune                   # Intel VTune profiling
```

With custom resource allocation:
```bash
./submit.sh --likwid -N 1 -n 1 -c 32     # 1 node, 1 MPI rank, 32 threads
./submit.sh --advisor -N 1 -n 1 -c 16    # 1 node, 1 MPI rank, 16 threads
./submit.sh --inspector -N 1 -n 1 -c 8   # 1 node, 1 MPI rank, 8 threads
./submit.sh --scorep -N 1 -n 4 -c 8      # 1 node, 4 MPI ranks, 8 threads each
```

Partition-specific usage:
```bash
# COSMA5 (modern hardware, default)
./submit.sh --run -N 1 -n 1 -c 64        # Uses swift_intel2025 binary
./submit.sh --run -N 1 -n 1 -c 64 -p cosma5  # Explicit partition specification

# COSMA (legacy hardware)
./submit.sh --run -N 1 -n 1 -c 16 -p cosma   # Uses swift_cosma binary, max 16 cores/node
./submit.sh --vtune -N 2 -n 2 -c 8 -p cosma  # Distributed across legacy nodes
```

#### Key Features of the Unified Approach

* **Automatic partition detection**: Scripts automatically detect whether running on COSMA5 (modern) or COSMA (legacy) hardware and adapt accordingly
* **Dynamic directory naming**: Output directories are automatically named based on tool, resources, and partition (e.g., `data-likwid-1-2-32-cosma`)
* **Unified interface**: Same calling convention for all tools and partitions
* **Resource-aware configuration**: Automatic selection of appropriate SWIFT binaries and compiler versions based on partition

**Benefits of Unified Scripts:**
- Single maintenance point for each functionality
- Automatic adaptation to different hardware
- Consistent naming scheme across all tools
- Reduced complexity and improved maintainability
- Future-proof design for new partitions (cosma7, cosma8, etc.)

### Performance analysis

For strong scaling, I begin with a minimal simulation of human mobility on a square 100 x 100 humans on 10 x 10 km (for minimum of either 1000 steps or 10 minutes). For weak scaling, I'd like to begin on 2 ranks from a square 100 x 100 humans on 10 x 10 km and then on 8 and 18 ranks, increasing the problem size accordingly.

To conduct performance analysis, we have the following queues:

* The _cosma5_ queue is comprised of the new COSMA5 nodes, a total of 3 nodes each with 256 cores and 1.5TB RAM
* The _cosma_ queue is comprised of the old COSMA nodes, a total of ~160 nodes each with 16 cores and 126GB RAM

#### Strong scaling

**Performance Analysis for Strong Scaling**

For strong scaling analysis, we keep the problem size constant (100 x 100 humans on 10 x 10 km) and vary the number of computational resources:

**Base Configuration (Serial Baseline):**
```bash
./submit.sh --run -N 1 -n 1 -c 1     # 1 node, 1 MPI rank, 1 thread (serial)
```

**Intra-node Scaling (COSMA5 queue - threading analysis):**
```bash
# Test threading efficiency within a single node
./submit.sh --run -N 1 -n 1 -c 16    # 1 node, 1 MPI rank, 16 threads
./submit.sh --run -N 1 -n 1 -c 32    # 1 node, 1 MPI rank, 32 threads  
./submit.sh --run -N 1 -n 1 -c 64    # 1 node, 1 MPI rank, 64 threads
./submit.sh --run -N 1 -n 1 -c 128   # 1 node, 1 MPI rank, 128 threads
./submit.sh --run -N 1 -n 1 -c 256   # 1 node, 1 MPI rank, 256 threads (max COSMA5)
```

**Inter-node Scaling (COSMA5 queue - MPI analysis):**
```bash
# Test MPI scaling across multiple nodes with fixed threads per rank
./submit.sh --run -N 2 -n 2 -c 128   # 2 nodes, 2 MPI ranks, 128 threads each
./submit.sh --run -N 3 -n 3 -c 85    # 3 nodes, 3 MPI ranks, ~85 threads each
```

**COSMA queue comparison:**
```bash
# Compare with older COSMA hardware
./submit.sh --run -N 1 -n 1 -c 16 -p cosma   # 1 old COSMA node, 1 MPI rank, 16 threads (max)
./submit.sh --run -N 2 -n 2 -c 8 -p cosma    # 2 old COSMA nodes, 2 MPI ranks, 8 threads each
./submit.sh --run -N 4 -n 4 -c 4 -p cosma    # 4 old COSMA nodes, 4 MPI ranks, 4 threads each
```

**Performance Analysis Integration:**
```bash
# Run the same configurations with different profiling tools
./submit.sh --likwid -N 1 -n 1 -c 64    # LIKWID analysis on COSMA5
./submit.sh --maqao -N 1 -n 1 -c 16 -p cosma  # MAQAO analysis on COSMA
./submit.sh --vtune -N 2 -n 2 -c 128   # VTune analysis across 2 COSMA5 nodes
```

#### Weak scaling

**Performance Analysis for Weak Scaling**

For weak scaling, we maintain constant work per computational unit by increasing problem size proportionally with the number of nodes. For 2D problems like human mobility, the problem size should scale with the square root of the number of computing nodes.

**Scaling Strategy:**
- When we double the number of nodes, we should increase the linear dimension by ~1.414 (√2)
- For example: 100×100 humans → 141×141 humans → 200×200 humans (approximately)

**Baseline (2 nodes, 100 x 100 humans on 10 x 10 km):**
```bash
./submit.sh --run -N 2 -n 2 -c 16 -p cosma  # 2 nodes, 2 ranks (1 per node), 16 threads each
# or on COSMA5:
./submit.sh --run -N 2 -n 2 -c 128          # 2 nodes, 2 ranks, 128 threads each
```

**Scale to 8 nodes (200 x 200 humans on 20 x 20 km):**
```bash
./submit.sh --run -N 8 -n 8 -c 16 -p cosma  # 8 nodes, 8 ranks (1 per node), 16 threads each
# or on COSMA5:
./submit.sh --run -N 8 -n 8 -c 128          # 8 nodes, 8 ranks, 128 threads each
```

**Scale to 18 nodes (300 x 300 humans on 30 x 30 km):**
```bash
./submit.sh --run -N 18 -n 18 -c 16 -p cosma  # 18 nodes, 18 ranks (1 per node), 16 threads each
# or on COSMA5:
./submit.sh --run -N 18 -n 18 -c 128          # 18 nodes, 18 ranks, 128 threads each
```

**Important**: When submitting these jobs, you'll need to adjust your initial conditions file to match the problem size. Use the `makeIC.py` script with the appropriate `-n` parameter:

```bash
# For 2 nodes (100x100 humans)
python makeIC.py -n 100 -o humans-rivers-3.hdf5

# For 8 nodes (200x200 humans)
python makeIC.py -n 200 -o humans-rivers-3.hdf5

# For 18 nodes (300x300 humans)
python makeIC.py -n 300 -o humans-rivers-3.hdf5
```

**Profiling at Scale:**
```bash
# Profile weak scaling performance at different node counts
./submit.sh --likwid -N 2 -n 2 -c 16 -p cosma  # Baseline profiling (2 nodes)
./submit.sh --likwid -N 8 -n 8 -c 16 -p cosma  # 8-node scaling analysis
./submit.sh --vtune -N 18 -n 18 -c 16 -p cosma # 18-node analysis
```

**Expected Results:**
- Execution time should remain approximately constant across all configurations
- Memory usage per node should remain consistent
- Communication overhead should grow modestly with node count

_Note: For the most accurate weak scaling measurements, create separate initial condition files for each problem size and ensure the work per node remains constant across all configurations._
