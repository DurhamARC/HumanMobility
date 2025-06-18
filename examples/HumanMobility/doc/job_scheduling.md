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
./submit.sh --run -N 3 -n 3 -c 85    # 3 nodes, 3 MPI ranks, ~85 threads each (max COSMA5)
```

**COSMA queue comparison:**
```bash
# Compare with older COSMA hardware (can use many more nodes)
./submit.sh --run -N 1 -n 1 -c 16 -p cosma   # 1 old COSMA node, 1 MPI rank, 16 threads (max)
./submit.sh --run -N 2 -n 2 -c 8 -p cosma    # 2 old COSMA nodes, 2 MPI ranks, 8 threads each
./submit.sh --run -N 4 -n 4 -c 4 -p cosma    # 4 old COSMA nodes, 4 MPI ranks, 4 threads each
./submit.sh --run -N 8 -n 8 -c 2 -p cosma    # 8 old COSMA nodes, 8 MPI ranks, 2 threads each
```

**Performance Analysis Integration:**
```bash
# Run the same configurations with different profiling tools
./submit.sh --likwid -N 1 -n 1 -c 64    # LIKWID analysis on COSMA5
./submit.sh --maqao -N 1 -n 1 -c 16 -p cosma  # MAQAO analysis on COSMA
./submit.sh --vtune -N 2 -n 2 -c 128   # VTune analysis across 2 COSMA5 nodes
./submit.sh --vtune -N 3 -n 3 -c 85    # VTune analysis across 3 COSMA5 nodes (maximum)
```

#### Weak scaling

**Performance Analysis for Weak Scaling**

For weak scaling, we maintain constant work per computational unit by increasing problem size proportionally with the number of nodes. For 2D problems like human mobility, the problem size should scale with the square root of the number of computing nodes.

**Scaling Strategy:**
- When we double the number of nodes, we should increase the linear dimension by ~1.414 (√2)
- For example: 100×100 humans → 141×141 humans → 173×173 humans (approximately)
- **Important**: Each scaling configuration requires its own input files with appropriately sized domains and human populations

**Step 1: Generate Input Files for Each Configuration**

Before running weak scaling tests, generate the appropriate input files for each node count:

```bash
# For 2 nodes baseline (100x100 humans on 10x10 km)
./submit.sh --gen -- -n 100 -b 10000
# Creates: humans-rivers-100.hdf5, river-rivers-100.hdf5 (or with -cosma suffix)

# For 3 nodes (122x122 humans on 12.2x12.2 km) - COSMA5 maximum  
./submit.sh --gen -- -n 122 -b 12200
# Creates: humans-rivers-122.hdf5, river-rivers-122.hdf5

# For extended scaling on COSMA legacy hardware
# For 8 nodes (200x200 humans on 20x20 km)
./submit.sh --gen -- -n 200 -b 20000
# Creates: humans-rivers-200.hdf5, river-rivers-200.hdf5

# For 18 nodes (300x300 humans on 30x30 km)
./submit.sh --gen -- -n 300 -b 30000
# Creates: humans-rivers-300.hdf5, river-rivers-300.hdf5
```

**Step 2: Run Weak Scaling Tests**

Now run the simulations using the corresponding input files for each configuration:

**Baseline (2 nodes, 100 x 100 humans on 10 x 10 km):**
```bash
./submit.sh --run -N 2 -n 2 -c 16 -p cosma -- --num-humans 100  # Uses humans-rivers-100.hdf5
# or on COSMA5:
./submit.sh --run -N 2 -n 2 -c 128 -- --num-humans 100         # Uses humans-rivers-100.hdf5
```

**Scale to 3 nodes (122 x 122 humans on 12.2 x 12.2 km) - COSMA5 maximum:**
```bash
./submit.sh --run -N 3 -n 3 -c 16 -p cosma -- --num-humans 122  # Uses humans-rivers-122.hdf5
# or on COSMA5:
./submit.sh --run -N 3 -n 3 -c 64 -- --num-humans 122          # Uses humans-rivers-122.hdf5
```

**Extended scaling on COSMA (legacy hardware with more nodes available):**
```bash
./submit.sh --run -N 8 -n 8 -c 16 -p cosma -- --num-humans 200   # Uses humans-rivers-200.hdf5
./submit.sh --run -N 18 -n 18 -c 16 -p cosma -- --num-humans 300 # Uses humans-rivers-300.hdf5
```

**Step 3: Profile at Scale**

Run performance analysis with the appropriate input files:

```bash
# Profile weak scaling performance at different node counts
./submit.sh --likwid -N 2 -n 2 -c 64 -- --num-humans 100        # Baseline profiling (2 COSMA5 nodes)
./submit.sh --likwid -N 3 -n 3 -c 64 -- --num-humans 122         # 3-node scaling analysis (max COSMA5)
./submit.sh --vtune -N 8 -n 8 -c 16 -p cosma -- --num-humans 200 # 8-node analysis on COSMA legacy
./submit.sh --vtune -N 18 -n 18 -c 16 -p cosma -- --num-humans 300 # 18-node analysis on COSMA legacy
```

**Input File Naming Convention:**

The [`gen.sh`](examples/HumanMobility/gen.sh) script generates files with names based on the `--num-humans` parameter:
- `humans-rivers-{NUM_HUMANS}.hdf5` (or `humans-rivers-{NUM_HUMANS}-cosma.hdf5` on COSMA partition)
- `river-rivers-{NUM_HUMANS}.hdf5` (or `river-rivers-{NUM_HUMANS}-cosma.hdf5` on COSMA partition)

The [`run.sh`](examples/HumanMobility/run.sh) script automatically uses the correct input files when passed the `--num-humans` parameter, ensuring that:
- 2-node runs use the 100×100 human configuration
- 3-node runs use the 122×122 human configuration  
- 8-node runs use the 200×200 human configuration
- 18-node runs use the 300×300 human configuration

**Expected Results:**
- Execution time should remain approximately constant across all configurations
- Memory usage per node should remain consistent (proportional to local problem size)
- Communication overhead should grow modestly with node count
- Each configuration maintains the same work per computational unit

**Key Points for Weak Scaling:**
1. **Generate all input files first** before running scaling tests
2. **Use `--num-humans` parameter** to specify which input files to use
3. **Ensure domain size scales with human count** (maintain constant density)
4. **Verify file naming consistency** across partitions (COSMA vs COSMA5)

_Note: COSMA5 is limited to 3 nodes maximum, so extended weak scaling studies should use the legacy COSMA partition which has ~160 nodes available. The input file generation step ensures that each node configuration has appropriately sized problems while maintaining constant work density._

### Visualisation of performance results

The performance analysis results can be visualized using the `pa_vis.py` script, which creates comprehensive PDF reports showing memory usage, CPU performance, and thread utilization across MPI ranks.

#### Usage

```bash
python pa_vis.py -d <output_directory>
```

**Command-line options:**
- `-d, --directory`: Directory containing performance log files (default: current directory)
- `--use-balance-logs`: Force use of `rank_*_balance.log` files
- `--use-dat-files`: Force use of `*_report-rank*-step*.dat` files  
- `--use-thread-files`: Force use of `thread_info_MPI-step*.dat` files

The script automatically detects and processes the following data sources (in order of preference):
1. **Balance logs**: `rank_memory_balance.log` and `rank_cpu_balance.log`
2. **Individual report files**: `memuse_report-rank*-step*.dat` and `mpiuse_report-rank*-step*.dat`
3. **Thread timing data**: `thread_info_MPI-step*.dat` or `thread_info-step*.dat`

#### Output

The tool generates a two-page PDF report:

**Page 1: System-level Performance**
- Memory usage across MPI ranks (average and peak)
- CPU time distribution and parallel efficiency
- Load balancing metrics

**Page 2: Thread-level Analysis**
- Individual thread utilization by rank
- Thread balance ratios within each rank
- Thread efficiency visualization

#### Example Usage

For a simulation run with directory `data-rivers-1-4-4-cosma`:

```bash
python pa_vis.py -d data-rivers-1-4-4-cosma
```

**Sample Output:**
```
Detected partition: cosma
Searching for performance files in: .../SWIFT/examples/HumanMobility/data-rivers-1-4-4-cosma
Detected 1 nodes from directory name
Output will be saved as: pa_vis-1-4-4-cosma.pdf
Using balance log files: .../rank_memory_balance.log, .../rank_cpu_balance.log
Using MPI thread info files: 100 files found
Found 4 ranks and 17 balance steps

============================================================
PERFORMANCE ANALYSIS SUMMARY
============================================================
Simulation: 17 steps, 4 MPI ranks
Memory & CPU steps analyzed: 17, 25, 33, 41, 49, 57, 65, 75, 90, 126, 190, 252, 370, 524, 635, 745, 934
Thread steps analyzed: 2 to 992 (100 total steps)

MEMORY USAGE:
  Rank 0: Avg = 1196.63 GB, Max = 1242.03 GB
  Rank 1: Avg = 1167.23 GB, Max = 1201.75 GB
  Rank 2: Avg = 1165.93 GB, Max = 1233.33 GB
  Rank 3: Avg = 1151.01 GB, Max = 1208.53 GB
  Overall average: 1170.20 GB
  Memory balance: 1.4%

CPU USAGE:
  Rank 0: Total = 81.9 seconds
  Rank 1: Total = 81.0 seconds
  Rank 2: Total = 81.0 seconds
  Rank 3: Total = 79.8 seconds
  Overall total: 323.8 seconds
  Parallel efficiency: 98.8%
  Load balance ratio: 0.97

INDIVIDUAL THREAD USAGE:
  Rank 0-3 Thread 0-3: 97.7 - 246.8 s (per thread)
  Total threads active: 16
  Overall total thread time: 2114.4 s
  Thread balance ratio: 0.40

Output saved to: pa_vis-1-4-4-cosma.pdf
============================================================
```

The output PDF filename automatically includes the run configuration (e.g., `pa_vis-1-4-4-cosma.pdf`) for easy identification and comparison across different configurations.