# A SLURM job scheduling suite

### Scripts

Whenever you want to submit a job for a simulation (without or with profiling measurement) or for other applications (like visualisation), use the provided job submission script.

* `./submit.sh --${type} -N ${nodes} -n ${ntasks} -c ${threads}` .. a job submission script to submit a job of choice with possible types to select from: `gen` .. to generate data of humans and a map of rivers; `vis` .. to plot simulation results; `map` .. to plot only geography without humans; a set of other options to run a simulation (`run` .. without proifiling; `likwid` .. profiled with _likwid_; `maqao` .. profiled with _Maqao_; `vtune` .. profiled with _Intel VTune_; `advisor` .. profiled with _Intel Advisor_; `inspector` .. profiled with _Intel Inspector_; `scorep` .. profiled with _Score-P_)

* `job.sh` .. a selector script to pass SLURM parameters further down to specific simulation scripts

* `gen.sh` .. a script to generate data of humans and a map of rivers (no parameters are passed)

* `visualise.sh` .. a script to plot simulation results (no parameters are passed)

* `map.sh` .. a script to plot only geography without humans (no parameters are passed)

* `run.sh` .. a simulation script without profiling (with automatic selection of call mode between serial and parallel determined from additional parameters passed in `submit.sh`)

* `run_likwid.sh` .. a simulation script with profiling by means of _likwid_

* `run_maqao.sh` .. a simulation script with profiling by means of _Maqao_

* `run_vtune.sh` .. a simulation script with profiling by means of _Intel VTune_

* `run_advisor.sh` .. a simulation script with profiling by means of _Intel Advisor_

* `run_inspector.sh` .. a simulation script with profiling by means of _Intel Inspector_

* `run_scorep.sh` .. a simulation script with profiling by means of _Score-P_

### Usage examples

Basic usage with defaults (1 node, 1 MPI rank, 64 threads):
```bash
./submit.sh --run
./submit.sh --maqao
./submit.sh --vtune
```

With custom resource allocation:
```bash
./submit.sh --likwid -N 1 -n 1 -c 32     # 1 node, 1 MPI rank, 32 threads
./submit.sh --advisor -N 1 -n 1 -c 16    # 1 node, 1 MPI rank, 16 threads
./submit.sh --inspector -N 1 -n 1 -c 8   # 1 node, 1 MPI rank, 8 threads
./submit.sh --scorep -N 1 -n 4 -c 8      # 1 node, 4 MPI ranks, 8 threads
```

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

#### Weak scaling

**Performance Analysis for Weak Scaling**

For weak scaling, we maintain constant work per computational unit by increasing problem size proportionally with resources:

**Baseline (2 MPI ranks, 100 x 100 humans):**
```bash
./submit.sh --run -N 1 -n 2 -c 8     # Base: 100x100 humans, 2 ranks, 8 threads each
```

**Scale to 8 MPI ranks (200 x 200 humans):**
```bash
# COSMA5 queue (high core count nodes)
./submit.sh --run -N 1 -n 8 -c 32    # 1 node, 8 ranks, 32 threads each (256 total cores)

# COSMA queue (distributed across multiple nodes)
./submit.sh --run -N 2 -n 8 -c 2 -p cosma    # 2 nodes, 8 ranks, 2 threads each
./submit.sh --run -N 4 -n 8 -c 2 -p cosma    # 4 nodes, 8 ranks, 2 threads each  
./submit.sh --run -N 8 -n 8 -c 2 -p cosma    # 8 nodes, 8 ranks, 2 threads each
```

**Scale to 16 MPI ranks (300 x 300 humans):**
```bash
# COSMA5 queue (maximum utilization)
./submit.sh --run -N 1 -n 18 -c 14   # 1 node, 18 ranks, ~14 threads each

# COSMA queue (distributed)
./submit.sh --run -N 2 -n 16 -c 1 -p cosma   # 2 nodes, 16 ranks, 1 thread each
./submit.sh --run -N 9 -n 16 -c 1 -p cosma   # 9 nodes, 16 ranks, 1 thread each
./submit.sh --run -N 18 -n 16 -c 1 -p cosma  # 18 nodes, 16 ranks, 1 thread each
```

_Note: Ideally, 18 MPI ranks are needed for an accurate weak scaling estimate. However, the older version of Cosma is limited to a maximum of 16 cores per node._
