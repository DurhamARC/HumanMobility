# Performance analysis

Application: Human Mobility

Parameters:

- Number of people: 10k → 1M humans
- Area of territory: 100 sq km → 10000 sq km

## Setup

* Cosma

## Benchmarks

### Memory bandwidth (MByte/s) on a single core
```
likwid-bench -t triad -W S0:2GB:1
```
```
====
Starting job ... 
...
LIKWID MICRO BENCHMARK
Test: triad
...
Using 1 work groups
Using 1 threads
...
Cycles:                 3591308877
CPU Clock:              2599987191
Cycle Clock:            2599987191
Time:                   1.381279e+00 sec
Iterations:             10
Iterations per thread:  10
Inner loop executions:  15625000
Size (Byte):            2000000000
Size per thread:        2000000000
Number of Flops:        1250000000
MFlops/s:               904.96
Data volume (Byte):     20000000000
MByte/s:                14479.33
Cycles per update:      5.746094
Cycles per cacheline:   45.968754
Loads per update:       3
Stores per update:      1
Load bytes per element: 24
Store bytes per elem.:  8
Load/store ratio:       3.00
Instructions:           2968750016
UOPs:                   4687500000
--------------------------------------------------------------------------------
```

### Memory bandwidth (MByte/s) on a full node
```
likwid-bench -t triad -W N:2GB:16
```
```
====
Starting job ...
...
LIKWID MICRO BENCHMARK
Test: triad
...
Using 1 work groups
Using 16 threads
...
Cycles:                 25008707817
CPU Clock:              2599930277
Cycle Clock:            2599930277
Time:                   9.618992e+00 sec
Iterations:             2048
Iterations per thread:  128
Inner loop executions:  976562
Size (Byte):            1999998976
Size per thread:        124999936
Number of Flops:        15999991808
MFlops/s:               1663.38
Data volume (Byte):     255999868928
MByte/s:                26614.00
Cycles per update:      3.126090
Cycles per cacheline:   25.008721
Loads per update:       3
Stores per update:      1
Load bytes per element: 24
Store bytes per elem.:  8
Load/store ratio:       3.00
Instructions:           37999980560
UOPs:                   59999969280
--------------------------------------------------------------------------------
```

### Peak floating-point performance (MFLOPS) on a single core
```
likwid-bench -t peakflops -w S0:10000MB:1
```
```
====
Starting job ...
...
LIKWID MICRO BENCHMARK
Test: peakflops
...
Using 1 work groups
Using 1 threads
...
Cycles:                 105343024434
CPU Clock:              2599963429
Cycle Clock:            2599963429
Time:                   4.051712e+01 sec
Iterations:             10
Iterations per thread:  10
Inner loop executions:  1250000000
Size (Byte):            10000000000
Size per thread:        10000000000
Number of Flops:        200000000000
MFlops/s:               4936.19
Data volume (Byte):     100000000000
MByte/s:                2468.09
Cycles per update:      8.427442
Cycles per cacheline:   67.419536
Loads per update:       1
Stores per update:      0
Load bytes per element: 8
Store bytes per elem.:  0
Instructions:           250000000032
UOPs:                   237500000000
--------------------------------------------------------------------------------
```

### Peak floating-point performance (MFLOPS) on a full node
```
likwid-bench -t peakflops -w N:10000MB:16
```
```
====
Starting job ...
...
LIKWID MICRO BENCHMARK
Test: peakflops
...
Using 1 work groups
Using 16 threads
...
Cycles:                 11756247712
CPU Clock:              2599984198
Cycle Clock:            2599984198
Time:                   4.521661e+00 sec
Iterations:             160
Iterations per thread:  10
Inner loop executions:  78125000
Size (Byte):            10000000000
Size per thread:        625000000
Number of Flops:        200000000000
MFlops/s:               44231.53
Data volume (Byte):     100000000000
MByte/s:                22115.77
Cycles per update:      0.940500
Cycles per cacheline:   7.523999
Loads per update:       1
Stores per update:      0
Load bytes per element: 8
Store bytes per elem.:  0
Instructions:           250000000032
UOPs:                   237500000000
--------------------------------------------------------------------------------
```

## Core (caches)

### Serial

#### Memory performance (`--enable-debug`)

##### 100 x 100 humans on 10 x 10 km for minimum(1000 steps, 10 minutes):
```
likwid-perfctr -f -C 0 -g MEM swift -A -s -g -G --hm-river --hm-randomwalk --threads=1 -n 1000 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: MEM
+-----------------------+----------+--------------+
|         Event         |  Counter |  HWThread 0  |
+-----------------------+----------+--------------+
|   INSTR_RETIRED_ANY   |   FIXC0  | 123986843446 |
| CPU_CLK_UNHALTED_CORE |   FIXC1  |  70409818611 |
|  CPU_CLK_UNHALTED_REF |   FIXC2  |  37073751360 |
|     TOPDOWN_SLOTS     |   FIXC3  | 422458260954 |
|      CAS_COUNT_RD     |  MBOX0C0 |      2104365 |
|      CAS_COUNT_WR     |  MBOX0C1 |      1941294 |
|      CAS_COUNT_RD     |  MBOX1C0 |      2034649 |
|      CAS_COUNT_WR     |  MBOX1C1 |      1885696 |
|      CAS_COUNT_RD     |  MBOX2C0 |      1971302 |
|      CAS_COUNT_WR     |  MBOX2C1 |      1827391 |
|      CAS_COUNT_RD     |  MBOX3C0 |      2004834 |
|      CAS_COUNT_WR     |  MBOX3C1 |      1856898 |
|      CAS_COUNT_RD     |  MBOX4C0 |      2103657 |
|      CAS_COUNT_WR     |  MBOX4C1 |      1938628 |
|      CAS_COUNT_RD     |  MBOX5C0 |      2029864 |
|      CAS_COUNT_WR     |  MBOX5C1 |      1879532 |
|      CAS_COUNT_RD     |  MBOX6C0 |      1974026 |
|      CAS_COUNT_WR     |  MBOX6C1 |      1830926 |
|      CAS_COUNT_RD     |  MBOX7C0 |      2019975 |
|      CAS_COUNT_WR     |  MBOX7C1 |      1867911 |
+-----------------------+----------+--------------+

+-----------------------------------+------------+
|               Metric              | HWThread 0 |
+-----------------------------------+------------+
|        Runtime (RDTSC) [s]        |    19.4184 |
|        Runtime unhalted [s]       |    35.2053 |
|            Clock [MHz]            |  3798.3176 |
|                CPI                |     0.5679 |
|  Memory read bandwidth [MBytes/s] |    53.5334 |
|  Memory read data volume [GBytes] |     1.0395 |
| Memory write bandwidth [MBytes/s] |    49.5309 |
| Memory write data volume [GBytes] |     0.9618 |
|    Memory bandwidth [MBytes/s]    |   103.0643 |
|    Memory data volume [GBytes]    |     2.0013 |
+-----------------------------------+------------+
```

##### 1000 x 1000 humans on 100 x 100 km for minimum(100 steps, 10 minutes):
```
likwid-perfctr -f -C 0 -g MEM swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: MEM
+-----------------------+----------+----------------+
|         Event         |  Counter |   HWThread 0   |
+-----------------------+----------+----------------+
|   INSTR_RETIRED_ANY   |   FIXC0  | 10078700900771 |
| CPU_CLK_UNHALTED_CORE |   FIXC1  |  4651632179068 |
|  CPU_CLK_UNHALTED_REF |   FIXC2  |  2450068924800 |
|     TOPDOWN_SLOTS     |   FIXC3  | 27909716623362 |
|      CAS_COUNT_RD     |  MBOX0C0 |      296253668 |
|      CAS_COUNT_WR     |  MBOX0C1 |      198353712 |
|      CAS_COUNT_RD     |  MBOX1C0 |      295934061 |
|      CAS_COUNT_WR     |  MBOX1C1 |      198108469 |
|      CAS_COUNT_RD     |  MBOX2C0 |      296181293 |
|      CAS_COUNT_WR     |  MBOX2C1 |      198068319 |
|      CAS_COUNT_RD     |  MBOX3C0 |      296381234 |
|      CAS_COUNT_WR     |  MBOX3C1 |      198231128 |
|      CAS_COUNT_RD     |  MBOX4C0 |      296422014 |
|      CAS_COUNT_WR     |  MBOX4C1 |      198429254 |
|      CAS_COUNT_RD     |  MBOX5C0 |      296071633 |
|      CAS_COUNT_WR     |  MBOX5C1 |      198181758 |
|      CAS_COUNT_RD     |  MBOX6C0 |      296372650 |
|      CAS_COUNT_WR     |  MBOX6C1 |      198253434 |
|      CAS_COUNT_RD     |  MBOX7C0 |      296656051 |
|      CAS_COUNT_WR     |  MBOX7C1 |      198468244 |
+-----------------------+----------+----------------+

+-----------------------------------+------------+
|               Metric              | HWThread 0 |
+-----------------------------------+------------+
|        Runtime (RDTSC) [s]        |  1229.3598 |
|        Runtime unhalted [s]       |  2325.8498 |
|            Clock [MHz]            |  3797.0889 |
|                CPI                |     0.4615 |
|  Memory read bandwidth [MBytes/s] |   123.3955 |
|  Memory read data volume [GBytes] |   151.6974 |
| Memory write bandwidth [MBytes/s] |    82.5715 |
| Memory write data volume [GBytes] |   101.5100 |
|    Memory bandwidth [MBytes/s]    |   205.9670 |
|    Memory data volume [GBytes]    |   253.2075 |
+-----------------------------------+------------+
```

##### 1000 x 1000 humans on 100 x 100 km for minimum, multiple rivers (1000 steps, 10 minutes):
```
likwid-perfctr -f -C 0 -g MEM swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: MEM
+-----------------------+---------+----------------+
|         Event         | Counter |   HWThread 0   |
+-----------------------+---------+----------------+
|   INSTR_RETIRED_ANY   |  FIXC0  | 16237189087458 |
| CPU_CLK_UNHALTED_CORE |  FIXC1  | 13024552770775 |
|  CPU_CLK_UNHALTED_REF |  FIXC2  | 10265924898034 |
|      CAS_COUNT_RD     | MBOX0C0 |     4634215368 |
|      CAS_COUNT_WR     | MBOX0C1 |     1900293496 |
|      CAS_COUNT_RD     | MBOX1C0 |     3896604436 |
|      CAS_COUNT_WR     | MBOX1C1 |     1068569638 |
|      CAS_COUNT_RD     | MBOX2C0 |     3801181672 |
|      CAS_COUNT_WR     | MBOX2C1 |     1139534217 |
|      CAS_COUNT_RD     | MBOX3C0 |     3903990616 |
|      CAS_COUNT_WR     | MBOX3C1 |     1066476734 |
+-----------------------+---------+----------------+

+-----------------------------------+------------+
|               Metric              | HWThread 0 |
+-----------------------------------+------------+
|        Runtime (RDTSC) [s]        |  3994.6320 |
|        Runtime unhalted [s]       |  5009.5327 |
|            Clock [MHz]            |  3298.6052 |
|                CPI                |     0.8021 |
|  Memory read bandwidth [MBytes/s] |   260.1250 |
|  Memory read data volume [GBytes] |  1039.1035 |
| Memory write bandwidth [MBytes/s] |    82.9092 |
| Memory write data volume [GBytes] |   331.1919 |
|    Memory bandwidth [MBytes/s]    |   343.0342 |
|    Memory data volume [GBytes]    |  1370.2954 |
+-----------------------------------+------------+
```

#### Floating-point performance (`--enable-debug`)

##### 100 x 100 humans on 10 x 10 km for minimum(1000 steps, 10 minutes):
```
likwid-perfctr -f -C 0 -g FLOPS_DP swift -A -s -g -G --hm-river --hm-randomwalk --threads=1 -n 1000 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: FLOPS_DP
+------------------------------------------+---------+--------------+
|                   Event                  | Counter |  HWThread 0  |
+------------------------------------------+---------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 143861628268 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  79237524584 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  41722922320 |
|               TOPDOWN_SLOTS              |  FIXC3  | 475423666740 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |    418829623 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |   1872322584 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |    295344149 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |            0 |
+------------------------------------------+---------+--------------+

+----------------------+------------+
|        Metric        | HWThread 0 |
+----------------------+------------+
|  Runtime (RDTSC) [s] |    21.9525 |
| Runtime unhalted [s] |    39.6194 |
|      Clock [MHz]     |  3798.2077 |
|          CPI         |     0.5508 |
|     DP [MFLOP/s]     |   177.2623 |
|   AVX DP [MFLOP/s]   |    53.8150 |
|  AVX512 DP [MFLOP/s] |          0 |
|   Packed [MUOPS/s]   |    32.5326 |
|   Scalar [MUOPS/s]   |    85.2895 |
|  Vectorization ratio |    27.6116 |
+----------------------+------------+
```

##### 1000 x 1000 humans on 100 x 100 km for minimum(100 steps, 10 minutes):
```
likwid-perfctr -f -C 0 -g FLOPS_DP swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: FLOPS_DP
+------------------------------------------+---------+----------------+
|                   Event                  | Counter |   HWThread 0   |
+------------------------------------------+---------+----------------+
|             INSTR_RETIRED_ANY            |  FIXC0  |  9763512748911 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  4519527714286 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  2380559070800 |
|               TOPDOWN_SLOTS              |  FIXC3  | 27117099211224 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |    12350095127 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |    13593594077 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |    10064623960 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |              0 |
+------------------------------------------+---------+----------------+

+----------------------+------------+
|        Metric        | HWThread 0 |
+----------------------+------------+
|  Runtime (RDTSC) [s] |  1195.7233 |
| Runtime unhalted [s] |  2259.8068 |
|      Clock [MHz]     |  3796.9584 |
|          CPI         |     0.4629 |
|     DP [MFLOP/s]     |    65.6944 |
|   AVX DP [MFLOP/s]   |    33.6687 |
|  AVX512 DP [MFLOP/s] |          0 |
|   Packed [MUOPS/s]   |    18.7457 |
|   Scalar [MUOPS/s]   |    11.3685 |
|  Vectorization ratio |    62.2487 |
+----------------------+------------+
```

##### 1000 x 1000 humans on 100 x 100 km for minimum, multiple rivers (1000 steps, 10 minutes):
```
likwid-perfctr -f -C 0 -g FLOPS_DP swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: FLOPS_DP
+--------------------------------------+---------+----------------+
|                 Event                | Counter |   HWThread 0   |
+--------------------------------------+---------+----------------+
|           INSTR_RETIRED_ANY          |  FIXC0  | 16579362600938 |
|         CPU_CLK_UNHALTED_CORE        |  FIXC1  | 13307840039017 |
|         CPU_CLK_UNHALTED_REF         |  FIXC2  | 10488495793516 |
| FP_COMP_OPS_EXE_SSE_FP_PACKED_DOUBLE |   PMC0  |    52287698182 |
| FP_COMP_OPS_EXE_SSE_FP_SCALAR_DOUBLE |   PMC1  |    69433493065 |
|       SIMD_FP_256_PACKED_DOUBLE      |   PMC2  |    82768760959 |
+--------------------------------------+---------+----------------+

+-------------------------+------------+
|          Metric         | HWThread 0 |
+-------------------------+------------+
|   Runtime (RDTSC) [s]   |  4080.0487 |
|   Runtime unhalted [s]  |  5118.4970 |
|       Clock [MHz]       |  3298.8265 |
|           CPI           |     0.8027 |
|       DP [MFLOP/s]      |   123.7936 |
|     AVX DP [MFLOP/s]    |    81.1449 |
|     Packed [MUOPS/s]    |    33.1017 |
|     Scalar [MUOPS/s]    |    17.0178 |
| Vectorization ratio [%] |    66.0455 |
+-------------------------+------------+
```

#### Comparison wrt scaling and with microbenchmarks

---

1. **Scaling: 100x100 Humans (10km x 10km) vs 1000x1000 Humans (100km x 100km)**

**Memory Metrics** (`MEM` group)

| Metric                        | 100x100 Humans | 1000x1000 Humans | Change                |
|-------------------------------|----------------|------------------|-----------------------|
| Runtime (RDTSC) [s]           | 19.42          | 1229.36          | ~63.3x increase       |
| Memory BW [MB/s]              | 103.06         | 205.97           | ~2x increase          |
| Memory Data Volume [GB]       | 2.00           | 253.21           | ~126.6x increase      |
| CPI                           | 0.5679         | 0.4615           | ~18.7% better         |

**Floating-Point Metrics** (`FLOPS_DP` group)

| Metric           | 100x100 Humans | 1000x1000 Humans | Change                |
|------------------|----------------|------------------|-----------------------|
| DP MFLOP/s       | 177.26         | 65.69            | ~63% decrease         |
| Vectorization [%]| 27.61          | 62.25            | ~2.25x improvement    |
| CPI              | 0.5508         | 0.4629           | ~16% improvement      |

---

2. **Comparison with `likwid-bench triad`**

| Metric         | SWIFT (best case) | likwid-bench triad | Ratio      |
|----------------|------------------|--------------------|------------|
| Memory BW      | 205.97 MB/s      | 26614.00 MB/s      | ~0.77%     |
| Read/Write Ratio | 1.1:1          | 3:1                | Lower      |
| CPI            | 0.4615           | 3.13               | ~6.8x better |

---

3. **Comparison with `likwid-bench peakflops`**

| Metric         | SWIFT (best case) | likwid-bench peakflops (single core) | likwid-bench peakflops (full node) | Ratio (single core) | Ratio (full node) |
|----------------|------------------|--------------------------------------|------------------------------------|---------------------|-------------------|
| MFLOP/s        | 177.26           | 4936.19                              | 44231.53                           | ~3.6%               | ~0.4%             |
| Vectorization  | 62.25%           | ~100%                                | ~100%                              | ~0.62x              | ~0.62x            |
| CPI            | 0.4629           | 8.43                                 | 0.94                               | ~18x better         | ~2x better        |

---

**Summary Table**

| Metric                | SWIFT (large, best case) | triad (single core) | triad (full node) | peakflops (single core) | peakflops (full node) |
|-----------------------|-------------------------|---------------------|-------------------|-------------------------|-----------------------|
| Memory BW [MB/s]      | 205.97                  | 14479.33            | 26614.00          | 2468.09                 | 22115.77              |
| DP MFLOP/s            | 65.69                   | (not measured)      | (not measured)    | 4936.19                 | 44231.53              |
| CPI                   | 0.4615                  | 3.13                | 3.13              | 8.43                    | 0.94                  |
| Vectorization [%]     | 62.25                   | (not measured)      | (not measured)    | ~100                    | ~100                  |

---

**Key Takeaways**

- **Scaling:** Memory bandwidth and data volume increase with problem size, but bandwidth is still far from hardware peak.
- **Compared to `triad`:** SWIFT achieves only a small fraction of the streaming bandwidth, with a more balanced (less streaming) access pattern and better CPI.
- **Compared to `peakflops`:** SWIFT achieves only ~3.6% of single-core peak FLOPS, and about 0.4% of full-node `peakflops` MFLOPS, with much lower vectorization.
- **Optimization potential:** There is significant headroom for improving both memory access patterns and vectorization in SWIFT.

---

### Parallel

#### Memory performance (`--enable-debug`, 4 MPI ranks, 4 OMP threads)

##### 100 x 100 humans on 10 x 10 km for minimum(10000 steps, 10 minutes):

```
likwid-mpirun -np 4 -t 4 -omp intel -g MEM -- swift_mpi --threads=4 -A -s -g -G --hm-river --hm-randomwalk -n 10000 ${HUMANMOBILITY}.yml
```
```
Group: 1
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+
|         Event         | Counter |   m5228:0:0  |   m5228:0:1  |   m5228:0:2  |   m5228:0:3  |   m5228:1:4  |   m5228:1:5  |   m5228:1:6  |   m5228:1:7  |   m5228:2:8  |   m5228:2:9  |  m5228:2:10  |  m5228:2:11  |  m5228:3:12  |  m5228:3:13  |  m5228:3:14  |  m5228:3:15  |
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+
|   INSTR_RETIRED_ANY   |  FIXC0  | 341001383176 |  94476288753 | 108228489509 | 157139864508 | 359764424707 |  92216965479 | 104581778267 | 150788089375 | 371625089799 |  94591782131 | 109177270245 | 154429978651 | 337094080890 |  96098187390 | 110280986563 | 159900454880 |
| CPU_CLK_UNHALTED_CORE |  FIXC1  | 280895787571 | 108450104598 |  74982138021 | 140211569134 | 288127266919 | 101224781427 |  72250116148 | 132304678131 | 293130918907 | 104687083572 |  73625866015 | 134097259742 | 286173688419 | 109087093895 |  75794547982 | 142937378290 |
|  CPU_CLK_UNHALTED_REF |  FIXC2  | 234079183520 |  92927087344 |  64614831190 | 120183334544 | 240140408820 |  86879553592 |  62273975608 | 113511073546 | 242431346300 |  90063364742 |  63294626902 | 115109635966 | 236550605798 |  93637848070 |  65134683822 | 122506849556 |
|      CAS_COUNT_RD     | MBOX0C0 |    395104151 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    391047867 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX0C1 |    323044630 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    301155310 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_RD     | MBOX1C0 |    378297602 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    403500528 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX1C1 |    298478851 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    307279648 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_RD     | MBOX2C0 |    389514168 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    511476132 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX2C1 |    309200166 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    363461718 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_RD     | MBOX3C0 |    372696559 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    463002788 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX3C1 |    296464156 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    306713669 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+

+----------------------------+---------+---------------+-------------+--------------+--------------+
|            Event           | Counter |      Sum      |     Min     |      Max     |      Avg     |
+----------------------------+---------+---------------+-------------+--------------+--------------+
|   INSTR_RETIRED_ANY STAT   |  FIXC0  | 2841395114323 | 92216965479 | 371625089799 | 1.775872e+11 |
| CPU_CLK_UNHALTED_CORE STAT |  FIXC1  | 2417980278771 | 72250116148 | 293130918907 | 1.511238e+11 |
|  CPU_CLK_UNHALTED_REF STAT |  FIXC2  | 2043338409320 | 62273975608 | 242431346300 | 1.277087e+11 |
|      CAS_COUNT_RD STAT     | MBOX0C0 |     786152018 |           0 |    395104151 | 4.913450e+07 |
|      CAS_COUNT_WR STAT     | MBOX0C1 |     624199940 |           0 |    323044630 | 3.901250e+07 |
|      CAS_COUNT_RD STAT     | MBOX1C0 |     781798130 |           0 |    403500528 | 4.886238e+07 |
|      CAS_COUNT_WR STAT     | MBOX1C1 |     605758499 |           0 |    307279648 | 3.785991e+07 |
|      CAS_COUNT_RD STAT     | MBOX2C0 |     900990300 |           0 |    511476132 | 5.631189e+07 |
|      CAS_COUNT_WR STAT     | MBOX2C1 |     672661884 |           0 |    363461718 | 4.204137e+07 |
|      CAS_COUNT_RD STAT     | MBOX3C0 |     835699347 |           0 |    463002788 | 5.223121e+07 |
|      CAS_COUNT_WR STAT     | MBOX3C1 |     603177825 |           0 |    306713669 | 3.769861e+07 |
+----------------------------+---------+---------------+-------------+--------------+--------------+

+-----------------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|               Metric              | m5228:0:0 | m5228:0:1 | m5228:0:2 | m5228:0:3 | m5228:1:4 | m5228:1:5 | m5228:1:6 | m5228:1:7 | m5228:2:8 | m5228:2:9 | m5228:2:10 | m5228:2:11 | m5228:3:12 | m5228:3:13 | m5228:3:14 | m5228:3:15 |
+-----------------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|        Runtime (RDTSC) [s]        |  136.6862 |  136.6862 |  136.6862 |  136.6862 |  136.6821 |  136.6821 |  136.6821 |  136.6821 |  136.6847 |  136.6847 |   136.6847 |   136.6847 |   136.6827 |   136.6827 |   136.6827 |   136.6827 |
|        Runtime unhalted [s]       |  108.0391 |   41.7124 |   28.8399 |   53.9286 |  110.8196 |   38.9331 |   27.7889 |   50.8871 |  112.7437 |   40.2646 |    28.3179 |    51.5763 |   110.0680 |    41.9570 |    29.1521 |    54.9765 |
|            Clock [MHz]            | 3119.9445 | 3034.2551 | 3017.1023 | 3033.2209 | 3119.5129 | 3029.2627 | 3016.4743 | 3030.4327 | 3143.7093 | 3022.1391 |  3024.3572 |  3028.8498 |  3145.3893 |  3028.9404 |  3025.4810 |  3033.5714 |
|                CPI                |    0.8237 |    1.1479 |    0.6928 |    0.8923 |    0.8009 |    1.0977 |    0.6908 |    0.8774 |    0.7888 |    1.1067 |     0.6744 |     0.8683 |     0.8489 |     1.1352 |     0.6873 |     0.8939 |
|  Memory read bandwidth [MBytes/s] |  719.0133 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  828.3133 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
|  Memory read data volume [GBytes] |   98.2792 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  113.2177 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
| Memory write bandwidth [MBytes/s] |  574.6009 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  598.6849 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
| Memory write data volume [GBytes] |   78.5400 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |   81.8311 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
|    Memory bandwidth [MBytes/s]    | 1293.6143 |         0 |         0 |         0 |         0 |         0 |         0 |         0 | 1426.9981 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
|    Memory data volume [GBytes]    |  176.8192 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  195.0488 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
+-----------------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+

+----------------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|                 Metric                 |     Sum    |    Min    |    Max    |    Avg    |  %ile 25  |  %ile 50  |  %ile 75  |
+----------------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |  2186.9428 |  136.6821 |  136.6862 |  136.6839 |  136.6821 |  136.6827 |  136.6847 |
|        Runtime unhalted [s] STAT       |   930.0048 |   27.7889 |  112.7437 |   58.1253 |  112.7437 |   29.1521 |   41.9570 |
|            Clock [MHz] STAT            | 48852.6429 | 3016.4743 | 3145.3893 | 3053.2902 | 3024.3572 | 3029.2627 | 3034.2551 |
|                CPI STAT                |    14.0270 |    0.6744 |    1.1479 |    0.8767 |    0.6928 |    0.8489 |    0.8939 |
|  Memory read bandwidth [MBytes/s] STAT |  1547.3266 |         0 |  828.3133 |   96.7079 |         0 |         0 |         0 |
|  Memory read data volume [GBytes] STAT |   211.4969 |         0 |  113.2177 |   13.2186 |         0 |         0 |         0 |
| Memory write bandwidth [MBytes/s] STAT |  1173.2858 |         0 |  598.6849 |   73.3304 |         0 |         0 |         0 |
| Memory write data volume [GBytes] STAT |   160.3711 |         0 |   81.8311 |   10.0232 |         0 |         0 |         0 |
|    Memory bandwidth [MBytes/s] STAT    |  2720.6124 |         0 | 1426.9981 |  170.0383 |         0 |         0 |         0 |
|    Memory data volume [GBytes] STAT    |   371.8680 |         0 |  195.0488 |   23.2417 |         0 |         0 |         0 |
+----------------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
```

##### 1000 x 1000 humans on 100 x 100 km for minimum(1000 steps, 10 minutes):

```
likwid-mpirun -np 4 -t 4 -omp intel -g MEM -- swift_mpi --threads=4 -A -s -g -G --hm-river --hm-randomwalk -n 10000 ${HUMANMOBILITY}.yml
```
```
Group: 1
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+
|         Event         | Counter |   m5237:0:0  |   m5237:0:1  |   m5237:0:2  |   m5237:0:3  |   m5237:1:4  |   m5237:1:5  |   m5237:1:6  |   m5237:1:7  |   m5237:2:8  |   m5237:2:9  |  m5237:2:10  |  m5237:2:11  |  m5237:3:12  |  m5237:3:13  |  m5237:3:14  |  m5237:3:15  |
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+
|   INSTR_RETIRED_ANY   |  FIXC0  | 935027159945 | 918078869462 | 945972539732 | 957959870956 | 896847519524 | 874744379907 | 924260839133 | 927599423107 | 916968570718 | 890387375529 | 933863665539 | 943080011762 | 903224228888 | 883709778741 | 928253499166 | 933640852265 |
| CPU_CLK_UNHALTED_CORE |  FIXC1  | 804755146770 | 779172644940 | 781786388839 | 789189897021 | 768855237917 | 744265028808 | 764818386042 | 763329436089 | 787163890777 | 758183098203 | 775728499290 | 775663816515 | 773399863303 | 748645688585 | 767735222132 | 764163093009 |
|  CPU_CLK_UNHALTED_REF |  FIXC2  | 696778006886 | 675063194728 | 677139994726 | 683414274686 | 665560928110 | 644816364608 | 662459649124 | 661041261192 | 680841627024 | 656554391870 | 671170753370 | 671005159318 | 669071116116 | 648364992652 | 664244851062 | 660990384080 |
|      CAS_COUNT_RD     | MBOX0C0 |   1303101073 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |   1173952283 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX0C1 |    604088681 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    516028925 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_RD     | MBOX1C0 |   1253823219 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |   1181399098 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX1C1 |    518406987 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    493431066 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_RD     | MBOX2C0 |   1249469599 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |   1279012528 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX2C1 |    547663411 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    570219086 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_RD     | MBOX3C0 |   1251631196 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |   1274434211 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
|      CAS_COUNT_WR     | MBOX3C1 |    517600181 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    493283369 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+

+----------------------------+---------+----------------+--------------+--------------+--------------+
|            Event           | Counter |       Sum      |      Min     |      Max     |      Avg     |
+----------------------------+---------+----------------+--------------+--------------+--------------+
|   INSTR_RETIRED_ANY STAT   |  FIXC0  | 14713618584374 | 874744379907 | 957959870956 | 9.196012e+11 |
| CPU_CLK_UNHALTED_CORE STAT |  FIXC1  | 12346855338240 | 744265028808 | 804755146770 | 771678458640 |
|  CPU_CLK_UNHALTED_REF STAT |  FIXC2  | 10688516949552 | 644816364608 | 696778006886 | 668032309347 |
|      CAS_COUNT_RD STAT     | MBOX0C0 |     2477053356 |            0 |   1303101073 | 1.548158e+08 |
|      CAS_COUNT_WR STAT     | MBOX0C1 |     1120117606 |            0 |    604088681 | 7.000735e+07 |
|      CAS_COUNT_RD STAT     | MBOX1C0 |     2435222317 |            0 |   1253823219 | 1.522014e+08 |
|      CAS_COUNT_WR STAT     | MBOX1C1 |     1011838053 |            0 |    518406987 | 6.323988e+07 |
|      CAS_COUNT_RD STAT     | MBOX2C0 |     2528482127 |            0 |   1279012528 | 1.580301e+08 |
|      CAS_COUNT_WR STAT     | MBOX2C1 |     1117882497 |            0 |    570219086 | 6.986766e+07 |
|      CAS_COUNT_RD STAT     | MBOX3C0 |     2526065407 |            0 |   1274434211 | 1.578791e+08 |
|      CAS_COUNT_WR STAT     | MBOX3C1 |     1010883550 |            0 |    517600181 | 6.318022e+07 |
+----------------------------+---------+----------------+--------------+--------------+--------------+

+-----------------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|               Metric              | m5237:0:0 | m5237:0:1 | m5237:0:2 | m5237:0:3 | m5237:1:4 | m5237:1:5 | m5237:1:6 | m5237:1:7 | m5237:2:8 | m5237:2:9 | m5237:2:10 | m5237:2:11 | m5237:3:12 | m5237:3:13 | m5237:3:14 | m5237:3:15 |
+-----------------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|        Runtime (RDTSC) [s]        |  306.1496 |  306.1496 |  306.1496 |  306.1496 |  306.1241 |  306.1241 |  306.1241 |  306.1241 |  306.1456 |  306.1456 |   306.1456 |   306.1456 |   306.1473 |   306.1473 |   306.1473 |   306.1473 |
|        Runtime unhalted [s]       |  309.5251 |  299.6855 |  300.6908 |  303.5383 |  295.7213 |  286.2632 |  294.1686 |  293.5959 |  302.7623 |  291.6156 |   298.3640 |   298.3391 |   297.4647 |   287.9437 |   295.2859 |   293.9120 |
|            Clock [MHz]            | 3002.8752 | 3000.9393 | 3001.7712 | 3002.3784 | 3003.4387 | 3000.9142 | 3001.6559 | 3002.2404 | 3005.9542 | 3002.3871 |  3004.9692 |  3005.4602 |  3005.3877 |  3002.1022 |  3005.0518 |  3005.7968 |
|                CPI                |    0.8607 |    0.8487 |    0.8264 |    0.8238 |    0.8573 |    0.8508 |    0.8275 |    0.8229 |    0.8584 |    0.8515 |     0.8307 |     0.8225 |     0.8563 |     0.8472 |     0.8271 |     0.8185 |
|  Memory read bandwidth [MBytes/s] | 1057.3707 |         0 |         0 |         0 |         0 |         0 |         0 |         0 | 1026.1885 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
|  Memory read data volume [GBytes] |  323.7136 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  314.1631 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
| Memory write bandwidth [MBytes/s] |  457.3470 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  433.3546 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
| Memory write data volume [GBytes] |  140.0166 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  132.6696 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
|    Memory bandwidth [MBytes/s]    | 1514.7176 |         0 |         0 |         0 |         0 |         0 |         0 |         0 | 1459.5430 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
|    Memory data volume [GBytes]    |  463.7302 |         0 |         0 |         0 |         0 |         0 |         0 |         0 |  446.8327 |         0 |          0 |          0 |          0 |          0 |          0 |          0 |
+-----------------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+

+----------------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|                 Metric                 |     Sum    |    Min    |    Max    |    Avg    |  %ile 25  |  %ile 50  |  %ile 75  |
+----------------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |  4898.2664 |  306.1241 |  306.1496 |  306.1416 |  306.1241 |  306.1456 |  306.1473 |
|        Runtime unhalted [s] STAT       |  4748.8760 |  286.2632 |  309.5251 |  296.8048 |  293.5959 |  295.7213 |  299.6855 |
|            Clock [MHz] STAT            | 48053.3225 | 3000.9142 | 3005.9542 | 3003.3327 | 3001.7712 | 3002.3871 | 3005.0518 |
|                CPI STAT                |    13.4303 |    0.8185 |    0.8607 |    0.8394 |    0.8238 |    0.8307 |    0.8515 |
|  Memory read bandwidth [MBytes/s] STAT |  2083.5592 |         0 | 1057.3707 |  130.2224 |         0 |         0 |         0 |
|  Memory read data volume [GBytes] STAT |   637.8767 |         0 |  323.7136 |   39.8673 |         0 |         0 |         0 |
| Memory write bandwidth [MBytes/s] STAT |   890.7016 |         0 |  457.3470 |   55.6688 |         0 |         0 |         0 |
| Memory write data volume [GBytes] STAT |   272.6862 |         0 |  140.0166 |   17.0429 |         0 |         0 |         0 |
|    Memory bandwidth [MBytes/s] STAT    |  2974.2606 |         0 | 1514.7176 |  185.8913 |         0 |         0 |         0 |
|    Memory data volume [GBytes] STAT    |   910.5629 |         0 |  463.7302 |   56.9102 |         0 |         0 |         0 |
+----------------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
```

#### Floating-point performance (`--enable-debug`)

##### 100 x 100 humans on 10 x 10 km for minimum(10000 steps, 10 minutes):

```
likwid-mpirun -np 4 -t 4 -omp intel -g FLOPS_DP -- swift_mpi --threads=4 -A -s -g -G --hm-river --hm-randomwalk -n 10000 ${HUMANMOBILITY}.yml
```
```
Group: 1
+--------------------------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+--------------+--------------+--------------+-------------+--------------+--------------+
|                 Event                | Counter |   m5335:0:0  |   m5335:0:1  |   m5335:0:2  |   m5335:0:3  |   m5335:1:4  |   m5335:1:5  |   m5335:1:6  |   m5335:1:7  |   m5335:2:8  |  m5335:2:9  |  m5335:2:10  |  m5335:2:11  |  m5335:3:12  |  m5335:3:13 |  m5335:3:14  |  m5335:3:15  |
+--------------------------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+--------------+--------------+--------------+-------------+--------------+--------------+
|           INSTR_RETIRED_ANY          |  FIXC0  | 321639846446 |  95295603041 | 108910848024 | 155472825044 | 313546724897 |  94256234895 | 108124055942 | 153092993255 | 354102823655 | 90172462296 | 103862094786 | 147648439350 | 348594216141 | 89095844083 | 102221345209 | 147222018909 |
|         CPU_CLK_UNHALTED_CORE        |  FIXC1  | 270939072666 | 108795059797 |  73441745166 | 136871060681 | 270188518834 | 106475416429 |  74103623128 | 136847769255 | 279981732694 | 99721523895 |  71746250642 | 129586517363 | 279131669583 | 98854040792 |  70657886881 | 128905664172 |
|         CPU_CLK_UNHALTED_REF         |  FIXC2  | 226196112298 |  93182929476 |  63301619290 | 117267745842 | 225514052218 |  91209295294 |  63836360484 | 117253940180 | 231285837224 | 85652307312 |  61703251714 | 111148431316 | 230557302768 | 84904931202 |  60771101846 | 110549719124 |
| FP_COMP_OPS_EXE_SSE_FP_PACKED_DOUBLE |   PMC0  |   1233656017 |    299866399 |    207790973 |    267129364 |   1436475677 |    283846907 |    192877546 |    260322176 |    849148674 |   293379462 |    203124677 |    267482063 |    754057509 |   286105969 |    196003709 |    258701926 |
| FP_COMP_OPS_EXE_SSE_FP_SCALAR_DOUBLE |   PMC1  |   4928924598 |   1450613650 |    605275577 |   1301031710 |   5716699045 |   1423932221 |    568589565 |   1313471659 |   3411173582 |  1341345940 |    586246720 |   1242421796 |   3063556224 |  1343623328 |    567353486 |   1228806835 |
|       SIMD_FP_256_PACKED_DOUBLE      |   PMC2  |    210446767 |    197850878 |    113108074 |    140181433 |    196810155 |    187607559 |    104398736 |    134087686 |    205509225 |   192591120 |    107716587 |    140119426 |    197362833 |   186098897 |    102038631 |    132977922 |
+--------------------------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+--------------+--------------+--------------+-------------+--------------+--------------+

+-------------------------------------------+---------+---------------+-------------+--------------+--------------+
|                   Event                   | Counter |      Sum      |     Min     |      Max     |      Avg     |
+-------------------------------------------+---------+---------------+-------------+--------------+--------------+
|           INSTR_RETIRED_ANY STAT          |  FIXC0  | 2733258375973 | 89095844083 | 354102823655 | 1.708286e+11 |
|         CPU_CLK_UNHALTED_CORE STAT        |  FIXC1  | 2336247551978 | 70657886881 | 279981732694 | 1.460155e+11 |
|         CPU_CLK_UNHALTED_REF STAT         |  FIXC2  | 1974334937588 | 60771101846 | 231285837224 | 1.233959e+11 |
| FP_COMP_OPS_EXE_SSE_FP_PACKED_DOUBLE STAT |   PMC0  |    7289969048 |   192877546 |   1436475677 | 4.556231e+08 |
| FP_COMP_OPS_EXE_SSE_FP_SCALAR_DOUBLE STAT |   PMC1  |   30093065936 |   567353486 |   5716699045 |   1880816621 |
|       SIMD_FP_256_PACKED_DOUBLE STAT      |   PMC2  |    2548905929 |   102038631 |    210446767 | 1.593066e+08 |
+-------------------------------------------+---------+---------------+-------------+--------------+--------------+

+-------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|          Metric         | m5335:0:0 | m5335:0:1 | m5335:0:2 | m5335:0:3 | m5335:1:4 | m5335:1:5 | m5335:1:6 | m5335:1:7 | m5335:2:8 | m5335:2:9 | m5335:2:10 | m5335:2:11 | m5335:3:12 | m5335:3:13 | m5335:3:14 | m5335:3:15 |
+-------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|   Runtime (RDTSC) [s]   |  132.3468 |  132.3468 |  132.3468 |  132.3468 |  132.3353 |  132.3353 |  132.3353 |  132.3353 |  132.3355 |  132.3355 |   132.3355 |   132.3355 |   132.3330 |   132.3330 |   132.3330 |   132.3330 |
|   Runtime unhalted [s]  |  104.2074 |   41.8443 |   28.2469 |   52.6428 |  103.9211 |   40.9530 |   28.5021 |   52.6350 |  107.6876 |   38.3553 |    27.5953 |    49.8421 |   107.3587 |    38.0209 |    27.1762 |    49.5793 |
|       Clock [MHz]       | 3114.2927 | 3035.6084 | 3016.4845 | 3034.6317 | 3114.9884 | 3035.1028 | 3018.1059 | 3034.4042 | 3147.3464 | 3027.0091 |  3023.1176 |  3031.2403 |  3147.7637 |  3027.1462 |  3022.9810 |  3031.7001 |
|           CPI           |    0.8424 |    1.1417 |    0.6743 |    0.8804 |    0.8617 |    1.1296 |    0.6854 |    0.8939 |    0.7907 |    1.1059 |     0.6908 |     0.8777 |     0.8007 |     1.1095 |     0.6912 |     0.8756 |
|       DP [MFLOP/s]      |   62.2457 |   21.4720 |   11.1320 |   18.1041 |   70.8571 |   20.7205 |   10.3671 |   17.9126 |   44.8217 |   20.3911 |    10.7557 |    17.6662 |    40.5124 |    20.1026 |    10.3339 |    17.2151 |
|     AVX DP [MFLOP/s]    |    6.3605 |    5.9798 |    3.4185 |    4.2368 |    5.9488 |    5.6707 |    3.1556 |    4.0530 |    6.2118 |    5.8213 |     3.2559 |     4.2353 |     5.9656 |     5.6252 |     3.0843 |     4.0195 |
|     Packed [MUOPS/s]    |   10.9115 |    3.7607 |    2.4247 |    3.0776 |   12.3420 |    3.5626 |    2.2464 |    2.9804 |    7.9696 |    3.6723 |     2.3489 |     3.0801 |     7.1896 |     3.5683 |     2.2522 |     2.9598 |
|     Scalar [MUOPS/s]    |   37.2425 |   10.9607 |    4.5734 |    9.8305 |   43.1986 |   10.7600 |    4.2966 |    9.9253 |   25.7767 |   10.1359 |     4.4300 |     9.3884 |    23.1504 |    10.1534 |     4.2873 |     9.2857 |
| Vectorization ratio [%] |   22.6596 |   25.5458 |   34.6478 |   23.8425 |   22.2216 |   24.8738 |   34.3328 |   23.0935 |   23.6162 |   26.5948 |    34.6500 |    24.7028 |    23.6968 |    26.0049 |    34.4400 |    24.1705 |
+-------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+

+------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|            Metric            |     Sum    |    Min    |    Max    |    Avg    |  %ile 25  |  %ile 50  |  %ile 75  |
+------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|   Runtime (RDTSC) [s] STAT   |  2117.4024 |  132.3330 |  132.3468 |  132.3376 |  132.3330 |  132.3353 |  132.3355 |
|   Runtime unhalted [s] STAT  |   898.5680 |   27.1762 |  107.6876 |   56.1605 |  107.6876 |   28.5021 |   41.8443 |
|       Clock [MHz] STAT       | 48861.9230 | 3016.4845 | 3147.7637 | 3053.8702 | 3023.1176 | 3031.7001 | 3035.6084 |
|           CPI STAT           |    14.0515 |    0.6743 |    1.1417 |    0.8782 |    0.6912 |    0.8617 |    0.8939 |
|       DP [MFLOP/s] STAT      |   414.6098 |   10.3339 |   70.8571 |   25.9131 |   11.1320 |   18.1041 |   21.4720 |
|     AVX DP [MFLOP/s] STAT    |    77.0426 |    3.0843 |    6.3605 |    4.8152 |    3.4185 |    4.2368 |    5.9488 |
|     Packed [MUOPS/s] STAT    |    74.3467 |    2.2464 |   12.3420 |    4.6467 |    2.2522 |    2.9804 |    3.5683 |
|     Scalar [MUOPS/s] STAT    |   227.3954 |    4.2873 |   43.1986 |   14.2122 |   10.9607 |    4.2873 |   43.1986 |
| Vectorization ratio [%] STAT |   429.0934 |   22.2216 |   34.6500 |   26.8183 |   23.6162 |   24.7028 |   26.5948 |
+------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
```

##### 1000 x 1000 humans on 100 x 100 km for minimum(1000 steps, 10 minutes):

```
likwid-mpirun -np 4 -t 4 -omp intel -g FLOPS_DP -- swift_mpi --threads=4 -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
```
```
Group: 1
+--------------------------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+
|                 Event                | Counter |   m5334:0:0  |   m5334:0:1  |   m5334:0:2  |   m5334:0:3  |   m5334:1:4  |   m5334:1:5  |   m5334:1:6  |   m5334:1:7  |   m5334:2:8  |   m5334:2:9  |  m5334:2:10  |  m5334:2:11  |  m5334:3:12  |  m5334:3:13  |  m5334:3:14  |   m5334:3:15  |
+--------------------------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+
|           INSTR_RETIRED_ANY          |  FIXC0  | 928603457884 | 911573592017 | 969853890424 | 975888043593 | 990161259348 | 966007358758 | 980630430015 | 997507812056 | 949733963252 | 928128991523 | 980345658060 | 996101066135 | 964132995603 | 940658959069 | 998496468941 | 1007029168610 |
|         CPU_CLK_UNHALTED_CORE        |  FIXC1  | 793709132880 | 770955305231 | 799673078449 | 799040042510 | 841996850660 | 816201744613 | 807755583411 | 814456886578 | 810874685654 | 782474264230 | 808541053373 | 813750612603 | 823452300321 | 793851472176 | 812653860819 |  816543101779 |
|         CPU_CLK_UNHALTED_REF         |  FIXC2  | 687190362924 | 667981690480 | 692733283216 | 692082389804 | 728841882236 | 707198316214 | 699755009798 | 705453081450 | 701540400028 | 677807241384 | 699089757886 | 703420064984 | 712007211786 | 687443295084 | 702251906564 |  705448866694 |
| FP_COMP_OPS_EXE_SSE_FP_PACKED_DOUBLE |   PMC0  |   2446903175 |   2539497409 |   2595236120 |   3596634676 |   2543018613 |   2653218612 |   2881304512 |   3667504788 |   2426244399 |   2644511043 |   2551783983 |   3716724177 |   2502062051 |   2633380750 |   2635341656 |    3560837119 |
| FP_COMP_OPS_EXE_SSE_FP_SCALAR_DOUBLE |   PMC1  |   3392789115 |   3292957442 |   3205400447 |   4053329366 |   3470664633 |   3477100065 |   3488209717 |   4221351979 |   3407043697 |   3406085381 |   3167849851 |   4203688028 |   3415853821 |   3450611329 |   3248905361 |    4075688959 |
|       SIMD_FP_256_PACKED_DOUBLE      |   PMC2  |   3797006147 |   3891979016 |   3924927923 |   3314325646 |   4333749324 |   4184528274 |   4175025409 |   3810679152 |   4032276636 |   4071116005 |   4021960270 |   3476932297 |   4103566526 |   4129436531 |   4058977364 |    3650775699 |
+--------------------------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+

+-------------------------------------------+---------+----------------+--------------+---------------+--------------+
|                   Event                   | Counter |       Sum      |      Min     |      Max      |      Avg     |
+-------------------------------------------+---------+----------------+--------------+---------------+--------------+
|           INSTR_RETIRED_ANY STAT          |  FIXC0  | 15484853115288 | 911573592017 | 1007029168610 | 9.678033e+11 |
|         CPU_CLK_UNHALTED_CORE STAT        |  FIXC1  | 12905929975287 | 770955305231 |  841996850660 | 8.066206e+11 |
|         CPU_CLK_UNHALTED_REF STAT         |  FIXC2  | 11170244760532 | 667981690480 |  728841882236 | 6.981403e+11 |
| FP_COMP_OPS_EXE_SSE_FP_PACKED_DOUBLE STAT |   PMC0  |    45594203083 |   2426244399 |    3716724177 | 2.849638e+09 |
| FP_COMP_OPS_EXE_SSE_FP_SCALAR_DOUBLE STAT |   PMC1  |    56977529191 |   3167849851 |    4221351979 | 3.561096e+09 |
|       SIMD_FP_256_PACKED_DOUBLE STAT      |   PMC2  |    62977262219 |   3314325646 |    4333749324 | 3.936079e+09 |
+-------------------------------------------+---------+----------------+--------------+---------------+--------------+

+-------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|          Metric         | m5334:0:0 | m5334:0:1 | m5334:0:2 | m5334:0:3 | m5334:1:4 | m5334:1:5 | m5334:1:6 | m5334:1:7 | m5334:2:8 | m5334:2:9 | m5334:2:10 | m5334:2:11 | m5334:3:12 | m5334:3:13 | m5334:3:14 | m5334:3:15 |
+-------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+
|   Runtime (RDTSC) [s]   |  323.8998 |  323.8998 |  323.8998 |  323.8998 |  323.9200 |  323.9200 |  323.9200 |  323.9200 |  323.9356 |  323.9356 |   323.9356 |   323.9356 |   323.9037 |   323.9037 |   323.9037 |   323.9037 |
|   Runtime unhalted [s]  |  305.2741 |  296.5226 |  307.5680 |  307.3245 |  323.8471 |  313.9259 |  310.6773 |  313.2548 |  311.8845 |  300.9610 |   310.9870 |   312.9907 |   316.7138 |   305.3289 |   312.5606 |   314.0564 |
|       Clock [MHz]       | 3003.0024 | 3000.7928 | 3001.3579 | 3001.8024 | 3003.6377 | 3000.7285 | 3001.2650 | 3001.7213 | 3005.1141 | 3001.3992 |  3006.9696 |  3007.7136 |  3006.9449 |  3002.4362 |  3008.7360 |  3009.4351 |
|           CPI           |    0.8547 |    0.8457 |    0.8245 |    0.8188 |    0.8504 |    0.8449 |    0.8237 |    0.8165 |    0.8538 |    0.8431 |     0.8248 |     0.8169 |     0.8541 |     0.8439 |     0.8139 |     0.8108 |
|       DP [MFLOP/s]      |   72.4749 |   73.9113 |   74.3921 |   75.6527 |   79.9324 |   78.7900 |   80.1152 |   82.7336 |   75.2885 |   77.1128 |    75.1978 |    78.8579 |    76.6717 |    77.9093 |    76.4286 |    79.6547 |
|     AVX DP [MFLOP/s]    |   46.8911 |   48.0640 |   48.4709 |   40.9303 |   53.5163 |   51.6736 |   51.5563 |   47.0570 |   49.7911 |   50.2707 |    49.6637 |    42.9336 |    50.6764 |    50.9959 |    50.1257 |    45.0847 |
|     Packed [MUOPS/s]    |   19.2773 |   19.8564 |   20.1302 |   21.3367 |   21.2298 |   21.1094 |   21.7842 |   23.0865 |   19.9377 |   20.7314 |    20.2934 |    22.2071 |    20.3938 |    20.8791 |    20.6676 |    22.2647 |
|     Scalar [MUOPS/s]    |   10.4748 |   10.1666 |    9.8963 |   12.5141 |   10.7146 |   10.7344 |   10.7687 |   13.0321 |   10.5177 |   10.5147 |     9.7793 |    12.9769 |    10.5459 |    10.6532 |    10.0305 |    12.5830 |
| Vectorization ratio [%] |   64.7930 |   66.1373 |   67.0415 |   63.0315 |   66.4587 |   66.2903 |   66.9193 |   63.9186 |   65.4653 |   66.3487 |    67.4812 |    63.1170 |    65.9147 |    66.2150 |    67.3254 |    63.8914 |
+-------------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+------------+------------+------------+------------+------------+

+------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|            Metric            |     Sum    |    Min    |    Max    |    Avg    |  %ile 25  |  %ile 50  |  %ile 75  |
+------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
|   Runtime (RDTSC) [s] STAT   |  5182.6364 |  323.8998 |  323.9356 |  323.9148 |  323.8998 |  323.9037 |  323.9200 |
|   Runtime unhalted [s] STAT  |  4963.8772 |  296.5226 |  323.8471 |  310.2423 |  305.3289 |  310.9870 |  313.2548 |
|       Clock [MHz] STAT       | 48063.0567 | 3000.7285 | 3009.4351 | 3003.9410 | 3001.3579 | 3002.4362 | 3006.9449 |
|           CPI STAT           |    13.3405 |    0.8108 |    0.8547 |    0.8338 |    0.8169 |    0.8248 |    0.8457 |
|       DP [MFLOP/s] STAT      |  1235.1235 |   72.4749 |   82.7336 |   77.1952 |   75.1978 |   76.6717 |   78.8579 |
|     AVX DP [MFLOP/s] STAT    |   777.7013 |   40.9303 |   53.5163 |   48.6063 |   46.8911 |   49.6637 |   50.6764 |
|     Packed [MUOPS/s] STAT    |   335.1853 |   19.2773 |   23.0865 |   20.9491 |   20.1302 |   20.7314 |   21.3367 |
|     Scalar [MUOPS/s] STAT    |   175.9028 |    9.7793 |   13.0321 |   10.9939 |   10.5147 |   10.7146 |   12.5830 |
| Vectorization ratio [%] STAT |  1050.3489 |   63.0315 |   67.4812 |   65.6468 |   63.9186 |   66.1373 |   66.4587 |
+------------------------------+------------+-----------+-----------+-----------+-----------+-----------+-----------+
```

#### Comparison wrt scaling and with microbenchmarks

---

1. **Scaling**: 100x100 Humans (10km x 10km) vs 1000x1000 Humans (100km x 100km)

**Memory Metrics** (`MEM` group, best rank)

| Metric                        | 100x100 Humans | 1000x1000 Humans | Change                |
|-------------------------------|----------------|------------------|-----------------------|
| Runtime (RDTSC) [s]           | 136.69         | 306.15           | ~2.24x increase       |
| Memory BW [MB/s]              | 1427           | 1515             | Slight increase       |
| Memory Data Volume [GB]       | 195            | 464              | ~2.38x increase       |
| CPI                           | 0.82           | 0.84             | Similar               |

**Floating-Point Metrics** (`FLOPS_DP` group, best rank)

| Metric           | 100x100 Humans | 1000x1000 Humans | Change                |
|------------------|----------------|------------------|-----------------------|
| DP MFLOP/s       | ~62            | ~44              | ~29% decrease         |
| Vectorization [%]| ~22–34         | ~22–34           | Similar (low)         |
| CPI              | ~0.84          | ~0.84            | Similar               |

---

2. **Comparison with `likwid-bench triad`**

| Metric         | swift_mpi (best rank) | likwid-bench triad | Ratio      |
|----------------|----------------------|--------------------|------------|
| Memory BW      | 1515 MB/s            | 26614.00 MB/s      | ~5.7%      |
| Read/Write Ratio | ~1.8:1             | 3:1                | Lower      |
| CPI            | 0.84                 | 3.13               | ~3.7x better |

- **Observation:**  
  - SWIFT achieves only ~5.7% of the streaming memory bandwidth of `triad`.
  - Access pattern is less streaming (lower read/write ratio).
  - CPI is better in SWIFT, but this is expected for memory-bound code.

---

3. **Comparison with `likwid-bench peakflops`**

| Metric         | swift_mpi (best rank) | likwid-bench peakflops | Ratio      |
|----------------|----------------------|------------------------|------------|
| MFLOP/s        | ~82.7                | 4936.19                | ~1.7%      |
| Memory BW      | 1515 MB/s            | 2468.09 MB/s           | ~61%       |
| CPI            | 0.84                 | 8.43                   | ~10x better |
| Vectorization  | ~67%                 | ~100%                  | ~0.67x     |

- **Observation:**  
  - SWIFT achieves about 1.7% of peak FLOPS and about 61% of `peakflops` memory bandwidth.
  - Vectorization ratio is much lower than `peakflops`.
  - CPI is much better in SWIFT, but `peakflops` is compute-bound.

---

**Summary Table**

| Metric                | swift_mpi (large, best rank) | triad             | peakflops         |
|-----------------------|------------------------------|-------------------|-------------------|
| Memory BW [MB/s]      | 1515                         | 26614.00          | 2468.09           |
| DP MFLOP/s            | ~82.7                        | (not measured)    | 4936.19           |
| CPI                   | 0.84                         | 3.13              | 8.43              |
| Vectorization [%]     | ~67                          | (not measured)    | ~100              |

---

**Key Takeaways**

- **Scaling:** Memory bandwidth and data volume increase with problem size, but bandwidth is still far from hardware peak.
- **Compared to `triad`:** SWIFT achieves only a small fraction of the streaming bandwidth, with a more balanced (less streaming) access pattern and better CPI.
- **Compared to `peakflops`:** SWIFT achieves about 1.7% of peak FLOPS, and about 61% of `peakflops` memory bandwidth, with much lower vectorization.
- **Optimization potential:** There is significant headroom for improving both memory access patterns and vectorization in SWIFT.

#### Memory performance (`--enable-ipo`)

##### 1000 x 1000 humans on 100 x 100 km for minimum, multiple rivers (100 steps, 10 minutes):
```
likwid-perfctr -f -C 0-15 -g MEM swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 1000 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: MEM
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+--------------+--------------+--------------+--------------+--------------+---------------+--------------+
|         Event         | Counter |  HWThread 0  |  HWThread 1  |  HWThread 2  |  HWThread 3  |  HWThread 4  |  HWThread 5  |  HWThread 6  |  HWThread 7  |   HWThread 8  |  HWThread 9  |  HWThread 10 |  HWThread 11 |  HWThread 12 |  HWThread 13 |  HWThread 14  |  HWThread 15 |
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+--------------+--------------+--------------+--------------+--------------+---------------+--------------+
|   INSTR_RETIRED_ANY   |  FIXC0  | 975876943752 | 965752585140 | 993753825451 | 979073598842 | 970379226841 | 978056530864 | 963295174486 | 978504614942 | 1025397969559 | 999298716967 | 993478925405 | 998344253860 | 974925642167 | 986585240716 | 1008142361518 | 958684471716 |
| CPU_CLK_UNHALTED_CORE |  FIXC1  | 887566984243 | 875916485053 | 900089754704 | 886926244471 | 879167838974 | 884625668090 | 871514404512 | 880467057248 |  922567316942 | 902274334852 | 897231000281 | 899206842936 | 881539711541 | 890218795683 |  905108770926 | 849800952747 |
|  CPU_CLK_UNHALTED_REF |  FIXC2  | 768140274512 | 759043739714 | 779977106012 | 768588381990 | 761894173040 | 766590485856 | 755228198478 | 762985594280 |  799331217282 | 781861953288 | 777521797962 | 779256954554 | 763927306012 | 771449474224 |  784324764002 | 736409433578 |
|      CAS_COUNT_RD     | MBOX0C0 |   4261444731 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    3504353741 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
|      CAS_COUNT_WR     | MBOX0C1 |   1512708396 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |     943101388 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
|      CAS_COUNT_RD     | MBOX1C0 |   4726480585 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    4403350275 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
|      CAS_COUNT_WR     | MBOX1C1 |   1539877534 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    1415800887 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
|      CAS_COUNT_RD     | MBOX2C0 |   9218646904 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    5702862507 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
|      CAS_COUNT_WR     | MBOX2C1 |   4397646832 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    1439969824 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
|      CAS_COUNT_RD     | MBOX3C0 |  10524023451 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    8511873950 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
|      CAS_COUNT_WR     | MBOX3C1 |   2136606278 |            0 |            0 |            0 |            0 |            0 |            0 |            0 |    3398124012 |            0 |            0 |            0 |            0 |            0 |             0 |            0 |
+-----------------------+---------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+--------------+--------------+--------------+--------------+--------------+---------------+--------------+

+----------------------------+---------+----------------+--------------+---------------+--------------+
|            Event           | Counter |       Sum      |      Min     |      Max      |      Avg     |
+----------------------------+---------+----------------+--------------+---------------+--------------+
|   INSTR_RETIRED_ANY STAT   |  FIXC0  | 15749550082226 | 958684471716 | 1025397969559 | 9.843469e+11 |
| CPU_CLK_UNHALTED_CORE STAT |  FIXC1  | 14214222163203 | 849800952747 |  922567316942 | 8.883889e+11 |
|  CPU_CLK_UNHALTED_REF STAT |  FIXC2  | 12316530854784 | 736409433578 |  799331217282 | 769783178424 |
|      CAS_COUNT_RD STAT     | MBOX0C0 |     7765798472 |            0 |    4261444731 | 4.853624e+08 |
|      CAS_COUNT_WR STAT     | MBOX0C1 |     2455809784 |            0 |    1512708396 | 1.534881e+08 |
|      CAS_COUNT_RD STAT     | MBOX1C0 |     9129830860 |            0 |    4726480585 | 5.706144e+08 |
|      CAS_COUNT_WR STAT     | MBOX1C1 |     2955678421 |            0 |    1539877534 | 1.847299e+08 |
|      CAS_COUNT_RD STAT     | MBOX2C0 |    14921509411 |            0 |    9218646904 | 9.325943e+08 |
|      CAS_COUNT_WR STAT     | MBOX2C1 |     5837616656 |            0 |    4397646832 |    364851041 |
|      CAS_COUNT_RD STAT     | MBOX3C0 |    19035897401 |            0 |   10524023451 | 1.189744e+09 |
|      CAS_COUNT_WR STAT     | MBOX3C1 |     5534730290 |            0 |    3398124012 | 3.459206e+08 |
+----------------------------+---------+----------------+--------------+---------------+--------------+

+-----------------------------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-------------+-------------+-------------+-------------+-------------+-------------+
|               Metric              | HWThread 0 | HWThread 1 | HWThread 2 | HWThread 3 | HWThread 4 | HWThread 5 | HWThread 6 | HWThread 7 | HWThread 8 | HWThread 9 | HWThread 10 | HWThread 11 | HWThread 12 | HWThread 13 | HWThread 14 | HWThread 15 |
+-----------------------------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-------------+-------------+-------------+-------------+-------------+-------------+
|        Runtime (RDTSC) [s]        |  1025.7361 |  1025.7361 |  1025.7361 |  1025.7361 |  1025.7361 |  1025.7361 |  1025.7361 |  1025.7361 |  1025.7361 |  1025.7361 |   1025.7361 |   1025.7361 |   1025.7361 |   1025.7361 |   1025.7361 |   1025.7361 |
|        Runtime unhalted [s]       |   341.3757 |   336.8947 |   346.1922 |   341.1293 |   338.1452 |   340.2444 |   335.2016 |   338.6449 |   354.8375 |   347.0324 |    345.0927 |    345.8526 |    339.0575 |    342.3956 |    348.1226 |    326.8501 |
|            Clock [MHz]            |  3004.2020 |  3000.2982 |  3000.3539 |  3000.2830 |  3000.1686 |  3000.2997 |  3000.3015 |  3000.3043 |  3000.8191 |  3000.3854 |   3000.2692 |   3000.1809 |   3000.2564 |   3000.2524 |   3000.3600 |   3000.3118 |
|                CPI                |     0.9095 |     0.9070 |     0.9057 |     0.9059 |     0.9060 |     0.9045 |     0.9047 |     0.8998 |     0.8997 |     0.9029 |      0.9031 |      0.9007 |      0.9042 |      0.9023 |      0.8978 |      0.8864 |
|  Memory read bandwidth [MBytes/s] |  1792.6230 |          0 |          0 |          0 |          0 |          0 |          0 |          0 |  1380.3123 |          0 |           0 |           0 |           0 |           0 |           0 |           0 |
|  Memory read data volume [GBytes] |  1838.7581 |          0 |          0 |          0 |          0 |          0 |          0 |          0 |  1415.8362 |          0 |           0 |           0 |           0 |           0 |           0 |           0 |
| Memory write bandwidth [MBytes/s] |   598.1633 |          0 |          0 |          0 |          0 |          0 |          0 |          0 |   449.0509 |          0 |           0 |           0 |           0 |           0 |           0 |           0 |
| Memory write data volume [GBytes] |   613.5577 |          0 |          0 |          0 |          0 |          0 |          0 |          0 |   460.6078 |          0 |           0 |           0 |           0 |           0 |           0 |           0 |
|    Memory bandwidth [MBytes/s]    |  2390.7863 |          0 |          0 |          0 |          0 |          0 |          0 |          0 |  1829.3633 |          0 |           0 |           0 |           0 |           0 |           0 |           0 |
|    Memory data volume [GBytes]    |  2452.3158 |          0 |          0 |          0 |          0 |          0 |          0 |          0 |  1876.4439 |          0 |           0 |           0 |           0 |           0 |           0 |           0 |
+-----------------------------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-------------+-------------+-------------+-------------+-------------+-------------+

+----------------------------------------+------------+-----------+-----------+-----------+
|                 Metric                 |     Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+------------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        | 16411.7776 | 1025.7361 | 1025.7361 | 1025.7361 |
|        Runtime unhalted [s] STAT       |  5467.0690 |  326.8501 |  354.8375 |  341.6918 |
|            Clock [MHz] STAT            | 48009.0464 | 3000.1686 | 3004.2020 | 3000.5654 |
|                CPI STAT                |    14.4402 |    0.8864 |    0.9095 |    0.9025 |
|  Memory read bandwidth [MBytes/s] STAT |  3172.9353 |         0 | 1792.6230 |  198.3085 |
|  Memory read data volume [GBytes] STAT |  3254.5943 |         0 | 1838.7581 |  203.4121 |
| Memory write bandwidth [MBytes/s] STAT |  1047.2142 |         0 |  598.1633 |   65.4509 |
| Memory write data volume [GBytes] STAT |  1074.1655 |         0 |  613.5577 |   67.1353 |
|    Memory bandwidth [MBytes/s] STAT    |  4220.1496 |         0 | 2390.7863 |  263.7594 |
|    Memory data volume [GBytes] STAT    |  4328.7597 |         0 | 2452.3158 |  270.5475 |
+----------------------------------------+------------+-----------+-----------+-----------+
```

#### Floating-point performance (`--enable-ipo`)

##### 1000 x 1000 humans on 100 x 100 km for minimum, multiple rivers (1000 steps, 10 minutes):
```
likwid-perfctr -f -C 0 -g FLOPS_DP swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 humanMobility.yml
```
```
--------------------------------------------------------------------------------
Group 1: FLOPS_DP
+--------------------------------------+---------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+
|                 Event                | Counter |   HWThread 0  |   HWThread 1  |   HWThread 2  |   HWThread 3  |   HWThread 4  |   HWThread 5  |   HWThread 6  |   HWThread 7  |   HWThread 8  |   HWThread 9  |  HWThread 10  |  HWThread 11  |  HWThread 12  |  HWThread 13  |  HWThread 14  |  HWThread 15  |
+--------------------------------------+---------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+
|           INSTR_RETIRED_ANY          |  FIXC0  | 1062384074352 | 1036415071743 | 1052600966172 | 1054789692902 | 1067560654541 | 1047140772969 | 1063377721926 | 1070395853218 | 1039607704349 | 1039203944144 | 1021167287992 | 1053327830351 | 1029103907737 | 1065359127658 | 1035985100198 | 1058987431644 |
|         CPU_CLK_UNHALTED_CORE        |  FIXC1  |  962368120170 |  937971084671 |  951319896126 |  952615880874 |  962674274203 |  942721837444 |  958001948413 |  963713542269 |  941090357070 |  937104057117 |  924319305694 |  949240686047 |  928964172167 |  959757595020 |  933028997627 |  936503975787 |
|         CPU_CLK_UNHALTED_REF         |  FIXC2  |  832954836298 |  812828889132 |  824416687198 |  825565709150 |  834282085390 |  816969631036 |  830207617928 |  835157406044 |  815491942954 |  812137083732 |  801021984958 |  822586472682 |  805038710476 |  831556644412 |  808526438616 |  811535268284 |
| FP_COMP_OPS_EXE_SSE_FP_PACKED_DOUBLE |   PMC0  |    3590712378 |    3346231008 |    3406454917 |    3453499476 |    3557858528 |    3403706927 |    3497341900 |    3419434893 |    3557593240 |    3426335665 |    3473419347 |    3509474601 |    3294336813 |    3550051922 |    3434509625 |    3434216267 |
| FP_COMP_OPS_EXE_SSE_FP_SCALAR_DOUBLE |   PMC1  |    4898009353 |    4427932292 |    4503997056 |    4539411167 |    4682085388 |    4502844923 |    4601153705 |    4517165183 |    4576502676 |    4513198827 |    4522478313 |    4609811691 |    4314254144 |    4662173507 |    4508650566 |    4232780054 |
|       SIMD_FP_256_PACKED_DOUBLE      |   PMC2  |    5395188163 |    5189979707 |    5280222115 |    5175091452 |    5270021723 |    5318743192 |    5209055517 |    5231177899 |    5066647603 |    5242994284 |    5184493447 |    5258492332 |    5032011937 |    5335969911 |    5365858948 |    5006897884 |
+--------------------------------------+---------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+

+-------------------------------------------+---------+----------------+---------------+---------------+--------------+
|                   Event                   | Counter |       Sum      |      Min      |      Max      |      Avg     |
+-------------------------------------------+---------+----------------+---------------+---------------+--------------+
|           INSTR_RETIRED_ANY STAT          |  FIXC0  | 16797407141896 | 1021167287992 | 1070395853218 | 1.049838e+12 |
|         CPU_CLK_UNHALTED_CORE STAT        |  FIXC1  | 15141395730699 |  924319305694 |  963713542269 | 9.463372e+11 |
|         CPU_CLK_UNHALTED_REF STAT         |  FIXC2  | 13120277408290 |  801021984958 |  835157406044 | 8.200173e+11 |
| FP_COMP_OPS_EXE_SSE_FP_PACKED_DOUBLE STAT |   PMC0  |    55355177507 |    3294336813 |    3590712378 | 3.459699e+09 |
| FP_COMP_OPS_EXE_SSE_FP_SCALAR_DOUBLE STAT |   PMC1  |    72612448845 |    4232780054 |    4898009353 | 4.538278e+09 |
|       SIMD_FP_256_PACKED_DOUBLE STAT      |   PMC2  |    83562846114 |    5006897884 |    5395188163 | 5.222678e+09 |
+-------------------------------------------+---------+----------------+---------------+---------------+--------------+

+-------------------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-------------+-------------+-------------+-------------+-------------+-------------+
|          Metric         | HWThread 0 | HWThread 1 | HWThread 2 | HWThread 3 | HWThread 4 | HWThread 5 | HWThread 6 | HWThread 7 | HWThread 8 | HWThread 9 | HWThread 10 | HWThread 11 | HWThread 12 | HWThread 13 | HWThread 14 | HWThread 15 |
+-------------------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-------------+-------------+-------------+-------------+-------------+-------------+
|   Runtime (RDTSC) [s]   |  1115.6544 |  1115.6544 |  1115.6544 |  1115.6544 |  1115.6544 |  1115.6544 |  1115.6544 |  1115.6544 |  1115.6544 |  1115.6544 |   1115.6544 |   1115.6544 |   1115.6544 |   1115.6544 |   1115.6544 |   1115.6544 |
|   Runtime unhalted [s]  |   370.1459 |   360.7623 |   365.8965 |   366.3950 |   370.2637 |   362.5896 |   368.4666 |   370.6634 |   361.9621 |   360.4288 |    355.5116 |    365.0968 |    357.2981 |    369.1418 |    358.8615 |    360.1980 |
|       Clock [MHz]       |  3003.9178 |  3000.2580 |  3000.1853 |  3000.0911 |  3000.0930 |  3000.1705 |  3000.1844 |  3000.1841 |  3000.4053 |  3000.0380 |   3000.1700 |   3000.2888 |   3000.2019 |   3000.8065 |   3000.3312 |   3000.3402 |
|           CPI           |     0.9059 |     0.9050 |     0.9038 |     0.9031 |     0.9018 |     0.9003 |     0.9009 |     0.9003 |     0.9052 |     0.9018 |      0.9052 |      0.9012 |      0.9027 |      0.9009 |      0.9006 |      0.8843 |
|       DP [MFLOP/s]      |    30.1708 |    28.5754 |    29.0751 |    28.8143 |    29.4696 |    29.2073 |    29.0700 |    28.9344 |    28.6453 |    28.9855 |     28.8685 |     29.2767 |     27.8141 |     29.6742 |     29.4366 |     27.9018 |
|     AVX DP [MFLOP/s]    |    19.3436 |    18.6078 |    18.9314 |    18.5545 |    18.8948 |    19.0695 |    18.6762 |    18.7555 |    18.1657 |    18.7979 |     18.5882 |     18.8535 |     18.0415 |     19.1313 |     19.2384 |     17.9514 |
|     Packed [MUOPS/s]    |     8.0544 |     7.6513 |     7.7862 |     7.7341 |     7.9127 |     7.8182 |     7.8038 |     7.7538 |     7.7302 |     7.7706 |      7.7604 |      7.8590 |      7.4632 |      7.9649 |      7.8881 |      7.5661 |
|     Scalar [MUOPS/s]    |     4.3903 |     3.9689 |     4.0371 |     4.0688 |     4.1967 |     4.0361 |     4.1242 |     4.0489 |     4.1021 |     4.0453 |      4.0537 |      4.1319 |      3.8670 |      4.1789 |      4.0413 |      3.7940 |
| Vectorization ratio [%] |    64.7217 |    65.8448 |    65.8547 |    65.5270 |    65.3435 |    65.9528 |    65.4245 |    65.6953 |    65.3315 |    65.7638 |     65.6878 |     65.5413 |     65.8699 |     65.5882 |     66.1233 |     66.6024 |
+-------------------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-------------+-------------+-------------+-------------+-------------+-------------+

+------------------------------+------------+-----------+-----------+-----------+
|            Metric            |     Sum    |    Min    |    Max    |    Avg    |
+------------------------------+------------+-----------+-----------+-----------+
|   Runtime (RDTSC) [s] STAT   | 17850.4704 | 1115.6544 | 1115.6544 | 1115.6544 |
|   Runtime unhalted [s] STAT  |  5823.6817 |  355.5116 |  370.6634 |  363.9801 |
|       Clock [MHz] STAT       | 48007.6661 | 3000.0380 | 3003.9178 | 3000.4791 |
|           CPI STAT           |    14.4230 |    0.8843 |    0.9059 |    0.9014 |
|       DP [MFLOP/s] STAT      |   463.9196 |   27.8141 |   30.1708 |   28.9950 |
|     AVX DP [MFLOP/s] STAT    |   299.6012 |   17.9514 |   19.3436 |   18.7251 |
|     Packed [MUOPS/s] STAT    |   124.5170 |    7.4632 |    8.0544 |    7.7823 |
|     Scalar [MUOPS/s] STAT    |    65.0852 |    3.7940 |    4.3903 |    4.0678 |
| Vectorization ratio [%] STAT |  1050.8725 |   64.7217 |   66.6024 |   65.6795 |
+------------------------------+------------+-----------+-----------+-----------+
```

#### Comparison with previous results (`--enable-debug`) and with microbenchmarks

---

**Memory Metrics (`MEM` group, best rank):**

| Metric                        | `--enable-debug`              | `--enable-ipo`   | Change                |
|                               | (4 MPI ranks × 4 OMP threads) | (16 OMP threads) |                       |
|-------------------------------|-------------------------------|------------------|-----------------------|
| Runtime (RDTSC) [s]           | 306.15                        | 1025.74          | ~3.35x increase       |
| Memory BW [MB/s]              | 1515                          | 264              | ~5.7x decrease        |
| Memory Data Volume [GB]       | 464                           | 270              | ~42% decrease         |
| CPI                           | 0.84                          | 0.90             | Slightly worse        |

- **Observation:**  
  - With IPO enabled, runtime increases and memory bandwidth drops significantly compared to the debug build.
  - CPI is slightly worse, indicating less efficient instruction execution.
  - Data volume is lower, possibly due to different run lengths or internal optimizations.

---

**Floating-Point Metrics (`FLOPS_DP` group, best rank):**

| Metric           | SWIFT(`--enable-debug`)       | SWIFT(`--enable-ipo`) | Change                |
|                  | (4 MPI ranks × 4 OMP threads) | (16 OMP threads)      |                       |
|------------------|-------------------------------|-----------------------|-----------------------|
| DP MFLOP/s       | ~77                           | ~29                   | ~2.7x decrease        |
| Vectorization [%]| ~66                           | ~66                   | Similar               |
| CPI              | ~0.84                         | ~0.90                 | Slightly worse        |

- **Observation:**  
  - Floating-point throughput drops with IPO enabled.
  - Vectorization ratio remains similar.
  - CPI is slightly worse.

---

**Comparison with Microbenchmarks:**

| Metric         | SWIFT(`--enable-ipo`) | SWIFT(`--enable-debug`)       | likwid-bench triad | likwid-bench peakflops |
|                | (16 OMP threads)      | (4 MPI ranks × 4 OMP threads) |                    |                        |
|----------------|-----------------------|-------------------------------|--------------------|------------------------|
| Memory BW      | 264 MB/s              | 1515 MB/s                     | 26614 MB/s         | 2468 MB/s              |
| DP MFLOP/s     | ~29                   | ~77                           | (not measured)     | 4936                   |
| Vectorization  | ~66%                  | ~66%                          | (not measured)     | ~100%                  |
| CPI            | 0.90                  | 0.84                          | 3.13 (triad)       | 8.43 (peakflops)       |

- **Observation:**  
  - Both SWIFT builds achieve only a small fraction of the hardware's peak memory bandwidth and floating-point throughput.
  - IPO does not improve performance for this workload; in fact, it reduces both memory and compute throughput.
  - Vectorization is much lower than in microbenchmarks.
  - CPI is much better than microbenchmarks, but this is typical for memory-bound codes.

---

**Key Takeaways**

- **IPO did not improve performance** for this configuration; memory and floating-point throughput are lower than with `--enable-debug`.
- **SWIFT remains far from hardware peak** for both memory and compute, with vectorization much lower than microbenchmarks.
- **Optimization potential remains high** for memory access patterns and vectorization in SWIFT.

---

### MAQAO Performance Analysis

#### Runtime Comparison

| Tool         | Configuration                                  | Runtime [s] | Notes                   |
|--------------|------------------------------------------------|-------------|-------------------------|
| **MAQAO**    | (4 MPI ranks × 4 OMP threads)                  | 294.31      | 1000 humans, 1000 steps |
| **LIKWID**   | `--enable-debug` (4 MPI ranks × 4 OMP threads) | 306.15      | 1000 humans, 100 steps  |
| **LIKWID**   | `--enable-ipo` (16 OMP threads)                | 1025.74     | 1000 humans, 100 steps  |

#### MAQAO Analysis Results

- **Total Time**: 294.31 seconds
- **Profiled Time**: 239.35 seconds (81.3% coverage)
- **Configuration**: 4 MPI ranks with 4 OpenMP threads each

#### Key Differences in Methodology

| Aspect              | LIKWID                           | MAQAO                            |
|---------------------|----------------------------------|----------------------------------|
| **Scope**           | Hardware counter-based profiling | Static + dynamic analysis        |
| **Granularity**     | Thread/rank level                | Loop level                       |
| **Memory Analysis** | Aggregate bandwidth measurements | Detailed access pattern analysis |
| **Vectorization**   | Overall percentage               | Per-loop analysis                |
| **Output**          | Performance counters             | Interactive HTML reports         |

#### Performance Insights

1. **Runtime Efficiency**
  - **MAQAO MPI run (294s)** vs **LIKWID debug (306s)**: Similar performance
2. **Analysis Depth**
  - **LIKWID** provides quantitative metrics (1515 MB/s memory BW, ~77 MFLOP/s)
  - **MAQAO** provides qualitative loop-by-loop optimization guidance
3. **Profiling Coverage**
  - **MAQAO** achieved 81.3% profiling coverage
  - **LIKWID** measures specific hardware events with 100% coverage
4. **Memory Access Patterns** - MAQAO shows various loops with different memory access patterns: 
  - Some loops show spans of 652,539-171,205 bytes
  - Memory access patterns vary from simple (6 loads) to complex (18 mixed loads/stores)

#### Complementary Analysis

The tools complement each other:
- **LIKWID** gives you the **"what"** (quantitative performance)
- **MAQAO** gives you the **"where and why"** (code locations and optimization opportunities)

---

## Intra (shared memory, threads)

## Inter (MPI)

### Strong scaling

### Weak scaling

### Load balancing
