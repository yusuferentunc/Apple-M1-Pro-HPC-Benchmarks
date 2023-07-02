## Task: Use a microbenchmark to explore the memory hierarchy

In this hands-on exercise you will use variants of the Schoenauer Triad microbenchmark to investigate the performance impact of data being located in different levels of the memory hierarchy.
Besides a sequential version you will also look at threaded versions with and without OpenMP worksharing and experience the overhead introduced by OpenMP.

Time to finish: around 60 Minutes.

## Preparation

* Copy the files micro.c, bench.pl, and the bench plot files
* Get an interactive single node job on Emmy cluster: `$ qsub -I -lnodes=1:ppn=40:f2.2 -lwalltime=01:00:00`
* Load Intel compiler module: `$ module load intel64`
* Load Likwid module: `$ module load likwid`

## Investigate the benchmark code

Open micro.c with a text editor. The benchmark contains three variants of
the so called Schoenauer Triad (`a[i]=b[i]+c[i]*d[i]`). This benchmark is
data transfer bound on any known architecture.

* What is the difference between the three variants sequential, throughput and worksharing?
* How does the argument `N` for the vector length correlate with the resident memory size?

## Compile benchmark

* Compile a threaded OpenMP binary with optimizing flags:
`$ icc -Ofast -xHost -std=c99 -qopenmp -D_GNU_SOURCE -o micro-vec  micro.c`

* Compile a threaded OpenMP binary without SIMD vectorization:
`$ icc -O3 -no-vec -std=c99 -qopenmp -D_GNU_SOURCE -o micro-novec  micro.c`

Copy the desired variant to micro before using the bench.pl script.

## Run benchmark

* Test if the benchmark is working:
    * `$ ./micro`
    * `$ ./micro 0 5000`

* Use the helper script ./bench.pl to scan data set size.

Use it as follows:
`$ ./bench.pl <numcores> <seq|tp|ws>`

You can generate a png plot of the result using gnuplot with:
`$ gnuplot bench.plot`

The `bench.plot` configuration expects the output data in `bench.dat` (optimized) and `bench-novec.dat` (no SIMD vectorization)!

* Explore **sequential** performance across memory hierarchy:
On compute node:
1. Copy micro-vec to micro
2. Execute `$ ./bench.pl 1 seq > bench.dat`
3. Copy micro-novec to micro
2. Execute `$ ./bench.pl 1 seq > bench-novec.dat`
On Emmy frontend: `$ gnuplot bench.plot`
The result image is in micro.png.

Try to reconstruct the cache sizes of the different levels from the plot. How does this compare
to the values reported by likwid-topology?
Is the vectorized L1 performance the theoretical limit?
What is the difference in performance between the SIMD and non-SIMD version?
Is this to be expected?
Try to derive an upper performance bound for the L1 Schoenauer Triad (SIMD and non-SIMD).
Does the measurement meet your model?

**Optional** if used without pragma vector aligned in the beginning:
Have a look at the assembly code to check if your model assumptions were correct:
`$ icc -S -Ofast -xHost -std=c99 -qopenmp -D_GNU_SOURCE -o micro.s  micro.c`
How can you fix this?

* Explore parallel **throughput** performance across the memory hierarchy:
On compute node:
1. Copy micro-vec to micro
2. Execute `$ for i in 1 2 4 6 8 10 20; do ./bench.pl $i tp > bench-tp-$i.dat; done`
On Emmy frontend: `$ gnuplot bench-tp.plot`
All core count variants have to be available with naming scheme `bench-tp-<i>.dat`!
The result image is in micro-tp.png.
Does the performance scale in all hierarchy levels?
What aspect does this benchmark test?

* Explore parallel **work sharing** performance across the memory hierarchy:
On compute node:
1. Copy micro-vec to micro
2. Execute `$ for i in 1 10 20; do ./bench.pl $i ws > bench-ws-$i.dat; done`
On Emmy frontend: `$ gnuplot bench-tp.plot`
All core count variants have to be available with naming scheme `bench-ws-<i>.dat`, in addition you also need the sequential result in `bench.plot`!
The result image is in micro-ws.png.
What is different from the throughput case?
Can you explain the result?
