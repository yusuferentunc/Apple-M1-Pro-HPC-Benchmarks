# Apple-M1-Pro-HPC-Benchmarks

- These benchmarks are updated to run in Apple M1 Pro.
- Project is designated for the sequential benchmarking (one CPU core only )
- Project can also be used in other ARM architecture Apple Silicon products. 
Only thing that may be change is `-mcpu=apple-m1` optimization flag at some of the benchmarks.

## Requirements
- Apple Clang
- gnuplot
- fopenmp

## Build and Run
- You can run all the benchmarks with `make all` command in the parent directory
- Results can be found in `results` file
- For the specific benchmark, you can check the details of makefile

## Benchmarks
* [TheBandwidthBenchmark](https://github.com/RRZE-HPC/TheBandwidthBenchmark)
* [MD-Bench](https://github.com/RRZE-HPC/MD-Bench)
* MICRO
* NuSiF-Solver
* Poisson-Solver