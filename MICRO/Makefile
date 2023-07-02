vec:
	clang -O3 -std=c99 -fopenmp -mcpu=apple-m1 -D_GNU_SOURCE -o micro-vec  micro.c

novec:
	clang -O3 -mcpu=apple-m1 -fno-vectorize -fno-slp-vectorize -std=c99 -fopenmp -D_GNU_SOURCE -o micro-novec  micro.c

vec-bench:
	make vec
	./bench.pl 1 seq vec> micro_bench_vec_seq.dat

novec-bench:
	make novec
	./bench.pl 1 seq vec> micro_bench_vec_seq.dat

plot:
	gnuplot bench.plot

clean:
	rm micro-vec
	rm micro-novec

all:
	make novec-bench
	make novec-bench
	make plot
