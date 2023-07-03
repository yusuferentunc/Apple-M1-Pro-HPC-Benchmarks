all:
	make vec-bench
	make novec-bench
	make plot

vec:
	clang -O3 -std=c99 -fopenmp -mcpu=apple-m1 -D_GNU_SOURCE -o micro-vec  micro.c

novec:
	clang -O3 -mcpu=apple-m1 -fno-vectorize -fno-slp-vectorize -std=c99 -fopenmp -D_GNU_SOURCE -o micro-novec  micro.c

vec-bench:
	make vec
	./bench.pl 1 seq vec> m1-vec-seq.dat

novec-bench:
	make novec
	./bench.pl 1 seq vec> m1-novec-seq.dat

plot:
	gnuplot bench.plot

clean:
	rm micro-vec
	rm micro-novec



