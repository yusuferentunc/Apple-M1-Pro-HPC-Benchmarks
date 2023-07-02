NUSIF-2D = NuSiF-Solver/2D-seq
NUSIF-3D = NuSiF-Solver/3D-seq
POISSON = Poisson-Solver/2D-seq
BANDWIDTH = TheBandwidthBenchmark

micro:
	make all -C MICRO

nusif:
	make -C $(NUSIF-2D)
	$(NUSIF-2D)/exe-CLANG $(NUSIF-2D)/dcavity.par > results/M1-NuSiF-2D-dcavity.txt
	$(NUSIF-2D)/exe-CLANG $(NUSIF-2D)/canal.par > results/M1-NuSiF-2D-canal.txt
	make -C $(NUSIF-3D)
	$(NUSIF-3D)/exe-CLANG $(NUSIF-3D)/dcavity.par > results/M1-NuSiF-3D-dcavity.txt
	$(NUSIF-3D)/exe-CLANG $(NUSIF-3D)/canal.par > results/M1-NuSiF-3D-canal.txt

poisson:
	make -C $(POISSON)
	$(POISSON)/exe-CLANG $(POISSON)/poisson.par 1 > results/M1-Poisson-2D-poisson-plain-sor.txt
	$(POISSON)/exe-CLANG $(POISSON)/poisson.par 2 > results/M1-Poisson-2D-poisson-red-black-sor.txt
	$(POISSON)/exe-CLANG $(POISSON)/poisson.par 3 > results/M1-Poisson-2D-poisson-red-black-accelerated-sor.txt
	
bandwidth:
	make -C $(BANDWIDTH)
	$(BANDWIDTH)/bwbench-CLANG > results/M1-TheBandwidthBenchmark.txt


all:
	echo "MICRO Benchmark"
	make micro
	echo "NuSiF Solver Benchmark"
	make nusif 
	echo "Poisson Solver Benchmark"
	make poisson
	echo "The Bandwidth Benchmark"
	make bandwidth